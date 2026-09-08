#!/usr/bin/env python3
"""Larger-scale efficacy probe: does population tumor-count suppression (and the
25%<75% PD1+ ordering) emerge with more cells + more time than the tractable
studies?

Memory-flat: records only scalar observables per tick (never retains frames), so
it can run long without the per-frame snapshot memory blowup. Saves partial
results after every run so a kill still leaves data. Intended to run in the
background and be sized up toward paper scale.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from process_bigraph import Composite
from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import tumor_microenvironment_document
from viva_tumor_tcell.run import tumor_count, pdl1p_fraction, pd1p_fraction, ifng_max

# --- scale knobs (override via argv: n_tumors n_tcells steps bounds nbins nseeds) ---
N_TUMORS = int(sys.argv[1]) if len(sys.argv) > 1 else 150
N_TCELLS = int(sys.argv[2]) if len(sys.argv) > 2 else 40
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
BOUNDS = float(sys.argv[4]) if len(sys.argv) > 4 else 600.0
NBINS = int(sys.argv[5]) if len(sys.argv) > 5 else 60
NSEEDS = int(sys.argv[6]) if len(sys.argv) > 6 else 3
SEEDS = list(range(1, NSEEDS + 1))
INTERVAL = 60.0

CONDITIONS = [('no T cells', 0, 0.0), ('25% PD1+', N_TCELLS, 0.25), ('75% PD1+', N_TCELLS, 0.75)]
OUT = REPO / 'reports' / 'scale_efficacy_results.json'
OUT.parent.mkdir(parents=True, exist_ok=True)


def run_one(core, n_tcells, pd1p, seed):
    """Run one condition/seed, recording scalar observables per tick (flat memory)."""
    import random
    random.seed(seed)
    doc = tumor_microenvironment_document(
        n_tumors=N_TUMORS, n_tcells=n_tcells, pd1_positive_frac=pd1p,
        tumor_pdl1n_frac=0.9, bounds=(BOUNDS, BOUNDS), n_bins=(NBINS, NBINS), seed=seed)
    sim = Composite({'state': doc}, core=core)
    from viva_tumor_tcell.run import snapshot
    rec = {'tumor': [], 'pdl1p': [], 'pd1p': [], 'ifng': []}
    for step in range(STEPS + 1):
        if step:
            sim.run(INTERVAL)
        f = snapshot(sim)
        rec['tumor'].append(tumor_count(f))
        rec['pdl1p'].append(pdl1p_fraction(f) * 100)
        rec['pd1p'].append(pd1p_fraction(f) * 100)
        rec['ifng'].append(ifng_max(f))
        del f
    return rec


def main():
    core = build_core()
    results = {'config': dict(n_tumors=N_TUMORS, n_tcells=N_TCELLS, steps=STEPS,
                              bounds=BOUNDS, nbins=NBINS, seeds=SEEDS),
               'runs': {}}
    t0 = time.time()
    for label, n_tcells, pd1p in CONDITIONS:
        results['runs'][label] = {}
        for seed in SEEDS:
            t = time.time()
            rec = run_one(core, n_tcells, pd1p, seed)
            dt = time.time() - t
            results['runs'][label][str(seed)] = {
                'final_tumor': rec['tumor'][-1], 'final_pdl1p': rec['pdl1p'][-1],
                'final_pd1p': rec['pd1p'][-1], 'peak_ifng': max(rec['ifng']),
                'tumor_series': rec['tumor'], 'seconds': round(dt, 1)}
            OUT.write_text(json.dumps(results, indent=2))
            print(f"[{time.time()-t0:6.0f}s] {label:9s} seed{seed}: "
                  f"final_tumor={rec['tumor'][-1]:4d} PDL1p={rec['pdl1p'][-1]:.0f}% "
                  f"PD1p={rec['pd1p'][-1]:.0f}% IFNg={max(rec['ifng']):.1f} ({dt:.0f}s)",
                  flush=True)

    # --- analysis ---
    def finals(label):
        return np.array([results['runs'][label][str(s)]['final_tumor'] for s in SEEDS], float)
    noT, p25, p75 = finals('no T cells'), finals('25% PD1+'), finals('75% PD1+')
    s25 = (noT - p25) / noT * 100
    s75 = (noT - p75) / noT * 100
    summary = {
        'final_tumor_mean_std': {k: f'{np.mean(finals(k)):.1f} ± {np.std(finals(k)):.1f}'
                                 for k in ('no T cells', '25% PD1+', '75% PD1+')},
        'suppression_25_pct': f'{np.mean(s25):.1f} ± {np.std(s25):.1f}',
        'suppression_75_pct': f'{np.mean(s75):.1f} ± {np.std(s75):.1f}',
        'ordering_25_gt_75_seeds': f'{int(np.sum(s25 >= s75))}/{len(SEEDS)}',
        'suppression_25_positive_seeds': f'{int(np.sum(s25 > 0))}/{len(SEEDS)}',
    }
    results['summary'] = summary
    OUT.write_text(json.dumps(results, indent=2))
    print('=== SCALE EFFICACY SUMMARY ===')
    print(json.dumps(summary, indent=2), flush=True)
    print(f'total {time.time()-t0:.0f}s -> {OUT}', flush=True)


if __name__ == '__main__':
    main()
