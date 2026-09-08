#!/usr/bin/env python3
"""Canonical run: T-cell efficacy vs. PD1+ (exhausted) fraction.

Three CODEX-style conditions over a central tumor mass ringed by T cells:
no-T control, 25% PD1+, 75% PD1+. Reports tumor count over time, renders an
interactive time-series + a spatial animation, and records the runs.
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import tumor_microenvironment_document
from viva_tumor_tcell.run import run_condition, tumor_count
from viva_tumor_tcell import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.tumor_microenvironment'
STUDY_SLUG = 'tcell-efficacy-pd1'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'

N_TUMORS = 40
N_TCELLS = 12
N_STEPS = 600
BOUNDS = (300.0, 300.0)
N_BINS = (30, 30)
SEED = 5

CONDITIONS = [
    ('no T cells',  dict(n_tcells=0,        pd1_positive_frac=0.0)),
    ('25% PD1+',    dict(n_tcells=N_TCELLS, pd1_positive_frac=0.25)),
    ('75% PD1+',    dict(n_tcells=N_TCELLS, pd1_positive_frac=0.75)),
]


def main() -> int:
    import random
    core = build_core()
    results = {}
    frames_25 = None

    for label, kw in CONDITIONS:
        run_id = uuid.uuid4().hex
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'started', 'spec_id': SPEC_ID,
            'label': f'{STUDY_SLUG}: {label}', 'started_at': time.time(), 'status': 'running',
            'n_steps': N_STEPS, 'emitter': 'ram', 'origin': 'canonical_run',
            'study_slug': STUDY_SLUG, 'investigation_slug': INVESTIGATION_SLUG, 'params': kw,
        })
        try:
            random.seed(SEED)
            doc = tumor_microenvironment_document(
                n_tumors=N_TUMORS, bounds=BOUNDS, n_bins=N_BINS, tumor_pdl1n_frac=0.9,
                seed=SEED, **kw)
            frames = run_condition(core, doc, N_STEPS)
        except Exception:
            append_run_event(WORKSPACE_ROOT, {
                'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
                'n_steps': 0, 'status': 'failed'})
            raise
        results[label] = [tumor_count(f) for f in frames]
        if label == '25% PD1+':
            frames_25 = frames
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
            'n_steps': N_STEPS, 'status': 'completed'})

    times_h = [i * 60.0 / 3600.0 for i in range(N_STEPS + 1)]

    # --- verdict (directional; single seed) ---
    finals = {k: v[-1] for k, v in results.items()}
    no_t = finals['no T cells']
    supp_25 = (no_t - finals['25% PD1+']) / no_t if no_t else 0.0
    supp_75 = (no_t - finals['75% PD1+']) / no_t if no_t else 0.0
    passed = supp_25 > 0 and supp_75 >= 0 and supp_25 >= supp_75
    verdict = {
        'final_tumor_count': finals,
        'suppression_25pct': round(supp_25 * 100, 1),
        'suppression_75pct': round(supp_75 * 100, 1),
        'behavior': 'T cells suppress tumor growth; 25% PD1+ suppresses >= 75% PD1+',
        'passed': bool(passed),
    }
    print(json.dumps(verdict, indent=2))

    # --- visualizations ---
    viz_dir = STUDY_DIR / 'viz'
    viz_dir.mkdir(parents=True, exist_ok=True)
    fig1 = viz.timeseries_figure(
        times_h, results,
        title='T-cell efficacy: tumor count vs. PD1+ fraction',
        yaxis='tumor count')
    viz.write_html(fig1, viz_dir / 'tumor_count_by_condition.html', STUDY_SLUG)
    if frames_25 is not None:
        fig2 = viz.spatial_animation_figure(
            frames_25, BOUNDS,
            title='25% PD1+ condition — cells over the IFNg field')
        viz.write_html(fig2, viz_dir / 'spatial_25pct.html', STUDY_SLUG)

    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
