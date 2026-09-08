#!/usr/bin/env python3
"""Canonical run: IFNg-driven PDL1n->PDL1p phenotype conversion.

Compares the tumor mass with active T cells (25% PD1+) against a matched no-T
control. Reports the PDL1p tumor fraction and the IFNg field over time; the
conversion and cytokine field should be elevated only when active T cells are
present.
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
from viva_tumor_tcell.run import run_condition, pdl1p_fraction, ifng_max
from viva_tumor_tcell import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.tumor_microenvironment'
STUDY_SLUG = 'phenotype-conversion'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'

N_TUMORS = 40
N_STEPS = 600
BOUNDS = (300.0, 300.0)
N_BINS = (30, 30)
SEED = 5

CONDITIONS = [
    ('with active T cells', dict(n_tcells=12, pd1_positive_frac=0.25)),
    ('no T cells',          dict(n_tcells=0,  pd1_positive_frac=0.0)),
]


def main() -> int:
    import random
    core = build_core()
    pdl1p_series = {}
    ifng_series = {}
    frames_withT = None

    for label, kw in CONDITIONS:
        run_id = uuid.uuid4().hex
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'started', 'spec_id': SPEC_ID,
            'label': f'{STUDY_SLUG}: {label}', 'started_at': time.time(), 'status': 'running',
            'n_steps': N_STEPS, 'emitter': 'ram', 'origin': 'canonical_run',
            'study_slug': STUDY_SLUG, 'investigation_slug': INVESTIGATION_SLUG, 'params': kw})
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
        pdl1p_series[label] = [pdl1p_fraction(f) * 100 for f in frames]
        ifng_series[label] = [ifng_max(f) for f in frames]
        if label == 'with active T cells':
            frames_withT = frames
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
            'n_steps': N_STEPS, 'status': 'completed'})

    times_h = [i * 60.0 / 3600.0 for i in range(N_STEPS + 1)]

    conv_withT = pdl1p_series['with active T cells'][-1]
    conv_noT = pdl1p_series['no T cells'][-1]
    ifng_withT = max(ifng_series['with active T cells'])
    ifng_noT = max(ifng_series['no T cells'])
    passed = conv_withT > conv_noT and ifng_withT > 0 and ifng_noT == 0.0
    verdict = {
        'final_pdl1p_pct': {'with T': round(conv_withT, 1), 'no T': round(conv_noT, 1)},
        'peak_ifng_ngml': {'with T': round(ifng_withT, 2), 'no T': round(ifng_noT, 2)},
        'behavior': 'PDL1p conversion and IFNg elevated with active T cells vs. none',
        'passed': bool(passed),
    }
    print(json.dumps(verdict, indent=2))

    viz_dir = STUDY_DIR / 'viz'
    viz_dir.mkdir(parents=True, exist_ok=True)
    fig1 = viz.timeseries_figure(
        times_h, pdl1p_series,
        title='IFNg-driven conversion: PDL1p tumor fraction over time',
        yaxis='PDL1p tumors (%)')
    viz.write_html(fig1, viz_dir / 'pdl1p_fraction.html', STUDY_SLUG)
    if frames_withT is not None:
        fig2 = viz.spatial_animation_figure(
            frames_withT, BOUNDS,
            title='Active T cells — IFNg field builds, tumors convert to PDL1p')
        viz.write_html(fig2, viz_dir / 'spatial_conversion.html', STUDY_SLUG)

    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
