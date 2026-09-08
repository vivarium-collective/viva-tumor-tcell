#!/usr/bin/env python3
"""Canonical run: in-vitro-style killing / cytotoxicity assay.

Well-mixed tumor cells (50% PDL1+) with active T cells (25% PD1+) vs. a matched
no-T control. Cytotoxicity = (control - experiment) / control x 100, computed on
the final tumor count, mirroring the paper's assay readout.
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import killing_assay_document
from viva_tumor_tcell.run import run_condition, tumor_count
from viva_tumor_tcell import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.killing_assay'
STUDY_SLUG = 'killing-assay-cytotoxicity'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'

N_TUMORS = 30
N_STEPS = 700
BOUNDS = (180.0, 180.0)
N_BINS = (18, 18)
SEED = 5

CONDITIONS = [
    ('+ T cells (50% PDL1+)', dict(include_tcells=True,  pdl1_positive_frac=0.5,
                                   tumor_t_ratio=1.5, pd1_positive_frac=0.25)),
    ('control (no T)',        dict(include_tcells=False, pdl1_positive_frac=0.5)),
]


def main() -> int:
    import random
    core = build_core()
    counts = {}
    frames_exp = None

    for label, kw in CONDITIONS:
        run_id = uuid.uuid4().hex
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'started', 'spec_id': SPEC_ID,
            'label': f'{STUDY_SLUG}: {label}', 'started_at': time.time(), 'status': 'running',
            'n_steps': N_STEPS, 'emitter': 'ram', 'origin': 'canonical_run',
            'study_slug': STUDY_SLUG, 'investigation_slug': INVESTIGATION_SLUG, 'params': kw})
        try:
            random.seed(SEED)
            doc = killing_assay_document(
                n_tumors=N_TUMORS, bounds=BOUNDS, n_bins=N_BINS, seed=SEED, **kw)
            frames = run_condition(core, doc, N_STEPS)
        except Exception:
            append_run_event(WORKSPACE_ROOT, {
                'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
                'n_steps': 0, 'status': 'failed'})
            raise
        counts[label] = [tumor_count(f) for f in frames]
        if kw.get('include_tcells'):
            frames_exp = frames
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
            'n_steps': N_STEPS, 'status': 'completed'})

    times_h = [i * 60.0 / 3600.0 for i in range(N_STEPS + 1)]
    exp_final = counts['+ T cells (50% PDL1+)'][-1]
    ctl_final = counts['control (no T)'][-1]
    cytotox = (ctl_final - exp_final) / ctl_final * 100 if ctl_final else 0.0
    verdict = {
        'final_tumor_count': {'+T': exp_final, 'control': ctl_final},
        'cytotoxicity_pct': round(cytotox, 1),
        'behavior': 'net tumor reduction vs. matched no-T control is positive',
        'passed': bool(cytotox > 0),
    }
    print(json.dumps(verdict, indent=2))

    viz_dir = STUDY_DIR / 'viz'
    viz_dir.mkdir(parents=True, exist_ok=True)
    fig1 = viz.timeseries_figure(
        times_h, counts,
        title=f'Killing assay: tumor count (+T vs control) — cytotoxicity {cytotox:.0f}%',
        yaxis='tumor count')
    viz.write_html(fig1, viz_dir / 'tumor_count_vs_control.html', STUDY_SLUG)
    if frames_exp is not None:
        fig2 = viz.spatial_animation_figure(
            frames_exp, BOUNDS, title='Killing assay (+T) — well-mixed cytotoxicity')
        viz.write_html(fig2, viz_dir / 'spatial_killing.html', STUDY_SLUG)

    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
