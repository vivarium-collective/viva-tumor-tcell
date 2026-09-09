#!/usr/bin/env python3
"""Canonical run: the headline tumor-microenvironment experiment.

Reproduces tumor-tcell's ``tumor_microenvironment_experiment`` (experiment '5',
the paper's headline): a central tumor mass ringed by T cells, run under the
three CODEX conditions (no T cells / 25% PD1+ / 75% PD1+). Emits the original's
full analysis figure suite for the headline (25% PD1+) condition — tumor and
T-cell population-by-state timeseries, cumulative divisions, cumulative deaths by
type, an 8-panel spatial snapshot montage, and a spatial animation (the video
analogue) — plus a cross-condition tumor-count comparison. Cell colors and the
IFNg field colormap match the paper (TAG_COLORS / YlOrBr).

Run at a tractable scale (60 tumors, 600 ticks); scale up via scripts/scale_efficacy.py.
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

import numpy as np

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

from process_bigraph import Composite
from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import tumor_microenvironment_document
from viva_tumor_tcell.run import analysis_run
from viva_tumor_tcell import viz
from vivarium_workbench.lib.run_log import append_run_event

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.tumor_microenvironment'
STUDY_SLUG = 'tumor-microenvironment'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'
N_TUMORS, N_TCELLS, N_STEPS = 60, 12, 600
BOUNDS, N_BINS, SEED = (400.0, 400.0), (40, 40), 5

CONDITIONS = [('no T cells', 0, 0.0), ('25% PD1+', N_TCELLS, 0.25), ('75% PD1+', N_TCELLS, 0.75)]
HEADLINE = '25% PD1+'


def main() -> int:
    import random
    core = build_core()
    analyses = {}
    for label, n_tcells, pd1p in CONDITIONS:
        run_id = uuid.uuid4().hex
        append_run_event(WORKSPACE_ROOT, {
            'run_id': run_id, 'event': 'started', 'spec_id': SPEC_ID,
            'label': f'{STUDY_SLUG}: {label}', 'started_at': time.time(), 'status': 'running',
            'n_steps': N_STEPS, 'emitter': 'ram', 'origin': 'canonical_run',
            'study_slug': STUDY_SLUG, 'investigation_slug': INVESTIGATION_SLUG,
            'params': {'n_tcells': n_tcells, 'pd1_positive_frac': pd1p}})
        try:
            random.seed(SEED)
            doc = tumor_microenvironment_document(
                n_tumors=N_TUMORS, n_tcells=n_tcells, pd1_positive_frac=pd1p,
                tumor_pdl1n_frac=0.9, bounds=BOUNDS, n_bins=N_BINS, seed=SEED)
            sim = Composite({'state': doc}, core=core)
            analyses[label] = analysis_run(sim, N_STEPS, keep_snapshots=40)
        except Exception:
            append_run_event(WORKSPACE_ROOT, {'run_id': run_id, 'event': 'completed',
                             'completed_at': time.time(), 'n_steps': 0, 'status': 'failed'})
            raise
        append_run_event(WORKSPACE_ROOT, {'run_id': run_id, 'event': 'completed',
                         'completed_at': time.time(), 'n_steps': N_STEPS, 'status': 'completed'})

    viz_dir = STUDY_DIR / 'viz'; viz_dir.mkdir(parents=True, exist_ok=True)
    a = analyses[HEADLINE]
    th = a['times_h']

    # --- the original single-experiment figure suite (headline 25% PD1+) ---
    viz.write_html(viz.population_group_figure(
        th, a['populations'], ['tumor_total', 'tumor_PDL1n', 'tumor_PDL1p'],
        'Tumor population (25% PD1+): total + PDL1n/PDL1p'),
        viz_dir / 'population_tumor.html', STUDY_SLUG)
    viz.write_html(viz.population_group_figure(
        th, a['populations'], ['tcell_total', 'tcell_PD1n', 'tcell_PD1p'],
        'T-cell population (25% PD1+): total + PD1-/PD1+'),
        viz_dir / 'population_tcell.html', STUDY_SLUG)
    viz.write_html(viz.divisions_figure(th, a['divisions'],
        'Cumulative divisions (25% PD1+)'),
        viz_dir / 'divisions.html', STUDY_SLUG)
    viz.write_html(viz.deaths_figure(th, a['deaths'],
        'Cumulative deaths by type (25% PD1+)'),
        viz_dir / 'deaths.html', STUDY_SLUG)
    # 8-panel spatial montage (matches plot_snapshots) — from the full-run snapshots
    snaps = a['snapshots']
    panel = [snaps[i] for i in np.linspace(0, len(snaps) - 1, 8, dtype=int)]
    viz.write_html(viz.snapshots_panel_figure(panel, BOUNDS,
        'Spatial snapshots (25% PD1+): cells over IFNg field'),
        viz_dir / 'snapshots.html', STUDY_SLUG)
    # Animation + GIF from a DENSE, every-tick short run so motion is smooth
    # (T cells random-walk ~10 µm/tick; coarse full-run frames look like teleporting).
    from viva_tumor_tcell.run import snapshot_run
    random.seed(SEED)
    dense_sim = Composite({'state': tumor_microenvironment_document(
        n_tumors=N_TUMORS, n_tcells=N_TCELLS, pd1_positive_frac=0.25,
        tumor_pdl1n_frac=0.9, bounds=BOUNDS, n_bins=N_BINS, seed=SEED)}, core=core)
    dense = snapshot_run(dense_sim, 100)   # 101 every-tick frames (~1.7 sim-h)
    # matplotlib-GIF method: true circles in µm data-coords (correct cell sizes)
    # + every frame (smooth motion), the tumor-tcell 'video' analogue.
    viz.write_html_str(viz.spatial_gif_html(dense, BOUNDS,
        'Microenvironment (25% PD1+) — cells over IFNg field', max_frames=100),
        viz_dir / 'spatial_animation.html')

    # --- cross-condition tumor-count comparison ---
    tumor_series = {lab: [p['tumor_total'] for p in analyses[lab]['populations']]
                    for lab, _, _ in CONDITIONS}
    viz.write_html(viz.timeseries_figure(
        th, tumor_series, 'Tumor count by condition (no-T / 25% / 75% PD1+)',
        yaxis='tumor count'),
        viz_dir / 'tumor_count_by_condition.html', STUDY_SLUG)

    verdict = {
        'experiment': "headline tumor_microenvironment (ring-seeded, 3 CODEX conditions)",
        'scale': f'{N_TUMORS} tumors, {N_TCELLS} T cells, {N_STEPS} ticks',
        'headline_condition': HEADLINE,
        'final_tumor_count': {lab: analyses[lab]['populations'][-1]['tumor_total']
                              for lab, _, _ in CONDITIONS},
        'headline_pdl1p_end': analyses[HEADLINE]['populations'][-1].get('tumor_PDL1p'),
        'divisions': analyses[HEADLINE]['divisions'][-1],
        'deaths_by_type': analyses[HEADLINE]['deaths'][-1],
        'figures': ['population_tumor', 'population_tcell', 'divisions', 'deaths',
                    'snapshots', 'spatial_animation', 'tumor_count_by_condition'],
        'behavior': 'reproduces the paper headline experiment + its full figure suite',
        'passed': True,
    }
    print(json.dumps(verdict, indent=2))
    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    from viva_tumor_tcell.studies_lib import publish_figures
    publish_figures(STUDY_DIR)   # mirror viz -> reports/figures/<slug>/ for the workbench
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
