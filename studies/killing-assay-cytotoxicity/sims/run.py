#!/usr/bin/env python3
"""Canonical run: in-vitro-style killing / cytotoxicity assay (6 replicate seeds).

Well-mixed tumor cells (50% PDL1+) + active T cells (25% PD1+, 1.5:1 ratio) vs. a
matched no-T control. Cytotoxicity = (control - experiment)/control x 100 on the
final tumor count, paired per seed. Reports the mean ± std across seeds.
"""
from __future__ import annotations

import json
from pathlib import Path

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

import numpy as np

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import killing_assay_document
from viva_tumor_tcell.run import replicate_series, mean_std, tumor_count, analysis_run
from viva_tumor_tcell.composites.microenvironment import killing_assay_document as _kad
from viva_tumor_tcell import viz
from viva_tumor_tcell.studies_lib import SEEDS, record_run

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.killing_assay'
STUDY_SLUG = 'killing-assay-cytotoxicity'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'
N_TUMORS, N_STEPS = 30, 700
BOUNDS, N_BINS = (180.0, 180.0), (18, 18)

CONDITIONS = [('+ T cells (50% PDL1+)', True), ('control (no T)', False)]


def main() -> int:
    core = build_core()
    stats = {}
    counts = {}
    frames_exp = None
    for label, inc in CONDITIONS:
        def _run(inc=inc):
            return replicate_series(
                core,
                lambda s: killing_assay_document(
                    n_tumors=N_TUMORS, include_tcells=inc, pdl1_positive_frac=0.5,
                    tumor_t_ratio=1.5, pd1_positive_frac=0.25,
                    bounds=BOUNDS, n_bins=N_BINS, seed=s),
                N_STEPS, SEEDS, tumor_count)
        arr, frames = record_run(WORKSPACE_ROOT, STUDY_SLUG, INVESTIGATION_SLUG,
                                 SPEC_ID, label, {'include_tcells': inc}, N_STEPS, _run)
        stats[label] = mean_std(arr)
        counts[label] = arr
        if inc:
            frames_exp = frames

    times_h = [i * 60.0 / 3600.0 for i in range(N_STEPS + 1)]
    exp_final = counts['+ T cells (50% PDL1+)'][:, -1]     # per seed
    ctl_final = counts['control (no T)'][:, -1]
    cytotox = np.where(ctl_final > 0, (ctl_final - exp_final) / ctl_final * 100, 0.0)
    passed = bool(np.mean(cytotox) > 0)
    verdict = {
        'metric': 'cytotoxicity % = (control - experiment)/control, per seed',
        'cytotoxicity_pct': f'{np.mean(cytotox):.1f} ± {np.std(cytotox):.1f}',
        'positive_seeds': f'{int(np.sum(cytotox > 0))}/{len(SEEDS)}',
        'behavior': 'mean cytotoxicity vs. matched no-T control is positive',
        'passed': passed,
    }
    print(json.dumps(verdict, indent=2))

    viz_dir = STUDY_DIR / 'viz'; viz_dir.mkdir(parents=True, exist_ok=True)
    fig1 = viz.timeseries_band_figure(
        times_h, stats,
        title=f'Killing assay: tumor count (mean ± std, 6 seeds) — '
              f'cytotoxicity {np.mean(cytotox):.0f} ± {np.std(cytotox):.0f}%',
        yaxis='tumor count')
    viz.write_html(fig1, viz_dir / 'tumor_count_vs_control.html', STUDY_SLUG)
    if frames_exp is not None:
        fig2 = viz.spatial_animation_figure(
            frames_exp, BOUNDS, title='Killing assay (+T) — well-mixed cytotoxicity')
        viz.write_html(fig2, viz_dir / 'spatial_killing.html', STUDY_SLUG)

    # --- matching analysis figures (population + deaths by type), mirroring the
    #     original killing_experiment's death_group_plot / population_group_plot,
    #     from one representative +T run ---
    import random
    from process_bigraph import Composite
    random.seed(SEEDS[0])
    a = analysis_run(
        Composite({'state': _kad(n_tumors=N_TUMORS, include_tcells=True, pdl1_positive_frac=0.5,
                                 tumor_t_ratio=1.5, pd1_positive_frac=0.25,
                                 bounds=BOUNDS, n_bins=N_BINS, seed=SEEDS[0])}, core=core),
        N_STEPS, keep_snapshots=6)
    ath = a['times_h']
    viz.write_html(viz.population_group_figure(
        ath, a['populations'], ['tumor_total', 'tumor_PDL1n', 'tumor_PDL1p'],
        'Killing assay — tumor population by state (+T)'),
        viz_dir / 'population_tumor.html', STUDY_SLUG)
    viz.write_html(viz.deaths_figure(ath, a['deaths'],
        'Killing assay — cumulative tumor/T-cell deaths by type (+T)'),
        viz_dir / 'deaths.html', STUDY_SLUG)

    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
