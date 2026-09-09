#!/usr/bin/env python3
"""Canonical run: IFNg secretion and PDL1n->PDL1p conversion (6 replicate seeds).

Active T cells reliably build a measurable IFNg field (the driver of tumor
phenotype conversion) — robust across seeds. The downstream PDL1p conversion
*fraction* is directionally higher with T cells but does not separate from
growth noise at this tractable scale; both are reported honestly (mean ± std).
"""
from __future__ import annotations

import json
from pathlib import Path

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

import numpy as np

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import tumor_microenvironment_document
from viva_tumor_tcell.run import replicate_series, mean_std, pdl1p_fraction, ifng_max
from viva_tumor_tcell import viz
from viva_tumor_tcell.studies_lib import SEEDS, record_run

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.tumor_microenvironment'
STUDY_SLUG = 'phenotype-conversion'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'
N_TUMORS, N_STEPS = 40, 600
BOUNDS, N_BINS = (300.0, 300.0), (30, 30)

CONDITIONS = [('with active T cells', 12, 0.25), ('no T cells', 0, 0.0)]
OBS = {'ifng': ifng_max, 'pdl1p': lambda f: pdl1p_fraction(f) * 100}


def main() -> int:
    core = build_core()
    ifng_stats, pdl1p_stats = {}, {}
    per_seed_ifng = {}
    frames_withT = None
    for label, n_tcells, pd1p in CONDITIONS:
        def _run(n_tcells=n_tcells, pd1p=pd1p):
            return replicate_series(
                core,
                lambda s: tumor_microenvironment_document(
                    n_tumors=N_TUMORS, n_tcells=n_tcells, pd1_positive_frac=pd1p,
                    tumor_pdl1n_frac=0.9, bounds=BOUNDS, n_bins=N_BINS, seed=s),
                N_STEPS, SEEDS, OBS)
        arrays, frames = record_run(WORKSPACE_ROOT, STUDY_SLUG, INVESTIGATION_SLUG,
                                    SPEC_ID, label, {'n_tcells': n_tcells}, N_STEPS, _run)
        ifng_stats[label] = mean_std(arrays['ifng'])
        pdl1p_stats[label] = mean_std(arrays['pdl1p'])
        per_seed_ifng[label] = arrays['ifng']
        if n_tcells:
            frames_withT = frames

    times_h = [i * 60.0 / 3600.0 for i in range(N_STEPS + 1)]
    # PASS metric: IFNg field forms with active T cells in every seed, and stays 0 without.
    peak_withT = np.max(per_seed_ifng['with active T cells'], axis=1)   # per seed
    peak_noT = np.max(per_seed_ifng['no T cells'], axis=1)
    passed = bool(np.all(peak_withT > 0) and np.all(peak_noT == 0))
    conv_withT = pdl1p_stats['with active T cells'][0][-1]
    conv_noT = pdl1p_stats['no T cells'][0][-1]
    verdict = {
        'primary_metric': 'peak IFNg field with active T cells > 0 in every seed; 0 without',
        'ifng_positive_seeds': f'{int(np.sum(peak_withT > 0))}/{len(SEEDS)} with T; '
                               f'{int(np.sum(peak_noT > 0))}/{len(SEEDS)} without',
        'secondary_pdl1p_pct_end': {'with T': round(float(conv_withT), 1),
                                    'no T': round(float(conv_noT), 1)},
        'behavior': 'active T cells reliably secrete IFNg; PDL1p conversion is directional but scale-limited',
        'passed': passed,
    }
    print(json.dumps(verdict, indent=2))

    viz_dir = STUDY_DIR / 'viz'; viz_dir.mkdir(parents=True, exist_ok=True)
    fig1 = viz.timeseries_band_figure(
        times_h, ifng_stats,
        title='IFNg field secreted by active T cells (mean ± std, 6 seeds)',
        yaxis='peak IFNg (ng/mL)')
    viz.write_html(fig1, viz_dir / 'ifng_field.html', STUDY_SLUG)
    fig2 = viz.timeseries_band_figure(
        times_h, pdl1p_stats,
        title='PDL1p tumor fraction — directional, scale-limited (mean ± std, 6 seeds)',
        yaxis='PDL1p tumors (%)')
    viz.write_html(fig2, viz_dir / 'pdl1p_fraction.html', STUDY_SLUG)
    if frames_withT is not None:
        fig3 = viz.spatial_animation_figure(
            frames_withT, BOUNDS, title='Active T cells build a local IFNg field')
        viz.write_html(fig3, viz_dir / 'spatial_conversion.html', STUDY_SLUG)

    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    from viva_tumor_tcell.studies_lib import publish_figures
    publish_figures(STUDY_DIR)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
