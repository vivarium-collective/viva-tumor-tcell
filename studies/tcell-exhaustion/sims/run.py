#!/usr/bin/env python3
"""Canonical run: T-cell exhaustion vs. starting PD1+ fraction (6 replicate seeds).

Seeding more exhausted (PD1+) T cells yields a higher exhausted fraction over
time — a robust, T-cell-intrinsic dynamic. Reports the PD1+ fraction over time
(mean ± std across seeds) for 25%-start vs 75%-start conditions.
"""
from __future__ import annotations

import json
from pathlib import Path

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

import numpy as np

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import tumor_microenvironment_document
from viva_tumor_tcell.run import replicate_series, mean_std, pd1p_fraction
from viva_tumor_tcell import viz
from viva_tumor_tcell.studies_lib import SEEDS, record_run

SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.tumor_microenvironment'
STUDY_SLUG = 'tcell-exhaustion'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'
N_TUMORS, N_TCELLS, N_STEPS = 40, 12, 600
BOUNDS, N_BINS = (300.0, 300.0), (30, 30)

CONDITIONS = [('25% PD1+ start', 0.25), ('75% PD1+ start', 0.75)]


def main() -> int:
    core = build_core()
    stats = {}
    first_frames = {}
    for label, pd1p in CONDITIONS:
        def _run(pd1p=pd1p):
            return replicate_series(
                core,
                lambda s: tumor_microenvironment_document(
                    n_tumors=N_TUMORS, n_tcells=N_TCELLS, pd1_positive_frac=pd1p,
                    tumor_pdl1n_frac=0.9, bounds=BOUNDS, n_bins=N_BINS, seed=s),
                N_STEPS, SEEDS, lambda f: pd1p_fraction(f) * 100)
        arr, frames = record_run(WORKSPACE_ROOT, STUDY_SLUG, INVESTIGATION_SLUG,
                                 SPEC_ID, label, {'pd1_positive_frac': pd1p}, N_STEPS, _run)
        stats[label] = mean_std(arr)
        first_frames[label] = frames

    times_h = [i * 60.0 / 3600.0 for i in range(N_STEPS + 1)]
    end25 = stats['25% PD1+ start'][0][-1]
    end75 = stats['75% PD1+ start'][0][-1]
    std25 = stats['25% PD1+ start'][1][-1]
    std75 = stats['75% PD1+ start'][1][-1]
    passed = end75 > end25
    verdict = {
        'metric': 'exhausted (PD1+) T-cell fraction at end (%), mean ± std over 6 seeds',
        'pd1p_end': {'25% start': f'{end25:.0f} ± {std25:.0f}',
                     '75% start': f'{end75:.0f} ± {std75:.0f}'},
        'behavior': 'exhaustion(75% start) > exhaustion(25% start)',
        'passed': bool(passed),
    }
    print(json.dumps(verdict, indent=2))

    viz_dir = STUDY_DIR / 'viz'; viz_dir.mkdir(parents=True, exist_ok=True)
    fig = viz.timeseries_band_figure(
        times_h, stats,
        title='T-cell exhaustion vs. starting PD1+ fraction (mean ± std, 6 seeds)',
        yaxis='exhausted (PD1+) T cells (%)')
    viz.write_html(fig, viz_dir / 'pd1p_fraction.html', STUDY_SLUG)
    fig2 = viz.spatial_animation_figure(
        first_frames['75% PD1+ start'], BOUNDS,
        title='75% PD1+ start — T cells (green=active, orange=exhausted) over IFNg')
    viz.write_html(fig2, viz_dir / 'spatial_75pct.html', STUDY_SLUG)

    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
