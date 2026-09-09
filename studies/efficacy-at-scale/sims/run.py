#!/usr/bin/env python3
"""Canonical run: T-cell efficacy vs. PD1+ fraction, at larger scale.

The tractable studies (40 cells, 600 ticks) could not separate population
tumor-count suppression from growth noise. This study re-runs the three
conditions at 150 cells / 1500 ticks across 3 seeds, where the paper's
efficacy-vs-exhaustion ordering emerges cleanly.

The heavy simulation (~30 min CPU) is ``scripts/scale_efficacy.py``, which writes
``reports/scale_efficacy_results.json`` (committed as the data of record). This
runner renders the mean±std figure from that committed data and re-derives the
verdict — fast and deterministic. To regenerate the data:

    python scripts/scale_efficacy.py 150 40 1500 600 60 3
"""
from __future__ import annotations

import json
import time
import uuid
from pathlib import Path

import numpy as np

STUDY_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_DIR.parents[1]

from viva_tumor_tcell import viz
from vivarium_workbench.lib.run_log import append_run_event

STUDY_SLUG = 'efficacy-at-scale'
INVESTIGATION_SLUG = 'tumor-tcell-showcase'
SPEC_ID = 'viva_tumor_tcell.composites.microenvironment.tumor_microenvironment'
DATA = WORKSPACE_ROOT / 'reports' / 'scale_efficacy_results.json'


def main() -> int:
    data = json.loads(DATA.read_text())
    cfg = data['config']
    steps = cfg['steps']
    times_h = [i * 60.0 / 3600.0 for i in range(steps + 1)]

    stats, finals = {}, {}
    for label, runs in data['runs'].items():
        series = np.array([runs[s]['tumor_series'] for s in runs], float)
        stats[label] = (np.mean(series, axis=0), np.std(series, axis=0))
        finals[label] = series[:, -1]

    noT, p25, p75 = finals['no T cells'], finals['25% PD1+'], finals['75% PD1+']
    s25 = (noT.mean() - p25) / noT.mean() * 100
    s75 = (noT.mean() - p75) / noT.mean() * 100
    ordering = int(np.sum(s25 >= s75))
    passed = ordering == len(p25)   # 25% suppresses more than 75% in every seed
    verdict = {
        'scale': f"{cfg['n_tumors']} tumors, {steps} ticks, {len(cfg['seeds'])} seeds",
        'final_tumor_mean_std': {k: f'{v.mean():.0f} ± {v.std():.0f}' for k, v in finals.items()},
        'suppression_pct': {'25% PD1+': f'{s25.mean():.1f} ± {s25.std():.1f}',
                            '75% PD1+': f'{s75.mean():.1f} ± {s75.std():.1f}'},
        'ordering_25_suppresses_more_than_75': f'{ordering}/{len(p25)} seeds',
        'behavior': 'active (25% PD1+) T cells suppress tumor growth more than exhausted (75% PD1+)',
        'passed': bool(passed),
    }

    run_id = uuid.uuid4().hex
    append_run_event(WORKSPACE_ROOT, {
        'run_id': run_id, 'event': 'started', 'spec_id': SPEC_ID,
        'label': f'{STUDY_SLUG}: render from committed scale data', 'started_at': time.time(),
        'status': 'running', 'n_steps': steps, 'emitter': 'ram', 'origin': 'canonical_run',
        'study_slug': STUDY_SLUG, 'investigation_slug': INVESTIGATION_SLUG,
        'params': cfg, 'replicates': len(cfg['seeds'])})

    viz_dir = STUDY_DIR / 'viz'; viz_dir.mkdir(parents=True, exist_ok=True)
    fig = viz.timeseries_band_figure(
        times_h, stats,
        title=f"Efficacy at scale: tumor count "
              f"({cfg['n_tumors']} tumors, {steps} ticks, {len(cfg['seeds'])} seeds, mean ± std)",
        yaxis='tumor count')
    viz.write_html(fig, viz_dir / 'tumor_count_at_scale.html', STUDY_SLUG)

    append_run_event(WORKSPACE_ROOT, {
        'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
        'n_steps': steps, 'status': 'completed'})

    print(json.dumps(verdict, indent=2))
    (STUDY_DIR / 'results.json').write_text(json.dumps(verdict, indent=2))
    from viva_tumor_tcell.studies_lib import publish_figures
    publish_figures(STUDY_DIR)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
