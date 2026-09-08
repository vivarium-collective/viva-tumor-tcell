#!/usr/bin/env python3
"""Render the larger-scale efficacy result: a mean±std tumor-count band chart
per condition, from reports/scale_efficacy_results.json."""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from viva_tumor_tcell import viz

data = json.loads((REPO / 'reports' / 'scale_efficacy_results.json').read_text())
cfg = data['config']
steps = cfg['steps']
times_h = [i * 60.0 / 3600.0 for i in range(steps + 1)]

stats = {}
for label, runs in data['runs'].items():
    series = np.array([runs[s]['tumor_series'] for s in runs])  # (n_seeds, n_frames)
    stats[label] = (np.mean(series, axis=0), np.std(series, axis=0))

fig = viz.timeseries_band_figure(
    times_h, stats,
    title=f"Larger-scale efficacy: tumor count "
          f"({cfg['n_tumors']} tumors, {steps} ticks, {len(cfg['seeds'])} seeds, mean ± std)",
    yaxis='tumor count')
out = REPO / 'reports' / 'scale_efficacy.html'
viz.write_html(fig, out, 'scale-efficacy')
print('wrote', out)
if 'summary' in data:
    print(json.dumps(data['summary'], indent=2))
