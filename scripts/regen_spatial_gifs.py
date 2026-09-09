#!/usr/bin/env python3
"""Regenerate just the spatial GIFs for the killing-assay + tumor-microenvironment
studies, without re-running the heavy 6-seed replicate series.

Memory-safe: each GIF comes from ONE single-seed sim, snapshotted at a STRIDE
(so only ~n_frames snapshots are held, not every tick), and each sim's frames are
freed before the next. This is the lean path used when a full study re-run
(studies/<slug>/sims/run.py) is too heavy for the machine.

Configs/seeds match the study runners so the GIFs are consistent with the
committed figures. Run from the workspace root:

    python scripts/regen_spatial_gifs.py [killing|microenv|all]
"""
from __future__ import annotations

import gc
import random
import sys
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]

from process_bigraph import Composite
from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import (
    killing_assay_document, tumor_microenvironment_document)
from viva_tumor_tcell import viz
from viva_tumor_tcell.run import snapshot, _cell_snapshot  # noqa: F401


def _strided_frames(sim, n_steps, n_frames=100, interval=60.0):
    """Run n_steps ticks but keep only ~n_frames snapshots (every stride-th),
    so memory stays bounded regardless of run length. Returns
    ``(frames, tumor_kills)`` — kills are counted EVERY tick (a killed tumor is
    flagged death='Tcell_death' then removed within a tick, so strided frames
    miss most death flags; the per-tick scan is cheap and accurate)."""
    stride = max(1, int(n_steps) // int(n_frames))
    frames = [snapshot(sim)]
    killed: set = set()
    for i in range(int(n_steps)):
        sim.run(float(interval))
        for cid, c in sim.state['cells'].items():
            if (isinstance(c, dict) and c.get('cell_type') == 'tumor'
                    and c.get('death') == 'Tcell_death'):
                killed.add(cid)
        if (i + 1) % stride == 0:
            frames.append(snapshot(sim))
    return frames, len(killed)


def _write(study, name, frames, bounds, title, field='IFNg'):
    out = WORKSPACE_ROOT / 'studies' / study / 'viz' / f'{name}.html'
    viz.write_html_str(viz.spatial_gif_html(frames, bounds, title, field=field,
                                            max_frames=len(frames) + 1), out)
    print(f"  wrote {out.relative_to(WORKSPACE_ROOT)} ({len(frames)} frames, "
          f"{out.stat().st_size // 1024} KB)")



def do_killing(core):
    print("killing-assay-cytotoxicity:")
    # spatial_killing — full +T assay run (matches runner: seed 1, N_STEPS 700)
    random.seed(1)
    sim = Composite({'state': killing_assay_document(
        n_tumors=30, include_tcells=True, pdl1_positive_frac=0.5, tumor_t_ratio=1.5,
        pd1_positive_frac=0.25, bounds=(180.0, 180.0), n_bins=(18, 18), seed=1)}, core=core)
    frames, _ = _strided_frames(sim, 700)
    _write('killing-assay-cytotoxicity', 'spatial_killing', frames, (180.0, 180.0),
           'Killing assay (+T) — well-mixed cytotoxicity (full run)')
    del sim, frames; gc.collect()
    # spatial_killing_demo — the tuned killing_demo composite
    random.seed(1)
    sim = Composite({'state': killing_assay_document(
        n_tumors=30, tumor_t_ratio=2.0, pdl1_positive_frac=1.0, include_tcells=True,
        pd1_positive_frac=0.0, bounds=(200.0, 200.0), n_bins=(20, 20), seed=1)}, core=core)
    frames, nk = _strided_frames(sim, 500)
    _write('killing-assay-cytotoxicity', 'spatial_killing_demo', frames, (200.0, 200.0),
           f'Killing demo — PDL1+ tumors killed by active T cells ({nk} kills)')
    del sim, frames; gc.collect()


def do_microenv(core):
    print("tumor-microenvironment:")
    # spatial_animation — headline 25% PD1+ (matches runner: seed 5, N_STEPS 600)
    random.seed(5)
    sim = Composite({'state': tumor_microenvironment_document(
        n_tumors=60, n_tcells=12, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
        bounds=(400.0, 400.0), n_bins=(40, 40), seed=5)}, core=core)
    frames, _ = _strided_frames(sim, 600)
    _write('tumor-microenvironment', 'spatial_animation', frames, (400.0, 400.0),
           'Microenvironment (25% PD1+) — cells over IFNg field (full run)')
    del sim, frames; gc.collect()
    # spatial_large_scale — 150 tumors + 30 T in a 650 um chamber
    random.seed(5)
    sim = Composite({'state': tumor_microenvironment_document(
        n_tumors=150, n_tcells=30, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
        bounds=(650.0, 650.0), n_bins=(65, 65), seed=5)}, core=core)
    frames, _ = _strided_frames(sim, 600)
    _write('tumor-microenvironment', 'spatial_large_scale', frames, (650.0, 650.0),
           'Larger-scale microenvironment (150 tumors, 30 T cells) over IFNg field')
    del sim, frames; gc.collect()


def do_phenotype(core):
    print("phenotype-conversion:")
    random.seed(1)
    sim = Composite({'state': tumor_microenvironment_document(
        n_tumors=40, n_tcells=12, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
        bounds=(300.0, 300.0), n_bins=(30, 30), seed=1)}, core=core)
    frames, _ = _strided_frames(sim, 600)
    _write('phenotype-conversion', 'spatial_conversion', frames, (300.0, 300.0),
           'Active T cells build a local IFNg field (full run)')
    del sim, frames; gc.collect()


def do_exhaustion(core):
    print("tcell-exhaustion:")
    random.seed(1)
    sim = Composite({'state': tumor_microenvironment_document(
        n_tumors=40, n_tcells=12, pd1_positive_frac=0.75, tumor_pdl1n_frac=0.9,
        bounds=(300.0, 300.0), n_bins=(30, 30), seed=1)}, core=core)
    frames, _ = _strided_frames(sim, 600)
    _write('tcell-exhaustion', 'spatial_75pct', frames, (300.0, 300.0),
           '75% PD1+ start — T cells (green=active, orange=exhausted) over IFNg (full run)')
    del sim, frames; gc.collect()


def main() -> int:
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    core = build_core()
    if which in ('killing', 'all'):
        do_killing(core)
    if which in ('microenv', 'all'):
        do_microenv(core)
    if which in ('phenotype', 'all'):
        do_phenotype(core)
    if which in ('exhaustion', 'all'):
        do_exhaustion(core)
    # mirror to reports/figures/<slug>/ for the read-only dashboard
    from viva_tumor_tcell.studies_lib import publish_figures
    for slug in ('killing-assay-cytotoxicity', 'tumor-microenvironment',
                 'phenotype-conversion', 'tcell-exhaustion'):
        publish_figures(WORKSPACE_ROOT / 'studies' / slug)
    print("done")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
