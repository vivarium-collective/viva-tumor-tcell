#!/usr/bin/env python3
"""Persist one representative simulation per study to ``studies/<slug>/runs.db``.

The workbench study **Results** tab previews the latest run's persisted emitter
store (per-store sparkline + first/last/min/max) and lists downloadable raw
stores. The replicate study runners (``studies/<slug>/sims/run.py``) report
mean ± std but emit to RAM — they persist nothing, so the Results tab is empty.

This script runs ONE representative simulation of each study's headline
condition with an Observables step + a SQLiteEmitter (via
``studies_lib.persist_run``), writing a scalar-per-step history the Results tab
reads. Deterministic (seed 1). Run from the workspace root:

    python scripts/persist_study_runs.py
"""
from __future__ import annotations

import random
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import (
    tumor_microenvironment_document, killing_assay_document)
from viva_tumor_tcell.studies_lib import persist_run

INV = 'tumor-tcell-showcase'

# (slug, doc builder, n_steps, run label)
STUDIES = [
    ('tumor-microenvironment',
     lambda: tumor_microenvironment_document(
         n_tumors=60, n_tcells=6, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
         bounds=(400., 400.), n_bins=(40, 40), seed=1),
     300, 'tumor-microenvironment: 25% PD1+ (representative)'),
    ('phenotype-conversion',
     lambda: tumor_microenvironment_document(
         n_tumors=40, n_tcells=12, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
         bounds=(300., 300.), n_bins=(30, 30), seed=1),
     400, 'phenotype-conversion: with active T cells (representative)'),
    ('tcell-exhaustion',
     lambda: tumor_microenvironment_document(
         n_tumors=40, n_tcells=12, pd1_positive_frac=0.75, tumor_pdl1n_frac=0.9,
         bounds=(300., 300.), n_bins=(30, 30), seed=1),
     400, 'tcell-exhaustion: 75% PD1+ start (representative)'),
    ('killing-assay-cytotoxicity',
     lambda: killing_assay_document(
         n_tumors=40, tumor_t_ratio=1.5, pdl1_positive_frac=0.5, include_tcells=True,
         pd1_positive_frac=0.25, bounds=(250., 250.), n_bins=(25, 25), seed=1),
     300, 'killing-assay (+T): well-mixed cytotoxicity (representative)'),
    ('efficacy-at-scale',
     lambda: tumor_microenvironment_document(
         n_tumors=80, n_tcells=16, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
         bounds=(450., 450.), n_bins=(45, 45), seed=1),
     300, 'efficacy-at-scale: 25% PD1+ larger scale (representative)'),
]


def main() -> int:
    core = build_core()
    for slug, build_doc, n_steps, label in STUDIES:
        study_dir = WORKSPACE_ROOT / 'studies' / slug
        db = study_dir / 'runs.db'
        if db.exists():
            db.unlink()   # one representative run per study; replace prior
        random.seed(1)
        sim_id = persist_run(study_dir, core, build_doc(), n_steps, label,
                             investigation=INV)
        size = db.stat().st_size if db.exists() else 0
        print(f'  {slug}: persisted {sim_id[:8]} ({n_steps} steps) -> runs.db ({size} bytes)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
