"""Shared helpers for the showcase study runners: a fixed replicate-seed set
and a run-recording wrapper that logs started/completed events to
``.pbg/runs.jsonl`` around each condition."""
from __future__ import annotations

import time
import uuid

from vivarium_workbench.lib.run_log import append_run_event

SEEDS = (1, 2, 3, 4, 5, 6)


def publish_figures(study_dir):
    """Mirror a study's viz/*.html into reports/figures/<slug>/ so the workbench
    discovers them WITHOUT a runs.db (studies/<slug>/viz needs a workbench run;
    reports/figures/<slug> does not — see study_spec.discover_viz_html_files)."""
    import shutil
    from pathlib import Path
    study_dir = Path(study_dir)
    slug = study_dir.name
    ws_root = study_dir.parents[1]
    dst = ws_root / 'reports' / 'figures' / slug
    dst.mkdir(parents=True, exist_ok=True)
    for f in (study_dir / 'viz').glob('*.html'):
        shutil.copy2(f, dst / f.name)


def record_run(root, slug, investigation, spec_id, label, params, n_steps, fn):
    """Emit started/completed events around ``fn()`` and return its result."""
    run_id = uuid.uuid4().hex
    append_run_event(root, {
        'run_id': run_id, 'event': 'started', 'spec_id': spec_id,
        'label': f'{slug}: {label}', 'started_at': time.time(), 'status': 'running',
        'n_steps': n_steps, 'emitter': 'ram', 'origin': 'canonical_run',
        'study_slug': slug, 'investigation_slug': investigation,
        'params': params, 'replicates': len(SEEDS)})
    try:
        result = fn()
    except Exception:
        append_run_event(root, {
            'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
            'n_steps': 0, 'status': 'failed'})
        raise
    append_run_event(root, {
        'run_id': run_id, 'event': 'completed', 'completed_at': time.time(),
        'n_steps': n_steps, 'status': 'completed'})
    return result


# Scalar observables persisted for the Results-tab preview. global_time is added
# automatically; these come from the Observables step's `observables` store.
_OBS_KEYS = (
    'tumor_total', 'tumor_PDL1n', 'tumor_PDL1p',
    'tcell_total', 'tcell_PD1n', 'tcell_PD1p',
    'dendritic_total', 'dendritic_inactive', 'dendritic_active',
    'pdl1p_fraction', 'pd1p_fraction', 'ifng_max',
)


def _emitter_node(study_dir, simulation_id, name):
    """A SQLiteEmitter step writing global_time + observables to
    ``<study_dir>/runs.db`` — the store the workbench Results tab reads."""
    from pathlib import Path
    study_dir = Path(study_dir)
    emit = {'global_time': 'float',
            'observables': {'_type': 'map', '_value': 'set_float'}}
    return {
        '_type': 'step', 'address': 'local:SQLiteEmitter',
        'config': {'emit': emit, 'file_path': str(study_dir),
                   'db_file': 'runs.db', 'simulation_id': simulation_id, 'name': name},
        'inputs': {'global_time': ['global_time'], 'observables': ['observables']},
    }


def persist_run(study_dir, core, doc, n_steps, name, interval=60.0,
                investigation=None):
    """Run one representative simulation of ``doc`` with an Observables step + a
    SQLiteEmitter, persisting a scalar per-step history to ``studies/<slug>/runs.db``.

    This is what populates the study Results tab (per-store sparklines +
    downloadable raw store). Returns the simulation_id. Kept separate from the
    replicate study runs (which report mean ± std) — one representative run is
    enough for the Results preview.
    """
    from pathlib import Path
    from process_bigraph import Composite
    try:
        from pbg_emitters.sqlite_emitter import (
            SQLiteEmitter, save_simulation_metadata)
    except ImportError:  # legacy location
        from process_bigraph.emitter import SQLiteEmitter  # type: ignore
        save_simulation_metadata = None
    core.register_link('SQLiteEmitter', SQLiteEmitter)

    study_dir = Path(study_dir)
    slug = study_dir.name
    sim_id = uuid.uuid4().hex

    doc = dict(doc)
    doc['observables'] = {k: 0.0 for k in _OBS_KEYS}
    doc['obs_step'] = {
        '_type': 'step', 'address': 'local:Observables', 'config': {},
        'inputs': {'cells': ['cells'], 'fields': ['fields']},
        'outputs': {'observables': ['observables']},
    }
    doc['emitter'] = _emitter_node(study_dir, sim_id, name)

    sim = Composite({'state': doc}, core=core)
    for _ in range(int(n_steps)):
        sim.run(float(interval))

    if save_simulation_metadata is not None:
        db_path = str(study_dir / 'runs.db')
        meta = {'study_slug': slug, 'n_steps': int(n_steps), 'interval': interval}
        if investigation:
            meta['investigation_slug'] = investigation
        save_simulation_metadata(db_path, sim_id, name=name, metadata=meta)
    return sim_id
