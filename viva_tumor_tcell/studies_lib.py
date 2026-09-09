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
