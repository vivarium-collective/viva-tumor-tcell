"""Snapshot-based run helper.

The tumor_tcell_agent map holds embedded (realized) per-cell behavior processes,
which makes a framework emitter over the map recurse into process instances. So
we drive the composite tick-by-tick and snapshot a lightweight, JSON-able view
of the cells + fields each tick. This is the data source for tests and study viz
until the M3 dashboard snapshot-emitter lands.
"""
from __future__ import annotations

import numpy as np

# Scalar / vector fields worth recording per cell for analysis + viz.
_SNAPSHOT_KEYS = (
    'cell_type', 'cell_state', 'location', 'radius',
    'internal_IFNg', 'receive_cytotoxic', 'total_cytotoxic_packets',
    'transfer_cytotoxic', 'accept_MHCI', 'accept_PDL1', 'accept_TCR',
    'present_MHCI', 'present_PDL1', 'present_TCR', 'TCR_timer', 'velocity_timer',
    'PD1n_divide_count', 'refractory_count', 'death',
)


def _cell_snapshot(cell):
    snap = {}
    for k in _SNAPSHOT_KEYS:
        if k in cell:
            v = cell[k]
            snap[k] = list(v) if isinstance(v, tuple) else v
    return snap


def snapshot(sim):
    """Return a JSON-able snapshot of the current composite state."""
    cells = {cid: _cell_snapshot(c) for cid, c in sim.state['cells'].items()}
    fields = {m: np.asarray(f).tolist() for m, f in sim.state.get('fields', {}).items()}
    return {
        'time': float(sim.state.get('global_time', 0.0)),
        'n_cells': len(cells),
        'n_tumor': sum(1 for c in cells.values() if c.get('cell_type') == 'tumor'),
        'n_tcell': sum(1 for c in cells.values() if c.get('cell_type') == 't-cell'),
        'cells': cells,
        'fields': fields,
    }


def snapshot_run(sim, n_steps, interval=60.0):
    """Run the composite for n_steps ticks, snapshotting each tick.

    Returns a list of per-tick snapshots (including the initial state).
    """
    frames = [snapshot(sim)]
    for _ in range(int(n_steps)):
        sim.run(float(interval))
        frames.append(snapshot(sim))
    return frames


# ------------------------------------------------------------------ observables
def _tumors(frame):
    return [c for c in frame['cells'].values() if c.get('cell_type') == 'tumor']


def _tcells(frame):
    return [c for c in frame['cells'].values() if c.get('cell_type') == 't-cell']


def tumor_count(frame):
    return len(_tumors(frame))


def tcell_count(frame):
    return len(_tcells(frame))


def pdl1p_fraction(frame):
    tum = _tumors(frame)
    if not tum:
        return 0.0
    return sum(1 for c in tum if c.get('cell_state') == 'PDL1p') / len(tum)


def pd1p_fraction(frame):
    tc = _tcells(frame)
    if not tc:
        return 0.0
    return sum(1 for c in tc if c.get('cell_state') == 'PD1p') / len(tc)


def ifng_max(frame):
    import numpy as np
    arr = frame['fields'].get('IFNg')
    return float(np.max(np.asarray(arr))) if arr is not None else 0.0


def run_condition(core, doc, n_steps, interval=60.0):
    """Build a Composite from a document and snapshot-run it."""
    from process_bigraph import Composite
    sim = Composite({'state': doc}, core=core)
    return snapshot_run(sim, n_steps, interval)


def replicate_series(core, doc_fn, n_steps, seeds, observables, interval=60.0):
    """Run one condition across replicate seeds.

    doc_fn(seed) -> composite document. ``observables`` is either a single
    callable ``frame -> value`` or a dict ``{name: callable}``. Returns
    ``(arrays, first_frames)`` where ``arrays`` matches the ``observables`` shape:
    a single (n_seeds, n_frames) array for a callable, or ``{name: array}`` for a
    dict. ``first_frames`` are the first seed's frames (for a spatial animation).
    """
    import random
    import numpy as np
    single = callable(observables)
    obs = {'_': observables} if single else dict(observables)
    per_seed = {k: [] for k in obs}
    first_frames = None
    for s in seeds:
        random.seed(s)
        frames = run_condition(core, doc_fn(s), n_steps, interval)
        for k, fn in obs.items():
            per_seed[k].append([fn(f) for f in frames])
        if first_frames is None:
            first_frames = frames
    arrays = {k: np.array(v, dtype=float) for k, v in per_seed.items()}
    return (arrays['_'] if single else arrays), first_frames


def mean_std(arr):
    """(mean, std) over the seed axis of an (n_seeds, n_frames) array."""
    import numpy as np
    return np.mean(arr, axis=0), np.std(arr, axis=0)
