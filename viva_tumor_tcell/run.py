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
