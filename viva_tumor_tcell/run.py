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


# ---- per-state population buckets (matches tumor-tcell population_group_plot) ----
TUMOR_STATES = ('PDL1n', 'PDL1p')
TCELL_STATES = ('PD1n', 'PD1p')
DC_STATES = ('inactive', 'active')


def population_counts(frame):
    """Return counts by cell_type and cell_state for one frame."""
    c = {'tumor_total': 0, 'tcell_total': 0, 'dendritic_total': 0}
    for st in TUMOR_STATES:
        c[f'tumor_{st}'] = 0
    for st in TCELL_STATES:
        c[f'tcell_{st}'] = 0
    for st in DC_STATES:
        c[f'dendritic_{st}'] = 0
    for cell in frame['cells'].values():
        ct, cs = cell.get('cell_type'), cell.get('cell_state')
        if ct == 'tumor':
            c['tumor_total'] += 1
            if cs in TUMOR_STATES:
                c[f'tumor_{cs}'] += 1
        elif ct == 't-cell':
            c['tcell_total'] += 1
            if cs in TCELL_STATES:
                c[f'tcell_{cs}'] += 1
        elif ct == 'dendritic':
            c['dendritic_total'] += 1
            if cs in DC_STATES:
                c[f'dendritic_{cs}'] += 1
    return c


def analysis_run(sim, n_steps, interval=60.0, keep_snapshots=8):
    """Run tick-by-tick, tracking the observables the tumor-tcell analyses need:

      * per-state population counts over time (population_group_plot),
      * cumulative deaths by death-reason type over time (death_group_plot) —
        detected before the DiffusionField removes the dead cell,
      * cumulative divisions per cell type over time (division_plot) — detected
        from the mother -> mother+'A'/'B' phylogeny id scheme,
      * peak IFNg over time,
      * a handful of evenly-spaced full snapshots for spatial panels / animation.

    Returns a dict of time series + snapshots.
    """
    import numpy as np
    times, pops, ifng = [], [], []
    death_cum = {}                 # reason -> running count
    death_series = []              # list of dict copies over time
    div_cum = {'tumor': 0, 't-cell': 0, 'dendritic': 0}
    div_series = []
    counted_deaths = set()         # cell ids already counted as dead
    prev_ids = set(sim.state['cells'].keys())
    snap_idx = set(np.linspace(0, n_steps, keep_snapshots, dtype=int).tolist())
    snapshots = []

    def _record(step):
        cells = sim.state['cells']
        t = float(sim.state.get('global_time', step * interval))
        times.append(t / 3600.0)
        pops.append(population_counts(snapshot(sim)))
        ifng.append(ifng_max(snapshot(sim)))
        # deaths: cells currently flagged dead, not yet counted
        for cid, c in cells.items():
            if c.get('death') and cid not in counted_deaths:
                counted_deaths.add(cid)
                reason = str(c.get('death'))
                death_cum[reason] = death_cum.get(reason, 0) + 1
        death_series.append(dict(death_cum))
        # divisions: new A/B pairs since last tick
        cur_ids = set(cells.keys())
        new_ids = cur_ids - prev_ids
        for nid in new_ids:
            if nid.endswith('A') and (nid[:-1] + 'B') in new_ids:
                ct = cells[nid].get('cell_type', 'tumor')
                div_cum[ct] = div_cum.get(ct, 0) + 1
        div_series.append(dict(div_cum))
        prev_ids.clear(); prev_ids.update(cur_ids)
        if step in snap_idx:
            snapshots.append(snapshot(sim))

    _record(0)
    for step in range(1, n_steps + 1):
        sim.run(float(interval))
        _record(step)

    return {'times_h': times, 'populations': pops, 'ifng': ifng,
            'deaths': death_series, 'divisions': div_series, 'snapshots': snapshots}
