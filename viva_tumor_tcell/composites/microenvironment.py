"""Composite generators for the tumor microenvironment.

``tumor_tcell_basic`` — a small well-mixed chamber of tumor cells and T cells
sharing a diffusing IFNg field, with collisions via viva-munk. This is the M1
core-seam composite: it exercises the physics swap, the per-cell behaviour
processes (embedded + realized through division), the neighbor-exchange, and
the IFNg field end to end.
"""
import math
import random as _random

import numpy as np
from viva_superpowers.composite_generator import composite_generator

from ..core import build_core  # noqa: F401  (ensures type/link registration on import)
from ..types import cells_store, MOLECULES

TIMESTEP = 60.0


def _zero_mol():
    return {m: 0.0 for m in MOLECULES}


def _behavior(address, interval=TIMESTEP, config=None):
    return {
        '_type': 'process',
        'address': f'local:{address}',
        'config': config or {},
        'interval': interval,
        'inputs': {'agent_id': ['id'], 'agents': ['..', '..', 'cells']},
        'outputs': {'agents': ['..', '..', 'cells']},
    }


def _tumor_cell(cid, location, state):
    return {
        'id': cid, 'type': 'circle', 'cell_type': 'tumor', 'cell_state': state,
        'death': '', 'mass': 8.0, 'radius': 7.5,
        'location': (float(location[0]), float(location[1])), 'speed': 0.0,
        'present_MHCI': 1000.0, 'present_PDL1': 0.0,
        'internal_IFNg': 0.0, 'receive_cytotoxic': 0.0,
        'exchange': _zero_mol(), 'local': _zero_mol(),
        'behavior': _behavior('TumorCellProcess'),
    }


def _t_cell(cid, location, state):
    return {
        'id': cid, 'type': 'circle', 'cell_type': 't-cell', 'cell_state': state,
        'death': '', 'mass': 2.0, 'radius': 3.75,
        'location': (float(location[0]), float(location[1])), 'speed': 10.0 / 60.0,
        'present_TCR': 50000.0, 'present_PD1': 0.0,
        'total_cytotoxic_packets': 0.0, 'transfer_cytotoxic': 0.0,
        'exchange': _zero_mol(), 'local': _zero_mol(),
        'behavior': _behavior('TCellProcess'),
    }


def tumor_tcell_basic_document(
    n_tumors=8, n_tcells=3, tumor_pdl1n_frac=0.9, tcell_pd1n_frac=0.75,
    bounds=(120.0, 120.0), n_bins=(12, 12), depth=15.0, seed=1,
):
    rng = _random.Random(seed)
    bx, by = float(bounds[0]), float(bounds[1])
    cells = {}

    for i in range(int(n_tumors)):
        cid = f"tumor_{i}"
        loc = (rng.uniform(0.1 * bx, 0.9 * bx), rng.uniform(0.1 * by, 0.9 * by))
        state = 'PDL1n' if rng.uniform(0, 1) < tumor_pdl1n_frac else 'PDL1p'
        cells[cid] = _tumor_cell(cid, loc, state)

    for i in range(int(n_tcells)):
        cid = f"tcell_{i}"
        loc = (rng.uniform(0.1 * bx, 0.9 * bx), rng.uniform(0.1 * by, 0.9 * by))
        state = 'PD1n' if rng.uniform(0, 1) < tcell_pd1n_frac else 'PD1p'
        cells[cid] = _t_cell(cid, loc, state)

    nx, ny = int(n_bins[0]), int(n_bins[1])
    ifng = np.zeros((nx, ny), dtype=float)

    return {
        # Pin the cells map to the explicit AGENT_SCHEMA (set_float for
        # set-semantics fields, float for accumulators) — see viva_tumor_tcell.types.
        'cells': cells_store(cells),
        'fields': {'IFNg': ifng},
        'physics': {
            '_type': 'process',
            'address': 'local:TumorTcellPhysics',
            'config': {'bounds_x': bx, 'bounds_y': by, 'neighbor_distance': 1.0},
            'interval': TIMESTEP,
            'inputs': {'cells': ['cells']},
            'outputs': {'cells': ['cells']},
        },
        'field': {
            '_type': 'process',
            'address': 'local:DiffusionField',
            'config': {'bounds_x': bx, 'bounds_y': by, 'n_bins_x': nx, 'n_bins_y': ny,
                       'depth': depth, 'molecules': ['IFNg']},
            'interval': TIMESTEP,
            'inputs': {'cells': ['cells'], 'fields': ['fields']},
            'outputs': {'cells': ['cells'], 'fields': ['fields']},
        },
        # NOTE: no in-document emitter. The opaque tumor_tcell_agent map holds
        # embedded (realized) behavior processes; a framework emitter over it
        # re-infers a schema from the live state and recurses into those process
        # instances. Results are gathered with a manual per-tick snapshot loop
        # (see snapshot_run in viva_tumor_tcell.run); an M3 dashboard emitter
        # will emit a serialized snapshot store instead.
    }


@composite_generator(
    name='tumor_tcell_basic',
    description='Tumor cells + T cells in a well-mixed chamber sharing a diffusing IFNg '
                'field; collisions via viva-munk. IFNg drives the tumor PDL1n->PDL1p '
                'switch; contacting T cells transfer cytotoxic packets and kill tumors.',
    default_n_steps=200,
)
def tumor_tcell_basic(core=None, *, n_tumors=8, n_tcells=3,
                      tumor_pdl1n_frac=0.9, tcell_pd1n_frac=0.75, seed=1):
    return tumor_tcell_basic_document(
        n_tumors=n_tumors, n_tcells=n_tcells,
        tumor_pdl1n_frac=tumor_pdl1n_frac, tcell_pd1n_frac=tcell_pd1n_frac, seed=seed)
