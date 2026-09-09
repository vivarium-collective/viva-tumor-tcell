"""Composite generators for the tumor microenvironment.

Generators
----------
* ``tumor_tcell_basic`` — small well-mixed chamber (M1 core-seam demo).
* ``tumor_microenvironment`` — CODEX-style layout: a central tumor mass with a
  ring of T cells around it, sharing a diffusing IFNg field. Parameterized by
  the PD1+ fraction of T cells and the number of T cells so the three headline
  conditions (no-T / 25% PD1+ / 75% PD1+) are just parameter choices.
* ``killing_assay`` — a well-mixed in-vitro cytotoxicity layout: tumor cells at a
  chosen PDL1+ fraction, with or without T cells (the 4-condition assay).
* ``lymph_node`` — the microenvironment plus dendritic cells and a diffusing
  tumor_debris field (dendritic cells take up debris and activate).

All use viva-munk collisions via ``TumorTcellPhysics`` and the ported
``DiffusionField``. Per-cell-type behavior processes (tumor/tcell/dendritic) operate on the whole cells map; division
adds daughters in-place. No in-document emitter (see ``viva_tumor_tcell.run``).
"""
import math
import random as _random

import numpy as np
from viva_superpowers.composite_generator import composite_generator

from ..core import build_core  # noqa: F401  (ensures type/link registration on import)
from ..types import cells_store, MOLECULES
from ..visualizations import viz_steps

TIMESTEP = 60.0


def _zero_mol():
    return {m: 0.0 for m in MOLECULES}


def _behavior_process(address, interval=TIMESTEP):
    """A top-level per-cell-type behavior process operating on the whole cells
    map (visible as its own node in the Composite Explorer)."""
    return {
        '_type': 'process',
        'address': f'local:{address}',
        'config': {},
        'interval': interval,
        'inputs': {'cells': ['cells']},
        'outputs': {'cells': ['cells']},
    }


def _tumor_cell(cid, location, state):
    return {
        'id': cid, 'type': 'circle', 'cell_type': 'tumor', 'cell_state': state,
        'death': '', 'mass': 8.0, 'radius': 7.5,
        'location': (float(location[0]), float(location[1])), 'speed': 0.0,
        'present_MHCI': 5e4 if state == 'PDL1p' else 1000.0,
        'present_PDL1': 5e4 if state == 'PDL1p' else 0.0,
        'internal_IFNg': 0.0, 'receive_cytotoxic': 0.0,
        'exchange': _zero_mol(), 'local': _zero_mol(),
    }


def _t_cell(cid, location, state):
    return {
        'id': cid, 'type': 'circle', 'cell_type': 't-cell', 'cell_state': state,
        'death': '', 'mass': 2.0, 'radius': 3.75,
        'location': (float(location[0]), float(location[1])),
        'speed': (5.0 if state == 'PD1p' else 10.0) / 60.0,
        'present_TCR': 50000.0, 'present_PD1': 5e4 if state == 'PD1p' else 0.0,
        'total_cytotoxic_packets': 0.0, 'transfer_cytotoxic': 0.0,
        'exchange': _zero_mol(), 'local': _zero_mol(),
    }


def _dendritic_cell(cid, location, state='inactive'):
    return {
        'id': cid, 'type': 'circle', 'cell_type': 'dendritic', 'cell_state': state,
        'death': '', 'mass': 2.0, 'radius': 5.0,
        'location': (float(location[0]), float(location[1])), 'speed': 3.0 / 60.0,
        'present_MHCI': 0.0, 'present_PDL1': 0.0, 'internal_tumor_debris': 0.0,
        'exchange': _zero_mol(), 'local': _zero_mol(),
    }


# ------------------------------------------------------------------ placement
def _disc(n, cx, cy, radius, rng):
    out = []
    for _ in range(int(n)):
        r = radius * math.sqrt(rng.uniform(0, 1))
        th = rng.uniform(0, 2 * math.pi)
        out.append((cx + r * math.cos(th), cy + r * math.sin(th)))
    return out


def _ring(n, cx, cy, r_inner, r_outer, rng):
    out = []
    for _ in range(int(n)):
        r = math.sqrt(rng.uniform(r_inner ** 2, r_outer ** 2))
        th = rng.uniform(0, 2 * math.pi)
        out.append((cx + r * math.cos(th), cy + r * math.sin(th)))
    return out


def _uniform(n, bx, by, rng, margin=0.1):
    return [(rng.uniform(margin * bx, (1 - margin) * bx),
             rng.uniform(margin * by, (1 - margin) * by)) for _ in range(int(n))]


def _physics(bx, by):
    return {
        '_type': 'process', 'address': 'local:TumorTcellPhysics',
        'config': {'bounds_x': bx, 'bounds_y': by, 'neighbor_distance': 1.0},
        'interval': TIMESTEP,
        'inputs': {'cells': ['cells']}, 'outputs': {'cells': ['cells']},
    }


def _field(bx, by, nx, ny, depth, molecules):
    return {
        '_type': 'process', 'address': 'local:DiffusionField',
        'config': {'bounds_x': bx, 'bounds_y': by, 'n_bins_x': nx, 'n_bins_y': ny,
                   'depth': depth, 'molecules': list(molecules)},
        'interval': TIMESTEP,
        'inputs': {'cells': ['cells'], 'fields': ['fields']},
        'outputs': {'cells': ['cells'], 'fields': ['fields']},
    }


# ============================================================ documents

def tumor_tcell_basic_document(
    n_tumors=8, n_tcells=3, tumor_pdl1n_frac=0.9, tcell_pd1n_frac=0.75,
    bounds=(120.0, 120.0), n_bins=(12, 12), depth=15.0, seed=1,
):
    rng = _random.Random(seed)
    bx, by = float(bounds[0]), float(bounds[1])
    cells = {}
    for i, loc in enumerate(_uniform(n_tumors, bx, by, rng)):
        st = 'PDL1n' if rng.uniform(0, 1) < tumor_pdl1n_frac else 'PDL1p'
        cells[f"tumor_{i}"] = _tumor_cell(f"tumor_{i}", loc, st)
    for i, loc in enumerate(_uniform(n_tcells, bx, by, rng)):
        st = 'PD1n' if rng.uniform(0, 1) < tcell_pd1n_frac else 'PD1p'
        cells[f"tcell_{i}"] = _t_cell(f"tcell_{i}", loc, st)
    nx, ny = int(n_bins[0]), int(n_bins[1])
    return {
        'cells': cells_store(cells),
        'fields': {'IFNg': np.zeros((nx, ny))},
        'physics': _physics(bx, by),
        'field': _field(bx, by, nx, ny, depth, ['IFNg']),
        'tumor_behavior': _behavior_process('TumorCellProcess'),
        'tcell_behavior': _behavior_process('TCellProcess'),
        **viz_steps(bounds=(bx, by), field='IFNg'),
    }


def tumor_microenvironment_document(
    n_tumors=60, n_tcells=6, pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9,
    bounds=(400.0, 400.0), n_bins=(40, 40), depth=15.0, seed=1,
):
    """Central tumor mass + a ring of T cells (CODEX layout).

    pd1_positive_frac is the fraction of T cells that START PD1+ (exhausted):
    the paper's 25%/75% PD1+ conditions are pd1_positive_frac 0.25 / 0.75, and
    n_tcells=0 is the no-T-cell control.
    """
    rng = _random.Random(seed)
    bx, by = float(bounds[0]), float(bounds[1])
    cx, cy = bx / 2, by / 2
    mass_radius = 0.28 * min(bx, by)

    cells = {}
    for i, loc in enumerate(_disc(n_tumors, cx, cy, mass_radius, rng)):
        st = 'PDL1n' if rng.uniform(0, 1) < tumor_pdl1n_frac else 'PDL1p'
        cells[f"tumor_{i}"] = _tumor_cell(f"tumor_{i}", loc, st)
    for i, loc in enumerate(_ring(n_tcells, cx, cy, mass_radius, mass_radius + 0.15 * min(bx, by), rng)):
        st = 'PD1p' if rng.uniform(0, 1) < pd1_positive_frac else 'PD1n'
        cells[f"tcell_{i}"] = _t_cell(f"tcell_{i}", loc, st)

    nx, ny = int(n_bins[0]), int(n_bins[1])
    return {
        'cells': cells_store(cells),
        'fields': {'IFNg': np.zeros((nx, ny))},
        'physics': _physics(bx, by),
        'field': _field(bx, by, nx, ny, depth, ['IFNg']),
        'tumor_behavior': _behavior_process('TumorCellProcess'),
        'tcell_behavior': _behavior_process('TCellProcess'),
        **viz_steps(bounds=(bx, by), field='IFNg'),
    }


def killing_assay_document(
    n_tumors=40, tumor_t_ratio=1.0, pdl1_positive_frac=0.0, include_tcells=True,
    pd1_positive_frac=0.25, bounds=(250.0, 250.0), n_bins=(25, 25), depth=15.0, seed=1,
):
    """Well-mixed in-vitro cytotoxicity layout. With include_tcells=False this is
    the matched no-T control used to compute % cytotoxicity."""
    rng = _random.Random(seed)
    bx, by = float(bounds[0]), float(bounds[1])
    cells = {}
    for i, loc in enumerate(_uniform(n_tumors, bx, by, rng)):
        st = 'PDL1p' if rng.uniform(0, 1) < pdl1_positive_frac else 'PDL1n'
        cells[f"tumor_{i}"] = _tumor_cell(f"tumor_{i}", loc, st)
    if include_tcells:
        n_tcells = max(1, int(round(n_tumors * tumor_t_ratio)))
        for i, loc in enumerate(_uniform(n_tcells, bx, by, rng)):
            st = 'PD1p' if rng.uniform(0, 1) < pd1_positive_frac else 'PD1n'
            cells[f"tcell_{i}"] = _t_cell(f"tcell_{i}", loc, st)
    nx, ny = int(n_bins[0]), int(n_bins[1])
    return {
        'cells': cells_store(cells),
        'fields': {'IFNg': np.zeros((nx, ny))},
        'physics': _physics(bx, by),
        'field': _field(bx, by, nx, ny, depth, ['IFNg']),
        'tumor_behavior': _behavior_process('TumorCellProcess'),
        'tcell_behavior': _behavior_process('TCellProcess'),
        **viz_steps(bounds=(bx, by), field='IFNg'),
    }


def lymph_node_document(
    n_tumors=40, n_tcells=6, n_dendritic=3, pd1_positive_frac=0.25,
    tumor_pdl1n_frac=0.9, bounds=(400.0, 400.0), n_bins=(40, 40), depth=15.0, seed=1,
):
    """Microenvironment + dendritic cells + a diffusing tumor_debris field.

    Dendritic cells take up debris released by dying tumors and activate — the
    antigen-presentation arm of the lymph-node extension. (Cross-compartment
    T-cell recirculation is a further extension; here all cells share one
    compartment.)"""
    doc = tumor_microenvironment_document(
        n_tumors=n_tumors, n_tcells=n_tcells, pd1_positive_frac=pd1_positive_frac,
        tumor_pdl1n_frac=tumor_pdl1n_frac, bounds=bounds, n_bins=n_bins, depth=depth, seed=seed)
    rng = _random.Random(seed + 999)
    bx, by = float(bounds[0]), float(bounds[1])
    cells = {cid: c for cid, c in doc['cells'].items() if not cid.startswith('_')}
    for i, loc in enumerate(_uniform(n_dendritic, bx, by, rng)):
        cells[f"dc_{i}"] = _dendritic_cell(f"dc_{i}", loc)
    nx, ny = int(n_bins[0]), int(n_bins[1])
    doc['cells'] = cells_store(cells)
    doc['fields'] = {'IFNg': np.zeros((nx, ny)), 'tumor_debris': np.zeros((nx, ny))}
    doc['field'] = _field(bx, by, nx, ny, depth, ['IFNg', 'tumor_debris'])
    doc['dendritic_behavior'] = _behavior_process('DendriticCellProcess')
    return doc


# ============================================================ generators

@composite_generator(
    name='tumor_tcell_basic',
    description='Small well-mixed chamber of tumor + T cells sharing a diffusing '
                'IFNg field; collisions via viva-munk (M1 core-seam demo).',
    parameters={
        'n_tumors': {'type': 'integer', 'default': 8, 'description': 'number of tumor cells'},
        'n_tcells': {'type': 'integer', 'default': 3, 'description': 'number of T cells'},
        'tumor_pdl1n_frac': {'type': 'float', 'default': 0.9, 'description': 'fraction of tumors starting PDL1n'},
        'tcell_pd1n_frac': {'type': 'float', 'default': 0.75, 'description': 'fraction of T cells starting PD1-'},
        'seed': {'type': 'integer', 'default': 1, 'description': 'RNG seed'},
    },
    default_n_steps=200,
    core_extensions=[build_core],
)
def tumor_tcell_basic(core=None, *, n_tumors=8, n_tcells=3,
                      tumor_pdl1n_frac=0.9, tcell_pd1n_frac=0.75, seed=1):
    return tumor_tcell_basic_document(
        n_tumors=n_tumors, n_tcells=n_tcells,
        tumor_pdl1n_frac=tumor_pdl1n_frac, tcell_pd1n_frac=tcell_pd1n_frac, seed=seed)


@composite_generator(
    name='tumor_microenvironment',
    description='CODEX layout: a central tumor mass ringed by T cells over a '
                'diffusing IFNg field. pd1_positive_frac + n_tcells select the '
                'no-T / 25% PD1+ / 75% PD1+ headline conditions.',
    parameters={
        'n_tumors': {'type': 'integer', 'default': 60, 'description': 'tumor cells in the central mass'},
        'n_tcells': {'type': 'integer', 'default': 6, 'description': 'T cells in the ring (0 = no-T control)'},
        'pd1_positive_frac': {'type': 'float', 'default': 0.25, 'description': 'fraction of T cells starting PD1+ (exhausted)'},
        'tumor_pdl1n_frac': {'type': 'float', 'default': 0.9, 'description': 'fraction of tumors starting PDL1n'},
        'seed': {'type': 'integer', 'default': 1, 'description': 'RNG seed'},
    },
    default_n_steps=300,
    core_extensions=[build_core],
)
def tumor_microenvironment(core=None, *, n_tumors=60, n_tcells=6,
                           pd1_positive_frac=0.25, tumor_pdl1n_frac=0.9, seed=1):
    return tumor_microenvironment_document(
        n_tumors=n_tumors, n_tcells=n_tcells, pd1_positive_frac=pd1_positive_frac,
        tumor_pdl1n_frac=tumor_pdl1n_frac, seed=seed)


@composite_generator(
    name='killing_assay',
    description='Well-mixed in-vitro cytotoxicity assay: tumor cells at a chosen '
                'PDL1+ fraction, with or without T cells (matched no-T control).',
    parameters={
        'n_tumors': {'type': 'integer', 'default': 40, 'description': 'number of tumor cells'},
        'tumor_t_ratio': {'type': 'float', 'default': 1.0, 'description': 'T:tumor ratio'},
        'pdl1_positive_frac': {'type': 'float', 'default': 0.0, 'description': 'fraction of tumors starting PDL1+'},
        'include_tcells': {'type': 'boolean', 'default': True, 'description': 'include T cells (False = matched no-T control)'},
        'pd1_positive_frac': {'type': 'float', 'default': 0.25, 'description': 'fraction of T cells starting PD1+'},
        'seed': {'type': 'integer', 'default': 1, 'description': 'RNG seed'},
    },
    default_n_steps=300,
    core_extensions=[build_core],
)
def killing_assay(core=None, *, n_tumors=40, tumor_t_ratio=1.0,
                  pdl1_positive_frac=0.0, include_tcells=True,
                  pd1_positive_frac=0.25, seed=1):
    return killing_assay_document(
        n_tumors=n_tumors, tumor_t_ratio=tumor_t_ratio,
        pdl1_positive_frac=pdl1_positive_frac, include_tcells=include_tcells,
        pd1_positive_frac=pd1_positive_frac, seed=seed)


@composite_generator(
    name='lymph_node',
    description='Tumor microenvironment plus dendritic cells and a diffusing '
                'tumor_debris field — DCs take up debris and activate.',
    parameters={
        'n_tumors': {'type': 'integer', 'default': 40, 'description': 'number of tumor cells'},
        'n_tcells': {'type': 'integer', 'default': 6, 'description': 'number of T cells'},
        'n_dendritic': {'type': 'integer', 'default': 3, 'description': 'number of dendritic cells'},
        'pd1_positive_frac': {'type': 'float', 'default': 0.25, 'description': 'fraction of T cells starting PD1+'},
        'seed': {'type': 'integer', 'default': 1, 'description': 'RNG seed'},
    },
    default_n_steps=300,
    core_extensions=[build_core],
)
def lymph_node(core=None, *, n_tumors=40, n_tcells=6, n_dendritic=3,
               pd1_positive_frac=0.25, seed=1):
    return lymph_node_document(
        n_tumors=n_tumors, n_tcells=n_tcells, n_dendritic=n_dendritic,
        pd1_positive_frac=pd1_positive_frac, seed=seed)
