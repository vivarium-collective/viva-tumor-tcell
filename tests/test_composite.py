"""Composite / generator integration tests."""
import numpy as np

from process_bigraph import Composite

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import (
    tumor_tcell_basic_document, _tumor_cell, _t_cell)
from viva_tumor_tcell.types import cells_store
from viva_tumor_tcell.run import snapshot_run


import pytest


@pytest.mark.parametrize('name', [
    'tumor_tcell_basic', 'tumor_microenvironment', 'killing_assay', 'lymph_node'])
def test_generator_registered(name):
    from process_bigraph.composite_generator import _REGISTRY
    matches = [eid for eid in _REGISTRY if eid.endswith('.' + name)]
    assert matches, f'{name} missing; have {list(_REGISTRY)[:8]}'


def test_all_generators_build_and_step():
    from viva_tumor_tcell.composites.microenvironment import (
        tumor_microenvironment_document, killing_assay_document, lymph_node_document)
    core = build_core()
    for doc_fn in (
        lambda: tumor_microenvironment_document(n_tumors=10, n_tcells=3, bounds=(150., 150.), n_bins=(15, 15)),
        lambda: killing_assay_document(n_tumors=10, bounds=(150., 150.), n_bins=(15, 15)),
        lambda: lymph_node_document(n_tumors=8, n_tcells=2, n_dendritic=2, bounds=(150., 150.), n_bins=(15, 15)),
    ):
        sim = Composite({'state': doc_fn()}, core=core)
        sim.run(60.0 * 3)   # 3 ticks, no exception


def test_basic_composite_builds_and_runs():
    core = build_core()
    doc = tumor_tcell_basic_document(n_tumors=6, n_tcells=3, seed=2)
    sim = Composite({'state': doc}, core=core)
    frames = snapshot_run(sim, 20)
    assert len(frames) == 21
    assert frames[0]['n_cells'] == 9
    # runs stably; cells remain a mix of tumor and t-cell
    assert frames[-1]['n_tumor'] >= 1


def test_composites_declare_renderable_visualizations():
    """Every composite embeds Visualization steps whose render() produces HTML —
    this is what the Composite Explorer's Visualizations tab shows."""
    from process_bigraph.visualization import render_results
    from viva_tumor_tcell.composites.microenvironment import (
        tumor_microenvironment_document, killing_assay_document, lymph_node_document)
    core = build_core()
    docs = {
        'basic': lambda: tumor_tcell_basic_document(n_tumors=6, n_tcells=3, seed=1),
        'microenv': lambda: tumor_microenvironment_document(
            n_tumors=8, n_tcells=3, bounds=(150., 150.), n_bins=(15, 15)),
        'killing': lambda: killing_assay_document(
            n_tumors=8, bounds=(150., 150.), n_bins=(15, 15)),
        'lymph_node': lambda: lymph_node_document(
            n_tumors=8, n_tcells=2, n_dendritic=2, bounds=(150., 150.), n_bins=(15, 15)),
    }
    for name, fn in docs.items():
        sim = Composite({'state': fn()}, core=core)
        for _ in range(5):
            sim.run(60.0)
        rendered = render_results(sim)
        assert len(rendered) == 3, f'{name}: expected 3 viz, got {list(rendered)}'
        for path, payload in rendered.items():
            html = payload.get('html', '')
            assert 'plotly' in html.lower(), f'{name} {path}: not a Plotly figure'


def test_ifng_feedback_and_killing_in_contact():
    """A PDL1p tumor (MHCI-high) surrounded by non-migrating T-cells: the T-cells
    secrete IFNg into the field and transfer cytotoxic packets; the tumor's
    received-packet count climbs. Exercises the full physics + exchange + field +
    behavior loop end to end. Deterministic (seeded); the T-cells are placed just
    touching the tumor (no overlap) so collisions don't eject them."""
    import math
    import random as _random
    _random.seed(11); np.random.seed(11)

    core = build_core()
    doc = tumor_tcell_basic_document(n_tumors=0, n_tcells=0, seed=7)
    cx, cy = 60.0, 60.0
    tu = _tumor_cell('tumor_A', (cx, cy), 'PDL1p')  # PDL1p -> present_MHCI 5e4
    cells = {'tumor_A': tu}
    gap = tu['radius'] + 3.75 + 0.5   # just touching (contact, no overlap)
    for k in range(6):
        th = 2 * math.pi * k / 6
        tc = _t_cell(f'tcell_{k}', (cx + gap * math.cos(th), cy + gap * math.sin(th)), 'PD1n')
        tc['speed'] = 0.0             # hold position
        cells[f'tcell_{k}'] = tc
    doc['cells'] = cells_store(cells)
    sim = Composite({'state': doc}, core=core)

    received = 0.0
    ifng_max = 0.0
    for _ in range(25):
        sim.run(60.0)
        for c in sim.state['cells'].values():
            if c.get('cell_type') == 'tumor':
                received = max(received, float(c.get('receive_cytotoxic', 0.0)))
        ifng_max = max(ifng_max, float(np.max(np.asarray(sim.state['fields']['IFNg']))))

    assert received > 0.0, 'tumor never received cytotoxic packets in contact'
    assert ifng_max > 0.0, 'no IFNg was secreted into the field'
