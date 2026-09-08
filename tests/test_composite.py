"""Composite / generator integration tests."""
import numpy as np

from process_bigraph import Composite

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.composites.microenvironment import (
    tumor_tcell_basic_document, _tumor_cell, _t_cell)
from viva_tumor_tcell.types import cells_store
from viva_tumor_tcell.run import snapshot_run


def test_generator_registered():
    from process_bigraph.composite_generator import _REGISTRY
    matches = [eid for eid in _REGISTRY if eid.endswith('.tumor_tcell_basic')]
    assert matches, f'tumor_tcell_basic missing; have {list(_REGISTRY)[:5]}'


def test_basic_composite_builds_and_runs():
    core = build_core()
    doc = tumor_tcell_basic_document(n_tumors=6, n_tcells=3, seed=2)
    sim = Composite({'state': doc}, core=core)
    frames = snapshot_run(sim, 20)
    assert len(frames) == 21
    assert frames[0]['n_cells'] == 9
    # runs stably; cells remain a mix of tumor and t-cell
    assert frames[-1]['n_tumor'] >= 1


def test_ifng_feedback_and_killing_in_contact():
    """A PDL1p tumor (MHCI-high) touching a non-migrating T-cell: the T-cell
    produces IFNg into the field and transfers cytotoxic packets; the tumor's
    received-packet count climbs. Exercises the full physics + exchange + field
    + behavior loop end to end."""
    core = build_core()
    doc = tumor_tcell_basic_document(n_tumors=0, n_tcells=0, seed=7)
    tu = _tumor_cell('tumor_A', (60.0, 60.0), 'PDL1p')
    tu['present_MHCI'] = 5e4         # PDL1p tumor presents high MHCI
    tc = _t_cell('tcell_A', (60.0, 68.5), 'PD1n')
    tc['speed'] = 0.0                # hold contact (no migration)
    doc['cells'] = cells_store({'tumor_A': tu, 'tcell_A': tc})
    sim = Composite({'state': doc}, core=core)

    received = 0.0
    ifng_max = 0.0
    for _ in range(30):
        sim.run(60.0)
        cells = sim.state['cells']
        # follow whichever tumor id currently exists (division renames)
        for cid, c in cells.items():
            if c.get('cell_type') == 'tumor':
                received = max(received, float(c.get('receive_cytotoxic', 0.0)))
        ifng_max = max(ifng_max, float(np.max(np.asarray(sim.state['fields']['IFNg']))))

    assert received > 0.0, 'tumor never received cytotoxic packets in contact'
    assert ifng_max > 0.0, 'no IFNg was secreted into the field'
