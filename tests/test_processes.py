"""Unit tests for the ported tumor-tcell processes (deterministic).

Stochastic death/division are suppressed via tiny-probability configs so each
biological mechanism is exercised in isolation.
"""
import numpy as np
import pytest

from viva_tumor_tcell.core import build_core
from viva_tumor_tcell.processes.tumor import TumorCellProcess
from viva_tumor_tcell.processes.t_cell import TCellProcess
from viva_tumor_tcell.processes.physics import TumorTcellPhysics
from viva_tumor_tcell.processes.field import DiffusionField

TINY = 1e-12


@pytest.fixture(scope='module')
def core():
    return build_core()


def _tumor(**over):
    a = {'id': 'u', 'cell_type': 'tumor', 'cell_state': 'PDL1n', 'death': '',
         'receive_cytotoxic': 0.0, 'internal_IFNg': 0.0, 'radius': 7.5,
         'present_MHCI': 1000.0, 'present_PDL1': 0.0,
         'local': {'IFNg': 0.0, 'tumor_debris': 0.0}}
    a.update(over)
    return a


def _tcell(**over):
    a = {'id': 't', 'cell_type': 't-cell', 'cell_state': 'PD1n', 'death': '',
         'accept_MHCI': 0.0, 'accept_PDL1': 0.0, 'present_TCR': 50000.0,
         'TCR_timer': 0.0, 'velocity_timer': 0.0, 'refractory_count': 0.0,
         'PD1n_divide_count': 0.0, 'total_cytotoxic_packets': 0.0,
         'radius': 3.75, 'speed': 0.0}
    a.update(over)
    return a


# ---------------------------------------------------------------- tumor
def test_tumor_ifng_drives_pdl1p_switch(core):
    p = TumorCellProcess(config={'death_apoptosis': TINY, 'PDL1n_growth': TINY}, core=core)
    out = p.update({'agent_id': 'u', 'agents': {'u': _tumor(internal_IFNg=16000.0)}}, 60.0)['agents']['u']
    assert out['cell_state'] == 'PDL1p'
    assert out['present_MHCI'] == pytest.approx(5e4)
    assert out['present_PDL1'] == pytest.approx(5e4)


def test_tumor_killed_by_cytotoxic_packets(core):
    p = TumorCellProcess(config={'death_apoptosis': TINY}, core=core)
    out = p.update({'agent_id': 'u', 'agents': {'u': _tumor(receive_cytotoxic=13000.0)}}, 60.0)['agents']['u']
    assert out['death'] == 'Tcell_death'
    assert out['exchange']['tumor_debris'] > 0


def test_tumor_pdl1n_internalizes_available_ifng(core):
    # local IFNg present -> PDL1n tumor internalizes some counts (degrade branch)
    p = TumorCellProcess(config={'death_apoptosis': TINY, 'PDL1n_growth': TINY}, core=core)
    tu = _tumor(local={'IFNg': 5.0, 'tumor_debris': 0.0})
    out = p.update({'agent_id': 'u', 'agents': {'u': tu}}, 60.0)['agents']['u']
    assert out['internal_IFNg'] > 0
    assert out['exchange']['IFNg'] < 0  # removed from the environment


# ---------------------------------------------------------------- t cell
def test_tcell_max_production_against_pdl1p_tumor(core):
    p = TCellProcess(config={'death_PD1n_14hr': TINY, 'PD1n_growth_28hr': TINY}, core=core)
    out = p.update({'agent_id': 't', 'agents': {'t': _tcell(accept_MHCI=5e4)}}, 60.0)['agents']['t']
    assert out['exchange']['IFNg'] == 270           # 1.62e4/3600 * 60
    assert out['transfer_cytotoxic'] == pytest.approx(40.0)  # 40/60 * 60


def test_tcell_no_contact_migrates_not_kills(core):
    p = TCellProcess(config={'death_PD1n_14hr': TINY, 'PD1n_growth_28hr': TINY}, core=core)
    out = p.update({'agent_id': 't', 'agents': {'t': _tcell(accept_MHCI=0.0)}}, 60.0)['agents']['t']
    assert out.get('transfer_cytotoxic', 0.0) == 0.0
    assert out['speed'] > 0     # migrating


def test_tcell_exhausts_to_pd1p_after_refractory_cycles(core):
    p = TCellProcess(config={'death_PD1n_14hr': TINY, 'PD1n_growth_28hr': TINY}, core=core)
    # refractory_count above threshold triggers PD1n -> PD1p
    out = p.update({'agent_id': 't', 'agents': {'t': _tcell(accept_MHCI=5e4, refractory_count=4)}}, 60.0)['agents']['t']
    assert out['cell_state'] == 'PD1p'


# ---------------------------------------------------------------- field
def test_field_deposits_diffuses_decays_samples(core):
    f = DiffusionField(config={'bounds_x': 120.0, 'bounds_y': 120.0,
                               'n_bins_x': 12, 'n_bins_y': 12, 'depth': 15.0,
                               'molecules': ['IFNg']}, core=core)
    cells = {'u': {'cell_type': 'tumor', 'death': '', 'location': (60.0, 60.0),
                   'exchange': {'IFNg': 1e6}, 'local': {'IFNg': 0.0}}}
    fields = {'IFNg': np.zeros((12, 12))}
    out = f.update({'cells': cells, 'fields': fields}, 60.0)
    assert float(np.max(out['fields']['IFNg'])) > 0            # deposited concentration
    assert out['cells']['u']['local']['IFNg'] >= 0            # sampled back
    assert out['cells']['u']['exchange']['IFNg'] < 0          # exchange reset (delta)


# ---------------------------------------------------------------- physics
def test_physics_detects_contact_and_exchanges_ligands(core):
    ph = TumorTcellPhysics(config={'bounds_x': 120.0, 'bounds_y': 120.0,
                                   'neighbor_distance': 1.0}, core=core)
    cells = {
        'u': {'cell_type': 'tumor', 'location': (60.0, 60.0), 'radius': 7.5,
              'mass': 8.0, 'speed': 0.0, 'present_MHCI': 5e4, 'present_PDL1': 0.0,
              'transfer_cytotoxic': 0.0},
        't': {'cell_type': 't-cell', 'location': (60.0, 71.0), 'radius': 3.75,
              'mass': 2.0, 'speed': 0.0, 'present_TCR': 5e4, 'present_PD1': 0.0,
              'transfer_cytotoxic': 40.0},
    }
    out = ph.update({'cells': cells}, 60.0)['cells']
    # T-cell should have picked up the tumor's presented MHCI
    assert out['t']['accept_MHCI'] == pytest.approx(5e4)
    # tumor should receive the T-cell's transferred packets
    assert out['u']['receive_cytotoxic'] == pytest.approx(40.0)
