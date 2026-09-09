"""Dendritic cell process — faithful port of tumor_tcell/processes/dendritic_cell.py.

A dendritic cell takes up ``tumor_debris`` from the field; once its internalized
debris exceeds a threshold it becomes ``active`` and presents MHCI/PDL1 (used in
the lymph-node variant to prime T cells). Dies by apoptosis; divides only when
active.

Flat-agent port (units stripped): debris field concentration in ng/mL, counts as
counts. Diffusion/MW constants transcribed from the source.
"""
import math
import random

from process_bigraph import Process
from scipy import constants

NAME = 'DendriticCellProcess'
TIMESTEP = 60
AVOGADRO = constants.N_A
PI = constants.pi

# tumor_debris diffusion: 0.0864 cm^2/day -> µm^2/s
DEBRIS_DIFFUSION_UM2_PER_S = 0.0864 * 1e8 / 86400.0
DEBRIS_MW = 29000.0
DEBRIS_COUNTS_PER_NG = AVOGADRO / DEBRIS_MW / 1e9


def get_probability_timestep(probability_parameter, timescale, timestep):
    rate = -math.log(1 - probability_parameter)
    return 1 - math.exp(-rate * (timestep / timescale))


def _daughter_locations(loc, diameter):
    mx, my = loc
    return [(mx + random.gauss(0, 0.1) * diameter, my + random.gauss(0, 0.1) * diameter)
            for _ in range(2)]


class DendriticCellProcess(Process):
    config_schema = {
        'agents_key': {'_type': 'string', '_default': 'cells'},
        'diameter': {'_type': 'float', '_default': 10.0},          # µm
        'mass': {'_type': 'float', '_default': 2.0},               # ng
        'velocity': {'_type': 'float', '_default': 3.0 / 60.0},    # µm/s (3 µm/min)
        'death_apoptosis': {'_type': 'float', '_default': 0.5},
        'death_time': {'_type': 'float', '_default': 4 * 24 * 60 * 60},
        'divide_prob': {'_type': 'float', '_default': 0.6},
        'divide_time': {'_type': 'float', '_default': 5 * 24 * 60 * 60},
        'internal_tumor_debris_threshold': {'_type': 'float', '_default': 415000.0},
        'PDL1p_PDL1_equilibrium': {'_type': 'float', '_default': 5e4},
        'PDL1p_MHCI_equilibrium': {'_type': 'float', '_default': 5e4},
        'tumor_debris_uptake': {'_type': 'float', '_default': 300 / 60},
    }

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def outputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def update(self, state, interval):
        """Top-level process: step every dendritic cell in the shared map."""
        cells = state['cells']
        out = {}
        add, remove = {}, []
        for agent_id, agent in cells.items():
            if agent.get('cell_type') != 'dendritic' or agent.get('death'):
                continue
            res = self._step_agent(agent_id, agent, interval)
            if not res:
                continue
            if '_add' in res:
                add.update(res.pop('_add'))
            if '_remove' in res:
                remove += res.pop('_remove')
            out.update(res)
        if add:
            out['_add'] = add
        if remove:
            out['_remove'] = remove
        return {'cells': out}

    def _step_agent(self, agent_id, agent, interval):
        timestep = float(interval)
        p = self.config

        cell_state = agent.get('cell_state') or 'inactive'
        external_debris = float((agent.get('local') or {}).get('tumor_debris', 0.0) or 0.0)  # ng/mL
        internal_debris = float(agent.get('internal_tumor_debris', 0.0) or 0.0)
        diameter = float(agent.get('radius', p['diameter'] / 2.0)) * 2.0

        # available debris counts within the diffusion sphere this interval
        diffusion_radius = (DEBRIS_DIFFUSION_UM2_PER_S * timestep) ** 0.5
        sphere_radius = diameter / 2 + diffusion_radius
        avail_volume = 4 / 3 * PI * sphere_radius ** 3
        available_debris = external_debris * avail_volume * DEBRIS_COUNTS_PER_NG / 1e12
        if not math.isfinite(available_debris):
            available_debris = 0.0

        # death by apoptosis
        if random.uniform(0, 1) < get_probability_timestep(p['death_apoptosis'], p['death_time'], timestep):
            return {agent_id: {'death': 'apoptosis'}}

        # division (active only)
        if cell_state == 'active':
            if random.uniform(0, 1) < get_probability_timestep(p['divide_prob'], p['divide_time'], timestep):
                return self._divide(agent_id, agent)

        u = {}
        new_cell_state = cell_state
        if cell_state == 'inactive' and internal_debris >= p['internal_tumor_debris_threshold']:
            new_cell_state = 'active'
            u['cell_state'] = 'active'
            u['cell_state_count'] = 1

        # uptake locally available debris
        uptake = min(int(p['tumor_debris_uptake'] * timestep), int(available_debris))
        u['exchange'] = {'tumor_debris': -uptake}
        u['internal_tumor_debris'] = uptake

        if new_cell_state == 'active':
            u['present_PDL1'] = p['PDL1p_PDL1_equilibrium']
            u['present_MHCI'] = p['PDL1p_MHCI_equilibrium']

        return {agent_id: u}

    def _divide(self, agent_id, agent):
        diameter = float(agent.get('radius', self.config['diameter'] / 2.0)) * 2.0
        locs = _daughter_locations(agent.get('location', (0.0, 0.0)), diameter)
        base = {
            'type': 'circle', 'cell_type': 'dendritic', 'cell_state': agent.get('cell_state', 'active'),
            'mass': float(agent.get('mass', self.config['mass'])), 'radius': diameter / 2.0,
            'speed': float(agent.get('speed', self.config['velocity']) or 0.0),
            'present_MHCI': float(agent.get('present_MHCI', 0.0) or 0.0),
            'present_PDL1': float(agent.get('present_PDL1', 0.0) or 0.0),
            'internal_tumor_debris': float(agent.get('internal_tumor_debris', 0.0) or 0.0) / 2.0,
        }
        add = {}
        for suffix, dloc in zip(('A', 'B'), locs):
            did = f"{agent_id}{suffix}"
            d = dict(base); d['id'] = did
            d['location'] = (float(dloc[0]), float(dloc[1]))
            add[did] = d
        return {'_add': add, '_remove': [agent_id]}


# --- workbench viewer contract (per-port meanings + units) ---
DendriticCellProcess.contract = {
    'summary': "Dendritic cell — takes up tumor_debris from the field; once internalized debris exceeds "
               "415000 counts it activates and presents MHCI/PDL1 (the antigen-presentation arm of the "
               "lymph-node extension).",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); per dendritic cell reads cell_state "
                        "(inactive/active), internal_tumor_debris (counts), and local tumor_debris "
                        "concentration (ng/mL)."},
    'outputs': {'cells': "Per dendritic cell: debris uptake (internal_tumor_debris + / field exchange -); "
                         "inactive -> active at threshold; present_MHCI/present_PDL1 (5e4) when active; "
                         "apoptosis; division when active (_add/_remove)."},
    'assumptions': ["Only active DCs divide; uptake ~300 molecules/cell/min; activation threshold 415000 "
                    "counts; diameter 10 um, mass 2 ng."],
}
