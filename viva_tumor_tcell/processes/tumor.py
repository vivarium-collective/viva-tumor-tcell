"""Tumor cell process — faithful port of tumor_tcell/processes/tumor.py.

Two states: ``PDL1n`` (proliferative, low MHCI/PDL1) and ``PDL1p`` (quiescent,
high MHCI + PDL1). Internalizing enough IFNg flips PDL1n -> PDL1p. Dies by
apoptosis or when accumulated cytotoxic packets exceed a threshold; on death it
releases tumor_debris into its exchange.

Flat-agent port: reads the whole ``cells`` map, acts on ``agent_id``, returns
delta updates plus ``_add``/``_remove`` for division. Equations and parameters
are transcribed verbatim from the vivarium-1.0 source; only units are stripped
(length µm, time s, mass ng, IFNg field concentration ng/mL, counts as counts).
"""
import math
import random

from process_bigraph import Process
from scipy import constants

NAME = 'TumorCellProcess'
TIMESTEP = 60
AVOGADRO = constants.N_A
PI = constants.pi

# IFNg diffusion: 1.25e-3 cm^2/day -> µm^2/s  (1 cm^2 = 1e8 µm^2, 1 day = 86400 s)
IFNG_DIFFUSION_UM2_PER_S = 1.25e-3 * 1e8 / 86400.0
# counts per ng of IFNg: N_A / MW[g/mol] / 1e9(ng per g)
IFNG_MW = 17000.0
IFNG_COUNTS_PER_NG = AVOGADRO / IFNG_MW / 1e9


def get_probability_timestep(probability_parameter, timescale, timestep):
    """Transition probability as a function of time (verbatim)."""
    rate = -math.log(1 - probability_parameter)
    timestep_fraction = timestep / timescale
    return 1 - math.exp(-rate * timestep_fraction)


def _daughter_locations(loc, diameter):
    """Two daughters placed near the mother with Gaussian jitter (verbatim intent)."""
    mx, my = loc
    out = []
    for _ in range(2):
        out.append((mx + random.gauss(0, 0.1) * diameter,
                    my + random.gauss(0, 0.1) * diameter))
    return out


class TumorCellProcess(Process):
    config_schema = {
        'agents_key': {'_type': 'string', '_default': 'cells'},
        'diameter': {'_type': 'float', '_default': 15.0},          # µm
        'mass': {'_type': 'float', '_default': 8.0},               # ng
        'initial_PDL1n': {'_type': 'float', '_default': 0.9},
        'death_apoptosis': {'_type': 'float', '_default': 0.5},
        'cytotoxic_packet_threshold': {'_type': 'float', '_default': 128 * 100},
        'PDL1n_growth': {'_type': 'float', '_default': 0.6},
        'Max_IFNg_internalization': {'_type': 'float', '_default': 31 / 60},
        'IFNg_threshold': {'_type': 'float', '_default': 15000.0},
        'reduction_IFNg_internalization': {'_type': 'float', '_default': 2.0},
        'PDL1p_PDL1_equilibrium': {'_type': 'float', '_default': 5e4},
        'PDL1p_MHCI_equilibrium': {'_type': 'float', '_default': 5e4},
        'tumor_debris_amount': {'_type': 'float', '_default': 1.4e15},
    }

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def outputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def update(self, state, interval):
        """Top-level process: step every tumor cell in the shared map, merging
        per-cell deltas plus any division _add/_remove into one update."""
        cells = state['cells']
        out = {}
        add, remove = {}, []
        for agent_id, agent in cells.items():
            if agent.get('cell_type') != 'tumor' or agent.get('death'):
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
        """Per-cell biology (verbatim). Returns {agent_id: delta} and/or
        {'_add': {...}, '_remove': [...]} for a division, or {} for no-op."""
        timestep = float(interval)
        p = self.config

        cell_state = agent.get('cell_state') or 'PDL1n'
        cytotoxic_packets = float(agent.get('receive_cytotoxic', 0.0) or 0.0)
        external_IFNg = float((agent.get('local') or {}).get('IFNg', 0.0) or 0.0)  # ng/mL
        internal_IFNg = float(agent.get('internal_IFNg', 0.0) or 0.0)
        diameter = float(agent.get('radius', p['diameter'] / 2.0)) * 2.0

        # available IFNg counts within the diffusion sphere this interval
        diffusion_area = IFNG_DIFFUSION_UM2_PER_S * timestep
        diffusion_radius = diffusion_area ** 0.5
        sphere_radius = diameter / 2 + diffusion_radius
        avail_volume = 4 / 3 * PI * sphere_radius ** 3            # µm^3
        available_IFNg = external_IFNg * avail_volume * IFNG_COUNTS_PER_NG / 1e12
        if not math.isfinite(available_IFNg):
            available_IFNg = 0.0

        agent_update = {}

        # death by apoptosis (0.95 by 5 days -> death_apoptosis over 432000 s)
        prob_death = get_probability_timestep(p['death_apoptosis'], 432000, timestep)
        if random.uniform(0, 1) < prob_death:
            return {agent_id: {
                'exchange': {'tumor_debris': int(p['tumor_debris_amount'])},
                'death': 'apoptosis'}}

        # death by cytotoxic packets
        if cytotoxic_packets >= p['cytotoxic_packet_threshold']:
            return {agent_id: {
                'exchange': {'tumor_debris': int(p['tumor_debris_amount'])},
                'death': 'Tcell_death'}}

        # division (PDL1n only; PDL1p is arrested)
        if cell_state == 'PDL1n':
            prob_divide = get_probability_timestep(p['PDL1n_growth'], 86400, timestep)
            if random.uniform(0, 1) < prob_divide:
                return self._divide(agent_id, agent)

        # state transition PDL1n -> PDL1p driven by internalized IFNg
        new_cell_state = cell_state
        if cell_state == 'PDL1n' and internal_IFNg >= p['IFNg_threshold']:
            new_cell_state = 'PDL1p'
            agent_update['cell_state'] = 'PDL1p'
            agent_update['cell_state_count'] = 1

        # behavior: present ligands + internalize/degrade IFNg
        if new_cell_state == 'PDL1p':
            agent_update['present_PDL1'] = p['PDL1p_PDL1_equilibrium']
            agent_update['present_MHCI'] = p['PDL1p_MHCI_equilibrium']
            IFNg_degrade = min(
                int(p['Max_IFNg_internalization'] / p['reduction_IFNg_internalization'] * timestep),
                int(available_IFNg))
        else:  # PDL1n
            IFNg_degrade = min(
                int(p['Max_IFNg_internalization'] * timestep),
                int(available_IFNg))

        agent_update.setdefault('exchange', {})['IFNg'] = -IFNg_degrade
        agent_update['internal_IFNg'] = IFNg_degrade

        return {agent_id: agent_update}

    def _divide(self, agent_id, agent):
        diameter = float(agent.get('radius', self.config['diameter'] / 2.0)) * 2.0
        loc = agent.get('location', (0.0, 0.0))
        locs = _daughter_locations(loc, diameter)

        base_shared = {
            'type': 'circle',
            'cell_type': 'tumor',
            'cell_state': agent.get('cell_state', 'PDL1n'),
            'mass': float(agent.get('mass', self.config['mass'])),
            'radius': diameter / 2.0,
            'velocity': (0.0, 0.0),
            'speed': 0.0,
            'present_MHCI': float(agent.get('present_MHCI', 1000.0) or 1000.0),
            'present_PDL1': float(agent.get('present_PDL1', 0.0) or 0.0),
            # split accumulators
            'internal_IFNg': float(agent.get('internal_IFNg', 0.0) or 0.0) / 2.0,
            'receive_cytotoxic': float(agent.get('receive_cytotoxic', 0.0) or 0.0) / 2.0,
        }
        add = {}
        for suffix, dloc in zip(('A', 'B'), locs):
            did = f"{agent_id}{suffix}"
            d = dict(base_shared)
            d['id'] = did
            d['location'] = (float(dloc[0]), float(dloc[1]))
            add[did] = d
        return {'_add': add, '_remove': [agent_id]}
