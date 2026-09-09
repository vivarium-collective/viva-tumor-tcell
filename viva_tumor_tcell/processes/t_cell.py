"""T cell process — faithful port of tumor_tcell/processes/t_cell.py.

Two states: ``PD1n`` (active — high IFNg + cytotoxic-packet secretion, fast
migration) and ``PD1p`` (exhausted — ~10x less secretion, slower, dies faster
next to PDL1+ tumor). TCR down-regulates after 6 h of activation, then a
refractory period; repeated refractory cycles or division count drive the
PD1n -> PD1p transition.

Flat-agent port: velocities are stored as a scalar ``speed`` (µm/s); the physics
process gives that speed a fresh random direction each tick (tumor-tcell's
persistent-random-walk migration). Equations/parameters transcribed verbatim
(µm/min migration params converted to µm/s).
"""
import math
import random

from process_bigraph import Process

NAME = 'TCellProcess'
TIMESTEP = 60
_PER_MIN_TO_PER_S = 1.0 / 60.0


def get_probability_timestep(probability_parameter, timescale, timestep):
    rate = -math.log(1 - probability_parameter)
    return 1 - math.exp(-rate * (timestep / timescale))


def _daughter_locations(loc, diameter):
    mx, my = loc
    return [(mx + random.gauss(0, 0.1) * diameter, my + random.gauss(0, 0.1) * diameter)
            for _ in range(2)]


class TCellProcess(Process):
    config_schema = {
        'agents_key': {'_type': 'string', '_default': 'cells'},
        'diameter': {'_type': 'float', '_default': 7.5},           # µm
        'mass': {'_type': 'float', '_default': 2.0},               # ng
        'initial_PD1n': {'_type': 'float', '_default': 0.8},
        'refractory_count_threshold': {'_type': 'float', '_default': 3},
        'activation_time': {'_type': 'float', '_default': 21600},
        'activation_refractory_time': {'_type': 'float', '_default': 43200},
        'TCR_downregulated': {'_type': 'float', '_default': 0.0},
        'TCR_upregulated': {'_type': 'float', '_default': 50000.0},
        'PDL1_critical_number': {'_type': 'float', '_default': 1e4},
        'death_PD1p_14hr': {'_type': 'float', '_default': 0.35},
        'death_PD1n_14hr': {'_type': 'float', '_default': 0.1},
        'death_PD1p_next_to_PDL1p_14hr': {'_type': 'float', '_default': 0.475},
        'PD1n_IFNg_production': {'_type': 'float', '_default': 1.62e4 / 3600},
        'PD1p_IFNg_production': {'_type': 'float', '_default': 1.62e3 / 3600},
        'PD1p_PD1_equilibrium': {'_type': 'float', '_default': 5e4},
        'ligand_threshold': {'_type': 'float', '_default': 1e4},
        'PD1n_growth_28hr': {'_type': 'float', '_default': 0.90},
        'PD1p_growth_28hr': {'_type': 'float', '_default': 0.20},
        'PD1n_divide_threshold': {'_type': 'float', '_default': 5},
        'PD1n_migration': {'_type': 'float', '_default': 10.0 * _PER_MIN_TO_PER_S},   # µm/s
        'migration_MHCIp_tumor_dwell_velocity': {'_type': 'float', '_default': 0.0},
        'PD1n_migration_MHCIp_tumor_dwell_time': {'_type': 'float', '_default': 25.0 * 60},
        'PD1p_migration': {'_type': 'float', '_default': 5.0 * _PER_MIN_TO_PER_S},
        'PD1p_migration_MHCIp_tumor_dwell_time': {'_type': 'float', '_default': 10.0 * 60},
        'PD1n_migration_refractory_time': {'_type': 'float', '_default': 35.0 * 60},
        'PD1p_migration_refractory_time': {'_type': 'float', '_default': 20.0 * 60},
        'cytotoxic_packet_production': {'_type': 'float', '_default': 40 / 60},
        'PD1n_cytotoxic_packets_max': {'_type': 'float', '_default': 10000.0},
        'PD1p_cytotoxic_packets_max': {'_type': 'float', '_default': 1000.0},
        'MHCIn_reduction_production': {'_type': 'float', '_default': 400.0},
        'cytotoxic_transfer_rate': {'_type': 'float', '_default': 400.0},
    }

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def outputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def update(self, state, interval):
        """Top-level process: step every T cell in the shared map."""
        cells = state['cells']
        out = {}
        add, remove = {}, []
        for agent_id, agent in cells.items():
            if agent.get('cell_type') != 't-cell' or agent.get('death'):
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

        cell_state = agent.get('cell_state') or 'PD1n'
        PDL1 = float(agent.get('accept_PDL1', 0.0) or 0.0)
        MHCI = float(agent.get('accept_MHCI', 0.0) or 0.0)
        TCR_timer = float(agent.get('TCR_timer', 0.0) or 0.0)
        velocity_timer = float(agent.get('velocity_timer', 0.0) or 0.0)
        TCR = float(agent.get('present_TCR', 0.0) or 0.0)
        refractory_count = float(agent.get('refractory_count', 0.0) or 0.0)
        PD1n_divide_counts = float(agent.get('PD1n_divide_count', 0.0) or 0.0)
        stockpile = float(agent.get('total_cytotoxic_packets', 0.0) or 0.0)

        # ---- death ----
        if cell_state == 'PD1n':
            if random.uniform(0, 1) < get_probability_timestep(p['death_PD1n_14hr'], 50400, timestep):
                return {agent_id: {'death': 'PD1n_apoptosis'}}
        elif cell_state == 'PD1p':
            if PDL1 >= p['PDL1_critical_number']:
                if random.uniform(0, 1) < get_probability_timestep(p['death_PD1p_next_to_PDL1p_14hr'], 50400, timestep):
                    return {agent_id: {'death': 'PD1p_PDL1_death'}}
            else:
                if random.uniform(0, 1) < get_probability_timestep(p['death_PD1p_14hr'], 50400, timestep):
                    return {agent_id: {'death': 'PD1p_apoptosis'}}

        # ---- division ----
        if cell_state == 'PD1n':
            if random.uniform(0, 1) < get_probability_timestep(p['PD1n_growth_28hr'], 100800, timestep):
                return self._divide(agent_id, agent, count_daughter=False)
        elif cell_state == 'PD1p':
            if random.uniform(0, 1) < get_probability_timestep(p['PD1p_growth_28hr'], 100800, timestep):
                return self._divide(agent_id, agent, count_daughter='PD1p')

        u = {}

        # TCR down-regulation after 6 h of activation
        if TCR_timer > p['activation_time']:
            TCR = p['TCR_downregulated']
            u['present_TCR'] = TCR

        # advance TCR timer while in contact (or during refractory period)
        if MHCI > p['ligand_threshold']:
            u['TCR_timer'] = timestep
        elif TCR_timer > p['activation_time']:
            u['TCR_timer'] = timestep

        # after refractory period: restore TCR, clear stockpile, count a cycle
        if TCR_timer > p['activation_refractory_time']:
            TCR = p['TCR_upregulated']
            u['TCR_timer'] = -p['activation_refractory_time']
            u['total_cytotoxic_packets'] = -stockpile
            u['present_TCR'] = TCR
            u['refractory_count'] = 1
            stockpile = 0.0

        # ---- state transition ----
        new_cell_state = cell_state
        if cell_state == 'PD1n':
            if (refractory_count > p['refractory_count_threshold']
                    or PD1n_divide_counts > p['PD1n_divide_threshold']):
                new_cell_state = 'PD1p'
                u['cell_state'] = 'PD1p'
                u['cell_state_count'] = 1

        # ---- behavior ----
        if new_cell_state == 'PD1n':
            packets_max = p['PD1n_cytotoxic_packets_max']
            ifng_prod = p['PD1n_IFNg_production']
            migration = p['PD1n_migration']
            dwell_time = p['PD1n_migration_MHCIp_tumor_dwell_time']
            refractory_time = p['PD1n_migration_refractory_time']
        else:  # PD1p
            u['present_PD1'] = p['PD1p_PD1_equilibrium']
            packets_max = p['PD1p_cytotoxic_packets_max']
            ifng_prod = p['PD1p_IFNg_production']
            migration = p['PD1p_migration']
            dwell_time = p['PD1p_migration_MHCIp_tumor_dwell_time']
            refractory_time = p['PD1p_migration_refractory_time']

        if MHCI > 0:
            u['velocity_timer'] = u.get('velocity_timer', 0.0) + timestep

        produced = None
        if MHCI >= p['ligand_threshold'] and TCR >= p['ligand_threshold']:
            if stockpile < packets_max:
                u['total_cytotoxic_packets'] = u.get('total_cytotoxic_packets', 0.0) \
                    + p['cytotoxic_packet_production'] * timestep
            u.setdefault('exchange', {})['IFNg'] = int(ifng_prod * timestep)
            produced = True
        elif MHCI > 0 and TCR >= p['ligand_threshold']:
            # 4-fold reduction vs MHCI-low tumor
            if stockpile < packets_max:
                u['total_cytotoxic_packets'] = u.get('total_cytotoxic_packets', 0.0) \
                    + p['cytotoxic_packet_production'] / p['MHCIn_reduction_production'] * timestep
            u.setdefault('exchange', {})['IFNg'] = int(ifng_prod / p['MHCIn_reduction_production'] * timestep)
            produced = True

        if produced is None:
            # migration velocity management (dwell/refractory)
            if velocity_timer >= refractory_time:
                u['velocity_timer'] = -velocity_timer  # reset to 0
                u['speed'] = migration
            elif velocity_timer > dwell_time:
                u['speed'] = migration
            elif 0 < velocity_timer < dwell_time:
                u['speed'] = p['migration_MHCIp_tumor_dwell_velocity']
            elif velocity_timer == 0:
                u['speed'] = migration

        # ---- cytotoxic transfer to the contacted tumor ----
        stockpile_now = stockpile + u.get('total_cytotoxic_packets', 0.0)
        cytotoxic_transfer = min(p['cytotoxic_transfer_rate'], max(0.0, stockpile_now))
        if cytotoxic_transfer > 0:
            # NOTE: original t_cell.py wrote key `cytotoxic_packets` (a schema
            # typo) so the stockpile was never actually drawn down. We port the
            # intended semantics and decrement `total_cytotoxic_packets`.
            u['total_cytotoxic_packets'] = u.get('total_cytotoxic_packets', 0.0) - cytotoxic_transfer
        u['transfer_cytotoxic'] = cytotoxic_transfer

        return {agent_id: u}

    def _divide(self, agent_id, agent, count_daughter=False):
        diameter = float(agent.get('radius', self.config['diameter'] / 2.0)) * 2.0
        locs = _daughter_locations(agent.get('location', (0.0, 0.0)), diameter)

        cell_state = agent.get('cell_state', 'PD1n')
        stockpile = float(agent.get('total_cytotoxic_packets', 0.0) or 0.0)
        pd1n_count = float(agent.get('PD1n_divide_count', 0.0) or 0.0)

        base = {
            'type': 'circle', 'cell_type': 't-cell', 'cell_state': cell_state,
            'mass': float(agent.get('mass', self.config['mass'])),
            'radius': diameter / 2.0, 'velocity': (0.0, 0.0),
            'speed': float(agent.get('speed', self.config['PD1n_migration']) or 0.0),
            'present_TCR': 50000.0,
            'total_cytotoxic_packets': stockpile / 2.0,  # split
        }
        add = {}
        for i, (suffix, dloc) in enumerate(zip(('A', 'B'), locs)):
            did = f"{agent_id}{suffix}"
            d = dict(base)
            d['id'] = did
            d['location'] = (float(dloc[0]), float(dloc[1]))
            # asymmetric_division: one daughter's PD1n count increments
            if cell_state == 'PD1n':
                d['PD1n_divide_count'] = pd1n_count + 1 if i == 0 else pd1n_count
            if count_daughter == 'PD1p' and i == 0:
                d['PD1p_divide_count'] = float(agent.get('PD1p_divide_count', 0.0) or 0.0) + 1
            add[did] = d
        return {'_add': add, '_remove': [agent_id]}


# --- workbench viewer contract (per-port meanings + units) ---
TCellProcess.contract = {
    'summary': "CD8+ T cell behavior — active PD1- <-> exhausted PD1+. In contact with an MHCI+ tumor "
               "it secretes IFNg (which converts tumors) and cytotoxic packets (which kill); TCR "
               "down-regulates after ~6 h of activation then a refractory period; repeated cycles or "
               "division count drive PD1- -> PD1+ exhaustion.",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); per T cell reads cell_state, accept_MHCI / "
                        "accept_PDL1 (from the contacted tumor), present_TCR, TCR + velocity timers (s), "
                        "refractory / divide counts, and cytotoxic stockpile."},
    'outputs': {'cells': "Per T cell: IFNg secretion (field exchange, counts) and transfer_cytotoxic "
                         "(packets) in contact; TCR up/down-regulation; PD1- -> PD1+ exhaustion; migration "
                         "speed (um/s, with dwell logic); division into two daughters (_add/_remove)."},
    'assumptions': ["PD1- secretes ~10x more IFNg/cytotoxic than PD1+; production is 4-fold reduced vs an "
                    "MHCI-low tumor; migration 10 (PD1-) / 5 (PD1+) um/min; diameter 7.5 um, mass 2 ng. "
                    "Parameters from tumor-tcell."],
}
