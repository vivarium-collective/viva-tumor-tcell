"""TumorTcellPhysics — collisions via viva-munk + neighbor-exchange.

Replaces tumor-tcell's ``Neighbors`` process. The collision physics is
delegated to a real ``viva_munk.processes.multibody.PymunkProcess`` instance
(walls, jitter, substeps, circle-circle collisions); this process:

  1. gives each cell's stored ``speed`` a fresh random direction (tumor-tcell's
     persistent-random-walk migration) and feeds flat bodies to viva-munk,
  2. steps the physics for the interval and reads back position deltas,
  3. computes neighbors from the new positions (T-cell picks its single nearest
     tumor; a tumor collects all its T-cell neighbors), and
  4. performs the ligand exchange: each cell's ``accept_*`` <- neighbours'
     ``present_*`` (membrane), and a T-cell's ``transfer_cytotoxic`` -> the
     contacted tumor's ``receive_cytotoxic`` (soluble).
"""
import math
import random

from process_bigraph import Process
from viva_munk.processes.multibody import PymunkProcess


class TumorTcellPhysics(Process):
    # scalar config only (see DiffusionField note on list-config concatenation)
    config_schema = {
        'bounds_x': {'_type': 'float', '_default': 200.0},   # µm
        'bounds_y': {'_type': 'float', '_default': 200.0},   # µm
        'neighbor_distance': {'_type': 'float', '_default': 1.0},   # µm from membranes
        'jitter_per_second': {'_type': 'float', '_default': 1e-3},
        'substeps': {'_type': 'integer', '_default': 10},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self.physics = PymunkProcess(config={
            'env_size': float(self.config['bounds_x']),
            'env_height': float(self.config['bounds_y']),
            'substeps': int(self.config['substeps']),
            'damping_per_second': 1.0,   # ballistic over a tick (no velocity decay)
            'jitter_per_second': float(self.config['jitter_per_second']),
            'friction': 0.9,
            'elasticity': 0.0,
            'wall_thickness': 10.0,
        }, core=core)
        self.neighbor_distance = float(self.config['neighbor_distance'])

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def outputs(self):
        return {'cells': 'map[tumor_tcell_agent]'}

    def update(self, state, interval):
        cells = state['cells']
        if not cells:
            return {'cells': {}}

        # --- build flat bodies with a fresh random-direction velocity ---
        bodies = {}
        old_loc = {}
        for cid, c in cells.items():
            loc = tuple(c.get('location', (0.0, 0.0)))
            old_loc[cid] = loc
            speed = float(c.get('speed', 0.0) or 0.0)
            theta = random.uniform(0, 2 * math.pi)
            vel = (speed * math.cos(theta), speed * math.sin(theta))
            bodies[cid] = {
                'type': 'circle',
                'location': loc,
                'velocity': vel,
                'radius': float(c.get('radius', 3.75) or 3.75),
                'mass': float(c.get('mass', 1.0) or 1.0),
            }

        phys_out = self.physics.update({'circle_particles': bodies}, interval)
        deltas = phys_out.get('circle_particles', {})

        # new absolute positions
        new_loc = {}
        for cid, loc in old_loc.items():
            d = deltas.get(cid, {}).get('location', (0.0, 0.0))
            new_loc[cid] = (loc[0] + d[0], loc[1] + d[1])

        # --- neighbor detection ---
        radii = {cid: float(cells[cid].get('radius', 3.75) or 3.75) for cid in cells}
        tcells = {cid: new_loc[cid] for cid, c in cells.items() if c.get('cell_type') == 't-cell'}
        tumors = {cid: new_loc[cid] for cid, c in cells.items() if c.get('cell_type') == 'tumor'}

        def within(a_loc, a_r, b_loc, b_r):
            dist = math.hypot(a_loc[0] - b_loc[0], a_loc[1] - b_loc[1])
            return (dist - a_r - b_r) <= self.neighbor_distance

        neighbors = {cid: [] for cid in cells}
        # each T-cell polarizes to its single nearest tumor
        for tc, tloc in tcells.items():
            best, best_d = None, None
            for tu, uloc in tumors.items():
                if within(tloc, radii[tc], uloc, radii[tu]):
                    d = math.hypot(tloc[0] - uloc[0], tloc[1] - uloc[1])
                    if best_d is None or d < best_d:
                        best, best_d = tu, d
            if best is not None:
                neighbors[tc] = [best]
        # each tumor collects all its T-cell neighbors
        for tu, uloc in tumors.items():
            neighbors[tu] = [tc for tc, tloc in tcells.items()
                             if within(uloc, radii[tu], tloc, radii[tc])]

        # --- build updates: movement + ligand exchange ---
        update = {}
        for cid, c in cells.items():
            u = {}
            d = deltas.get(cid, {}).get('location', (0.0, 0.0))
            if d != (0.0, 0.0):
                u['location'] = (d[0], d[1])
            ns = neighbors[cid]
            ctype = c.get('cell_type')
            if ctype == 't-cell':
                if ns:
                    tu = cells[ns[0]]
                    u['accept_PDL1'] = float(tu.get('present_PDL1', 0.0) or 0.0)
                    u['accept_MHCI'] = float(tu.get('present_MHCI', 0.0) or 0.0)
                else:
                    u['accept_PDL1'] = 0.0
                    u['accept_MHCI'] = 0.0
            elif ctype == 'tumor':
                acc_TCR = sum(float(cells[n].get('present_TCR', 0.0) or 0.0) for n in ns)
                acc_PD1 = sum(float(cells[n].get('present_PD1', 0.0) or 0.0) for n in ns)
                u['accept_TCR'] = acc_TCR
                u['accept_PD1'] = acc_PD1
                # receive cytotoxic packets from all contacting T-cells this tick
                received = sum(float(cells[n].get('transfer_cytotoxic', 0.0) or 0.0) for n in ns)
                if received:
                    u['receive_cytotoxic'] = received
            if u:
                update[cid] = u

        return {'cells': update}
