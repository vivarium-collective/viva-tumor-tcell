"""DiffusionField — faithful consolidation of tumor-tcell's Fields + LocalField.

Each tick, for the 2D concentration fields (ng/mL) over the chamber:
  1. deposit every cell's ``exchange`` counts into its field bin
     (count -> concentration via Avogadro / MW / bin-volume),
  2. diffuse (explicit Laplacian) and exponentially decay each field,
  3. sample each cell's local concentration back into ``local``,
  4. release debris + remove cells flagged dead (their exchange, incl. debris,
     is deposited this same tick before removal).

Diffusion/decay constants are transcribed from tumor_tcell/processes/fields.py:
IFNg D = 1.25e-3 cm²/day, tumor_debris D = 0.0864 cm²/day; IFNg decay = 4.5 h
e-fold. MWs: IFNg 17000, tumor_debris 29000 g/mol.
"""
import numpy as np
from process_bigraph import Process
from scipy import constants

AVOGADRO = constants.N_A
CM2_PER_DAY_TO_UM2_PER_S = 1e8 / 86400.0

DIFFUSION_UM2_PER_S = {
    'IFNg': 1.25e-3 * CM2_PER_DAY_TO_UM2_PER_S,
    'tumor_debris': 0.0864 * CM2_PER_DAY_TO_UM2_PER_S,
}
DECAY_PER_S = {
    'IFNg': np.log(2) / (4.5 * 60 * 60),   # 4.5 h e-fold
}
MOLECULAR_WEIGHT = {'IFNg': 17000.0, 'tumor_debris': 29000.0}  # g/mol

# Concentration cap (ng/mL): tumor death dumps 1.4e15 debris molecules -> a huge
# point-source concentration that can overflow to inf/NaN under diffusion. Uptake
# is rate-capped anyway, so clamp the field to a large finite value.
_FIELD_CAP = 1e12

# _LAP: 5-point Laplacian via np.roll (reflect-free, adequate for the
# qualitative diffusion the studies exercise).
def _laplacian(f):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0)
            + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4.0 * f)


class DiffusionField(Process):
    description = (
        "Soluble-field manager (tumor-tcell Fields + LocalField consolidated): "
        "deposits each cell's secreted amounts into its field bin, diffuses and "
        "decays the 2D IFNg / tumor-debris concentration grids (ng/mL), and "
        "samples each cell's local concentration back onto it."
    )
    # NB: scalar config only — a `list`-typed config field concatenates the
    # default with the provided value (list apply is additive), so bounds/n_bins
    # are passed as x/y scalars. `molecules` defaults to [] (empty is
    # concat-safe: [] + ['IFNg'] == ['IFNg']).
    config_schema = {
        'bounds_x': {'_type': 'float', '_default': 200.0},   # µm
        'bounds_y': {'_type': 'float', '_default': 200.0},   # µm
        'n_bins_x': {'_type': 'integer', '_default': 20},
        'n_bins_y': {'_type': 'integer', '_default': 20},
        'depth': {'_type': 'float', '_default': 15.0},       # µm
        'molecules': {'_type': 'list[string]', '_default': []},
        'diffusion_dt': {'_type': 'float', '_default': 1.0}, # s
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self.bounds = [float(self.config['bounds_x']), float(self.config['bounds_y'])]
        self.n_bins = [int(self.config['n_bins_x']), int(self.config['n_bins_y'])]
        self.depth = float(self.config['depth'])
        self.molecules = list(self.config['molecules']) or ['IFNg']
        self.dx = self.bounds[0] / self.n_bins[0]
        self.dy = self.bounds[1] / self.n_bins[1]
        # bin volume in liters (1 L = 1e15 µm³)
        self.bin_volume_L = (self.dx * self.dy * self.depth) * 1e-15
        # ng/mL per (mol/L) = MW[g/mol] * 1e6
        self.conc_conversion = {m: MOLECULAR_WEIGHT.get(m, 1.0) * 1e6 for m in self.molecules}

    def inputs(self):
        arr = {'_type': 'array', '_shape': tuple(self.n_bins), '_data': 'float'}
        return {'cells': 'map[tumor_tcell_agent]', 'fields': {'_type': 'map', '_value': arr}}

    def outputs(self):
        arr = {'_type': 'array', '_shape': tuple(self.n_bins), '_data': 'float'}
        return {'cells': 'map[tumor_tcell_agent]', 'fields': {'_type': 'map', '_value': arr}}

    def _bin(self, location):
        ix = int(location[0] / self.bounds[0] * self.n_bins[0])
        iy = int(location[1] / self.bounds[1] * self.n_bins[1])
        ix = min(max(ix, 0), self.n_bins[0] - 1)
        iy = min(max(iy, 0), self.n_bins[1] - 1)
        return ix, iy

    def update(self, state, interval):
        timestep = float(interval)
        cells = state['cells']
        original = {m: np.array(state['fields'][m], dtype=float) for m in self.molecules}
        fields = {m: original[m].copy() for m in self.molecules}

        # 1. deposit exchange counts -> concentration into the cell's bin
        for cid, c in cells.items():
            exch = c.get('exchange') or {}
            if not exch:
                continue
            ix, iy = self._bin(c.get('location', (0.0, 0.0)))
            for mol, counts in exch.items():
                if mol not in fields or not counts:
                    continue
                conc = counts / (self.bin_volume_L * AVOGADRO) * self.conc_conversion[mol]
                fields[mol][ix, iy] += conc

        # 2. diffuse + decay
        for mol in self.molecules:
            f = fields[mol]
            D = DIFFUSION_UM2_PER_S.get(mol, 1.0) / (self.dx * self.dy)
            dt = min(self.config['diffusion_dt'], timestep)
            t = 0.0
            if len(np.unique(f)) != 1:  # skip if uniform
                while t < timestep:
                    f = f + D * dt * _laplacian(f)
                    t += dt
            if mol in DECAY_PER_S:
                f = f * np.exp(-DECAY_PER_S[mol] * timestep)
            f = np.nan_to_num(f, nan=0.0, posinf=_FIELD_CAP, neginf=0.0)
            fields[mol] = np.clip(f, 0.0, _FIELD_CAP)

        # 3/4. sample local, reset exchange, remove dead cells
        cell_update = {}
        removed = []
        for cid, c in cells.items():
            if c.get('death'):
                removed.append(cid)   # exchange already deposited above
                continue
            ix, iy = self._bin(c.get('location', (0.0, 0.0)))
            local = {mol: float(fields[mol][ix, iy]) for mol in self.molecules}
            u = {'local': local}
            exch = c.get('exchange') or {}
            if exch:  # zero it out via merge-sum
                u['exchange'] = {mol: -counts for mol, counts in exch.items()}
            cell_update[cid] = u

        if removed:
            cell_update['_remove'] = removed
        # The fields store apply is additive, so return DELTAS (new - original),
        # matching viva-munk's DiffusionAdvection. Returning absolute arrays would
        # double the field every tick.
        delta_fields = {m: fields[m] - original[m] for m in self.molecules}
        return {'fields': delta_fields, 'cells': cell_update}


# --- workbench viewer contract (per-port meanings + units) ---
DiffusionField.contract = {
    'summary': "Soluble-field manager (tumor-tcell Fields + LocalField consolidated): deposits each "
               "cell's exchange counts into its field bin, diffuses + decays the 2D concentration fields, "
               "and samples each cell's local concentration back onto it.",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); reads each cell's location (um) and exchange "
                        "amounts (molecule counts).",
               'fields': "2D concentration grids per molecule (ng/mL) on an n_bins grid over bounds (um)."},
    'outputs': {'cells': "Per cell: sampled local concentrations (local, ng/mL) and exchange reset; cells "
                         "flagged dead are removed (_remove).",
                'fields': "Field deltas (ng/mL) after deposit + diffuse + decay."},
    'assumptions': ["IFNg: D = 1.25e-3 cm^2/day, 4.5 h e-fold decay. tumor_debris: D = 0.0864 cm^2/day, no "
                    "decay. Counts <-> concentration via Avogadro / molecular weight (IFNg 17000, debris "
                    "29000 g/mol) / bin volume (depth um)."],
}
