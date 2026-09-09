"""In-composite Visualization Steps for the tumor-tcell showcase.

These are ``process_bigraph.visualization.Visualization`` subclasses (via
viva-superpowers) that live *inside* the composite document. During a run the
workbench's Composite Explorer calls ``render_results(composite)`` at the end,
which finds every embedded Visualization instance and calls its ``render()`` —
the resulting HTML is what the Explorer's **Visualizations** tab shows.

They are stateful (legacy streaming contract): each tick they accumulate the
``cells`` map (and ``fields``) into their own buffers, and ``update()`` returns
the freshly rendered HTML. They reuse the same Plotly figure builders as the
study figures (:mod:`viva_tumor_tcell.viz`), so the in-Explorer figures match
the paper's colors and forms.

Three views, wired into every composite via :func:`viz_steps`:

* ``PopulationTimeseries`` — total + per-phenotype cell counts over time.
* ``PhenotypeFractions``   — the paper's central readout: PDL1+ tumor fraction
  and PD1+ T-cell fraction over time (the IFNg-driven conversion).
* ``SpatialLayout``        — animated scatter of cells over the IFNg field, the
  interactive analogue of the paper's videos.
"""
from __future__ import annotations

from viva_superpowers.visualization import Visualization

from .. import viz as _viz
from ..run import population_counts

TIMESTEP = 60.0


def _fig_html(fig, height='420px'):
    """Serialize a Plotly figure to a standalone HTML document for the iframe."""
    return fig.to_html(include_plotlyjs='cdn', full_html=True,
                       config={'displayModeBar': False, 'responsive': True},
                       default_width='100%', default_height=height)


def _empty(msg):
    return f'<p style="color:#888;font:14px system-ui;padding:12px">{msg}</p>'


def _time_h(state, step):
    """Prefer the wired global_time; fall back to tick x TIMESTEP (hours)."""
    t = state.get('time')
    if t is None:
        t = step * TIMESTEP
    return float(t) / 3600.0


class PopulationTimeseries(Visualization):
    """Total + per-phenotype cell counts vs time (mirrors population_group_plot)."""

    description = (
        "Total and per-phenotype cell counts over time (tumor PDL1n/PDL1p, T cell "
        "PD1n/PD1p, dendritic) as an interactive Plotly line chart."
    )
    config_schema = {
        'title': {'_type': 'string', '_default': 'Cell populations over time'},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._times = []
        self._pops = []
        self._step = 0

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]', 'time': 'float'}

    def update(self, state):
        cells = state.get('cells') or {}
        self._times.append(_time_h(state, self._step))
        self._step += 1
        self._pops.append(population_counts({'cells': cells}))
        return {'html': self._render()}

    def _render(self):
        if not self._pops:
            return _empty('No population data yet.')
        last = self._pops[-1]
        keys = ['tumor_total', 'tumor_PDL1n', 'tumor_PDL1p',
                'tcell_total', 'tcell_PD1n', 'tcell_PD1p']
        if last.get('dendritic_total'):
            keys += ['dendritic_total', 'dendritic_inactive', 'dendritic_active']
        # drop count-series that are flat-zero across the run (keeps legend clean)
        keys = [k for k in keys if any(p.get(k, 0) for p in self._pops)]
        title = (self.config or {}).get('title', 'Cell populations over time')
        fig = _viz.population_group_figure(self._times, self._pops, keys, title)
        return _fig_html(fig)

    @classmethod
    def demo(cls):
        return {'cells': {'u': {'cell_type': 'tumor', 'cell_state': 'PDL1n'}},
                'time': 0.0}

    @classmethod
    def is_visualization(cls):
        return True


class PhenotypeFractions(Visualization):
    """PDL1+ tumor fraction and PD1+ T-cell fraction over time.

    This is the paper's headline readout: efficacy tracks the *rate* of the
    IFNg-driven PDL1n -> PDL1p tumor conversion (and T-cell PD1n -> PD1p
    exhaustion), not raw kill count."""

    description = (
        "PDL1+ tumor fraction and PD1+ T-cell fraction over time (Plotly) — the "
        "paper's headline readout, tracking the rate of IFNg-driven phenotype "
        "conversion rather than raw kill count."
    )
    config_schema = {
        'title': {'_type': 'string', '_default': 'Phenotype conversion over time'},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._times = []
        self._pdl1p = []
        self._pd1p = []
        self._step = 0

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]', 'time': 'float'}

    def update(self, state):
        cells = (state.get('cells') or {}).values()
        self._times.append(_time_h(state, self._step))
        self._step += 1
        tum = [c for c in cells if c.get('cell_type') == 'tumor']
        tc = [c for c in cells if c.get('cell_type') == 't-cell']
        self._pdl1p.append(
            sum(c.get('cell_state') == 'PDL1p' for c in tum) / len(tum) if tum else 0.0)
        self._pd1p.append(
            sum(c.get('cell_state') == 'PD1p' for c in tc) / len(tc) if tc else 0.0)
        return {'html': self._render()}

    def _render(self):
        if not self._times:
            return _empty('No phenotype data yet.')
        series = {'PDL1+ tumor fraction': self._pdl1p}
        if any(self._pd1p):
            series['PD1+ T-cell fraction'] = self._pd1p
        title = (self.config or {}).get('title', 'Phenotype conversion over time')
        fig = _viz.timeseries_figure(self._times, series, title, yaxis='fraction')
        return _fig_html(fig)

    @classmethod
    def demo(cls):
        return {'cells': {'u': {'cell_type': 'tumor', 'cell_state': 'PDL1n'}},
                'time': 0.0}

    @classmethod
    def is_visualization(cls):
        return True


class SpatialLayout(Visualization):
    """Animated GIF of cells over the IFNg field (the paper's video analogue).

    New-style: ``accumulate`` buffers a compact per-tick frame; ``render`` builds
    the animated GIF once at end-of-run via the matplotlib method — true circles
    in µm data-coordinates (correct cell sizes) with every frame drawn, so it
    matches the study spatial figures rather than a pixel-sized scatter."""

    description = (
        "Animated GIF of cells (true circles in µm) over the IFNg field, colored "
        "by type/state — the interactive analogue of the paper's videos."
    )
    config_schema = {
        'title': {'_type': 'string', '_default': 'Spatial layout'},
        'bounds_x': {'_type': 'float', '_default': 400.0},
        'bounds_y': {'_type': 'float', '_default': 400.0},
        'field': {'_type': 'string', '_default': 'IFNg'},
        'max_frames': {'_type': 'integer', '_default': 120},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._frames = []
        self._step = 0

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]',
                'fields': {'_type': 'map', '_value': {'_type': 'array',
                           '_data': 'float'}},
                'time': 'float'}

    def accumulate(self, state):
        cfg = self.config or {}
        field = cfg.get('field', 'IFNg')
        cells = state.get('cells') or {}
        fields = state.get('fields') or {}
        # keep a compact snapshot: only what the figure needs
        frame_cells = {}
        for cid, c in cells.items():
            if not isinstance(c, dict):
                continue
            frame_cells[cid] = {
                'cell_type': c.get('cell_type'),
                'cell_state': c.get('cell_state'),
                'location': list(c.get('location', (0.0, 0.0))),
                'radius': float(c.get('radius', 5.0) or 5.0),
            }
        arr = fields.get(field)
        frame_field = {}
        if arr is not None:
            import numpy as np
            frame_field[field] = np.asarray(arr, dtype=float).tolist()
        self._frames.append({'time': self._step * TIMESTEP,
                             'cells': frame_cells, 'fields': frame_field})
        self._step += 1

    def render(self):
        if not self._frames:
            return _empty('No spatial data yet.')
        cfg = self.config or {}
        bounds = (float(cfg.get('bounds_x', 400.0)), float(cfg.get('bounds_y', 400.0)))
        field = cfg.get('field', 'IFNg')
        title = cfg.get('title', 'Spatial layout')
        try:
            # matplotlib-GIF method — returns a full HTML page (autoplaying gif).
            return _viz.spatial_gif_html(
                self._frames, bounds, title, field=field,
                max_frames=int(cfg.get('max_frames', 120)))
        except Exception as e:  # noqa: BLE001
            return _empty(f'spatial render failed: {e}')

    @classmethod
    def demo(cls):
        return {'cells': {'u': {'cell_type': 'tumor', 'cell_state': 'PDL1n',
                                'location': [200.0, 200.0], 'radius': 7.5}},
                'fields': {'IFNg': [[0.0, 0.0], [0.0, 0.0]]},
                'time': 0.0}

    @classmethod
    def is_visualization(cls):
        return True


for _cls in (PopulationTimeseries, PhenotypeFractions, SpatialLayout):
    _cls.__pb_kind__ = 'visualization'
    _cls.__pb_aliases__ = [_cls.__name__]


# --- workbench loom contracts (per-port meanings) ---
PopulationTimeseries.contract = {
    'summary': "Total + per-phenotype cell counts over time, rendered as an interactive Plotly line chart.",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); counted by cell_type and cell_state each tick.",
               'time': "Simulation time (s); the x-axis (converted to hours)."},
    'outputs': {'html': "Rendered Plotly HTML (line per count series: tumor total/PDL1n/PDL1p, T cell "
                        "total/PD1n/PD1p, and dendritic when present)."},
    'assumptions': ["Accumulates one count snapshot per tick; flat-zero series are dropped from the legend."],
}
PhenotypeFractions.contract = {
    'summary': "PDL1+ tumor fraction and PD1+ T-cell fraction over time — the paper's IFNg-driven "
               "conversion readout.",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); the PDL1p tumor fraction and PD1p T-cell "
                        "fraction are computed each tick.",
               'time': "Simulation time (s); the x-axis (converted to hours)."},
    'outputs': {'html': "Rendered Plotly HTML (fraction in 0-1 vs time, one line per fraction)."},
    'assumptions': ["Fractions are 0 when the denominator population (tumors / T cells) is empty."],
}
SpatialLayout.contract = {
    'summary': "Animated GIF of cells over the IFNg field — true circles in um data-coordinates, colored "
               "by type/state; the interactive analogue of the paper's videos.",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); reads location (um), radius (um), cell_type, "
                        "cell_state per tick.",
               'fields': "2D concentration grids per molecule (ng/mL); the IFNg grid is drawn as the "
                         "background heatmap.",
               'time': "Simulation time (s); frame timestamps."},
    'outputs': {'html': "A standalone HTML page embedding an autoplaying, looping GIF (matplotlib)."},
    'assumptions': ["Buffers a compact frame per tick; the GIF is rendered once at end-of-run "
                    "(max_frames subsampling for a smooth, lightweight animation)."],
}


def register_visualizations(core):
    """Register the three viz Steps on a core."""
    core.register_link('PopulationTimeseries', PopulationTimeseries)
    core.register_link('PhenotypeFractions', PhenotypeFractions)
    core.register_link('SpatialLayout', SpatialLayout)
    return core


def viz_steps(bounds=(400.0, 400.0), field='IFNg'):
    """Return the three embedded Visualization step nodes for a composite doc.

    Wired to read the shared ``cells`` map (+ ``fields`` for the spatial view)
    and ``global_time``; each writes its rendered HTML to its own viz store."""
    bx, by = float(bounds[0]), float(bounds[1])
    return {
        'viz_populations': {
            '_type': 'step', 'address': 'local:PopulationTimeseries',
            'config': {'title': 'Cell populations over time'},
            'inputs': {'cells': ['cells'], 'time': ['global_time']},
            'outputs': {'html': ['viz', 'populations_html']},
        },
        'viz_phenotypes': {
            '_type': 'step', 'address': 'local:PhenotypeFractions',
            'config': {'title': 'Phenotype conversion over time'},
            'inputs': {'cells': ['cells'], 'time': ['global_time']},
            'outputs': {'html': ['viz', 'phenotypes_html']},
        },
        'viz_spatial': {
            '_type': 'step', 'address': 'local:SpatialLayout',
            'config': {'title': 'Spatial layout over the IFNg field',
                       'bounds_x': bx, 'bounds_y': by, 'field': field},
            'inputs': {'cells': ['cells'], 'fields': ['fields'],
                       'time': ['global_time']},
            'outputs': {'html': ['viz', 'spatial_html']},
        },
    }
