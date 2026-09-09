"""Observables — scalar aggregates for the emitter / study Results preview.

A Step that reads the whole ``cells`` map and the ``fields`` grids each tick and
writes a flat map of scalar floats into an ``observables`` store. An emitter
(e.g. SQLiteEmitter) wired to that store persists one scalar-per-step series per
key, which the workbench's study Results tab previews as sparklines
(first/last/min/max). Values are absolute each tick (map[float] apply replaces,
not accumulates), so the series reflect true counts/fractions.
"""
from __future__ import annotations

from process_bigraph import Step

from ..run import population_counts, pdl1p_fraction, pd1p_fraction, ifng_max


class Observables(Step):
    description = (
        "Scalar aggregates for the emitter + study Results/Readouts: per-tick cell "
        "counts by type/state, PDL1p / PD1p fractions, and peak IFNg, written to an "
        "observables store for a SQLiteEmitter to persist."
    )
    config_schema = {}

    def inputs(self):
        return {'cells': 'map[tumor_tcell_agent]',
                'fields': {'_type': 'map', '_value': {'_type': 'array', '_data': 'float'}}}

    def outputs(self):
        # set_float (replace semantics), NOT plain float — a map[float] store
        # apply is additive (fluxes accumulate), which would sum counts across
        # ticks. set_float overwrites, so each tick holds the true value.
        return {'observables': {'_type': 'map', '_value': 'set_float'}}

    def update(self, state):
        cells = state.get('cells') or {}
        fields = state.get('fields') or {}
        frame = {'cells': cells, 'fields': fields}
        counts = population_counts(frame)   # tumor_/tcell_/dendritic_ counts by state
        obs = {k: float(v) for k, v in counts.items()}
        obs['pdl1p_fraction'] = float(pdl1p_fraction(frame))
        obs['pd1p_fraction'] = float(pd1p_fraction(frame))
        obs['ifng_max'] = float(ifng_max(frame) or 0.0)
        return {'observables': obs}


# --- workbench viewer contract (per-port meanings + units) ---
Observables.contract = {
    'summary': "Scalar aggregates for the emitter + study Results preview: per-tick cell counts by "
               "type/state, PDL1p / PD1p fractions, and peak IFNg — written to an observables store "
               "for a SQLiteEmitter to persist.",
    'inputs': {'cells': "All cells (map[tumor_tcell_agent]); counted by cell_type and cell_state.",
               'fields': "2D concentration grids (ng/mL) — peak IFNg is read from the IFNg grid."},
    'outputs': {'observables': "Flat map of scalar floats (absolute each tick): tumor_total/PDL1n/PDL1p, "
                               "tcell_total/PD1n/PD1p, dendritic_total/inactive/active, pdl1p_fraction, "
                               "pd1p_fraction, ifng_max."},
    'assumptions': ["Absolute values each tick (map[float] apply replaces); fractions are 0 when the "
                    "denominator population is empty."],
}
