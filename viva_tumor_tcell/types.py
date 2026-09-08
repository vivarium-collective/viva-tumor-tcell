"""The ``tumor_tcell_agent`` schema.

Each cell in the ``cells`` map is a flat dict. We pin the map's value to an
explicit per-field schema (rather than letting inference guess) so that:

  * accumulator fields (timers, counts, stockpiles, internalized IFNg) are
    plain ``float`` — additive-delta apply;
  * set-semantics fields (presented/accepted membrane ligands, migration speed,
    the per-tick cytotoxic transfer) are ``set_float`` (viva-munk's replace-apply
    float) — the writer emits the target value and it replaces;
  * ``exchange`` is ``map[float]`` (additive; deposit then reset via a negative
    delta) and ``local`` is ``map[set_float]`` (field concentrations, replaced
    each tick). Molecule keys are pre-populated in the initial state so the
    generic map apply accepts them.

The embedded per-cell ``behavior`` process is NOT declared here — it rides along
in the state and resolves to a ProcessLink, exactly as viva-munk's ``grow_divide``
does. (Pinning the value to an opaque Node instead makes realize recurse through
that back-reference at build time — see PORT_PLAN.md.)

``set_float`` is registered by viva-munk (``register_pymunk_types``); it is
available on any core built from :func:`viva_tumor_tcell.core.build_core`.
"""

# field name -> bigraph-schema type string
AGENT_SCHEMA = {
    'id': 'string',
    'type': 'string',
    'cell_type': 'string',
    'cell_state': 'string',
    'death': 'string',                 # '' = alive; a reason string when dead
    'location': 'tuple[float,float]',
    'mass': 'float',
    'radius': 'float',
    'speed': 'set_float',              # migration speed magnitude (µm/s)
    # accumulators (additive delta)
    'internal_IFNg': 'float',
    'cell_state_count': 'float',
    'refractory_count': 'float',
    'total_cytotoxic_packets': 'float',
    'TCR_timer': 'float',
    'velocity_timer': 'float',
    'PDL1n_divide_count': 'float',
    'PD1n_divide_count': 'float',
    'PD1p_divide_count': 'float',
    'receive_cytotoxic': 'float',
    # membrane ligands, set each tick (replace)
    'present_TCR': 'set_float',
    'present_PD1': 'set_float',
    'present_PDL1': 'set_float',
    'present_MHCI': 'set_float',
    'accept_PDL1': 'set_float',
    'accept_MHCI': 'set_float',
    'accept_PD1': 'set_float',
    'accept_TCR': 'set_float',
    'transfer_cytotoxic': 'set_float',
    # soluble exchange (counts, additive) + sampled local concentrations (set)
    'exchange': 'map[float]',
    'local': 'map[set_float]',
}

# molecules whose exchange/local keys must exist from the start
MOLECULES = ('IFNg', 'tumor_debris')


def register_types(core):
    """Register ``tumor_tcell_agent`` as a normal object schema (NOT an opaque
    Node) so ``map[tumor_tcell_agent]`` resolves with a known field structure."""
    core.register_type('tumor_tcell_agent', dict(AGENT_SCHEMA))
    return core


def cells_store(cells):
    """Wrap a {cid: agent} dict as a schema-pinned map store."""
    return {'_type': 'map[tumor_tcell_agent]', **cells}
