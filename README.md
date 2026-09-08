# viva-tumor-tcell

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph) port
of [tumor-tcell](https://github.com/vivarium-collective/tumor-tcell) — a
multiscale agent-based model of the tumor microenvironment (Hickey, Agmon et
al., *Cell Systems*) — as a **viva- workspace**. Collisions are provided by
[viva-munk](https://github.com/vivarium-collective/viva-munk); the tumor and
T-cell biology, the diffusing IFNg field, and the cell–cell neighbor exchange are
ported faithfully from the original vivarium-1.0 processes.

See [PORT_PLAN.md](PORT_PLAN.md) for the architecture, decisions, and milestones.

## Status

Milestone 1 (core seam): `tumor_tcell_agent` type; `TumorCellProcess`,
`TCellProcess`, `DiffusionField`, and `TumorTcellPhysics` (collisions via
viva-munk); the `tumor_tcell_basic` composite generator.

## Install (development)

```bash
uv venv .venv && source .venv/bin/activate
uv pip install -e .        # requires the sibling viva-munk checkout
```
