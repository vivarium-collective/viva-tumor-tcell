# viva-tumor-tcell — port plan

Porting [tumor-tcell](https://github.com/vivarium-collective/tumor-tcell) (a
vivarium-1.0 / vivarium-core ABM of the tumor microenvironment, from Hickey,
Agmon et al., *Cell Systems*) to a process-bigraph "viva-" workspace, using
**viva-munk** for the collision physics and porting everything else faithfully.

## Decisions (locked with the user)

- **Scope: tractable-faithful.** Port every experiment *type* (the 3 PD1+
  conditions, the 4-condition in-vitro killing/cytotoxicity assay, the lymph
  node) plus all plots, and run the committed studies at reduced scale
  (fewer cells / shorter time) that reproduces the qualitative results. Each
  study documents the exact params to scale to the paper figure; full-scale
  runs are a follow-up (mini/overnight).
- **Fields: port tumor-tcell's own `Fields`/`LocalField` faithfully** — its
  exact diffusion rates (IFNg 1.25e-3 cm²/day, debris 0.0864 cm²/day) and
  decay (IFNg 4.5 h e-fold). viva-munk is used ONLY for multibody collisions.

## Architecture

tumor-tcell's `Neighbors` process did **two** jobs: pymunk circle collisions,
**and** the biological neighbor-exchange (`present`/`accept`/`transfer`/
`receive`, "T-cell picks nearest tumor, tumor collects all T-cells"). viva-munk's
`PymunkProcess` covers only the collisions. So:

- **`TumorTcellPhysics`** (this repo) wraps a `viva_munk` `PymunkProcess`
  instance as its collision engine — builds flat pymunk bodies from the cells,
  steps viva-munk's real pymunk space (walls, jitter, substeps, collisions),
  and **also** performs the ported neighbor-exchange bookkeeping.
- Migration model is preserved: each tick a cell's velocity is set to
  `speed × random_direction` (tumor-tcell's persistent-random-walk migration),
  then integrated by the physics engine.

### Store layout (deliberate translation)

tumor-tcell nested each cell into `internal` / `boundary` / `neighbors` /
`globals` sub-stores. Here each cell is a **flat `tumor_tcell_agent`** with those
ports flattened into namespaced keys (`cell_state`, `internal_IFNg`, `speed`,
`present_TCR`, `accept_MHCI`, `receive_cytotoxic`, `exchange`, `local`, …). The
biology (equations, parameters, transitions, rates) is preserved exactly; only
the store shape changes, matching viva-munk's flat-agent + per-cell-process idiom.

Per-cell behavior processes (`TumorCellProcess`, `TCellProcess`) are embedded in
each agent (key `behavior`) and read the whole `cells` map via `agent_id`, like
viva-munk's `GrowDivide`. Division is done by the process writing `_add`/`_remove`
into the `cells` map (replaces vivarium's `MetaDivision` + `Remove` + `_divider`).

### Field consolidation

tumor-tcell split field handling across `Fields` (diffuse/decay + sample local
into cells) and `LocalField` (deposit each cell's `exchange` counts into the
field bin, reset exchange). Ported as **one** `DiffusionField` process per tick:
(1) deposit exchange → field bins (count→concentration via Avogadro/MW/bin-vol),
reset exchange; (2) diffuse + decay; (3) sample each cell's local concentration.
Same math, fewer moving parts.

## Faithfulness notes / known quirks

- Original `t_cell.py` decrements the stockpile on transfer by writing key
  `internal.cytotoxic_packets` (line ~513) which does not match its schema key
  `total_cytotoxic_packets` — so in the original the stockpile is never actually
  drawn down on transfer. We port the **intended** semantics (decrement
  `total_cytotoxic_packets`) and note the difference.

## Milestones

- **M1 — core seam (this commit):** `tumor_tcell_agent` type; `TumorCellProcess`,
  `TCellProcess`, `DiffusionField`, `TumorTcellPhysics`; `build_core()`; a
  `tumor_tcell_basic` composite generator; a test that builds + runs it and shows
  IFNg-driven tumor state switch + T-cell-mediated tumor death. Physics via viva-munk.
- **M2 — full environment + remaining processes:** dendritic cell, lymph-node
  transport; the full `tumor_microenvironment` generator + CODEX-style seeding.
- **M3 — investigation + studies + viz + publish:** killing/cytotoxicity assay
  (4 conditions), the 3 PD1+ conditions, lymph node; population/death/division +
  spatial-snapshot animations; runs into `.pbg/runs.jsonl`; read-only workbench.
