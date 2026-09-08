# viva-tumor-tcell

<!-- BEGIN:dashboard -->
> **📊 Read-only dashboard** (after first push + Pages enable):
> **https://vivarium-collective.github.io/viva-tumor-tcell/dashboard/**
> — browse the showcase investigation, its studies, and their interactive
> figures with no server. Rebuilt on every push to `main` by
> `.github/workflows/publish-dashboard.yml`.
<!-- END:dashboard -->

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph) port
of [tumor-tcell](https://github.com/vivarium-collective/tumor-tcell) — a
multiscale agent-based model of the tumor microenvironment (Hickey, Agmon et
al., *Cell Systems*) — as a **viva- workspace**. Collisions are provided by
[viva-munk](https://github.com/vivarium-collective/viva-munk); the tumor,
T-cell, and dendritic-cell biology, the diffusing IFNg / tumor-debris fields, and
the cell–cell neighbor exchange are ported faithfully from the original
vivarium-1.0 processes.

See [PORT_PLAN.md](PORT_PLAN.md) for the architecture, decisions, and milestones.

## What it reproduces

The showcase investigation (`investigations/tumor-tcell-showcase/`) reconstructs
the paper's mechanism and reports every study across **6 replicate seeds** (mean
± std), separating what holds robustly at tractable scale from what is
scale-limited. Each study has an interactive mean±std time-series and a spatial
animation.

**Robust across seeds (the cell + cytokine mechanism):**

| Study | Result (6 seeds, ~10 h sim) |
|---|---|
| **tcell-exhaustion** | Exhausted (PD1+) fraction tracks the starting PD1+ fraction: **63 ± 22% (75% start) vs 16 ± 11% (25% start)** — 6/6 seeds |
| **phenotype-conversion** | Active T cells secrete a measurable **IFNg field in 6/6 seeds; 0/6 without** |
| **killing-assay-cytotoxicity** | **cytotoxicity ~11 ± 8%** vs matched no-T control (positive in 5/6 seeds) |

**Scale-limited (documented, not claimed):** population tumor-count suppression
and the 25%<75% PD1+ efficacy *ordering* do **not** separate from growth noise at
this scale (25%>75% in only 3/6 seeds), and the PDL1p conversion *fraction* is
directional only (2/6). These need paper-scale runs (1200 cells, 3 days) — a
mini/overnight job — to reach significance. See `PORT_PLAN.md` and the
investigation's caveats.

## Processes & composites

`TumorCellProcess`, `TCellProcess`, `DendriticCellProcess` (per-cell biology);
`DiffusionField` (IFNg / tumor_debris); `TumorTcellPhysics` (viva-munk collisions
+ neighbor exchange). Composite generators: `tumor_tcell_basic`,
`tumor_microenvironment`, `killing_assay`, `lymph_node`.

## Install (development)

```bash
uv venv .venv && source .venv/bin/activate
uv pip install -e .        # requires the sibling viva-munk checkout
pytest                     # 17 tests
```

## Run a study locally

```bash
python studies/tcell-exhaustion/sims/run.py   # prints the verdict, writes viz/*.html
```
