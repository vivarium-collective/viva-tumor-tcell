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

The showcase investigation (`investigations/tumor-tcell-showcase/`) recovers the
paper's central finding — IFNg-driven tumor phenotype transformation gates
therapeutic T-cell efficacy — at a tractable scale, across three studies (each
with an interactive population time-series and a spatial animation):

| Study | Result (single seed, ~10 h sim) |
|---|---|
| **tcell-efficacy-pd1** | T cells suppress tumor growth; **25% PD1+ suppresses more than 75% PD1+** (13.6% vs 10.2% vs no-T control) |
| **phenotype-conversion** | PDL1p tumor fraction **21.6% with active T cells vs 10.2% without**; IFNg field 6.7 vs 0 ng/mL |
| **killing-assay-cytotoxicity** | **25% cytotoxicity** vs a matched no-T control |

Effects are directional at reduced scale (tens of cells, tens of hours); the
studies document the parameters to reach the paper's full-scale figures
(1200 cells, 3 days).

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
python studies/tcell-efficacy-pd1/sims/run.py   # prints the verdict, writes viz/*.html
```
