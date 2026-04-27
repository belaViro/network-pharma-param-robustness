# Network Pharmacology Parameter Robustness Analysis

This repository contains a reproducible computational workflow for evaluating how robust network pharmacology targets and pathways are to analysis choices.

The project focuses on oxidative-stress-related network pharmacology of strawberry leaf compounds and is being developed toward a **dual-perturbation robustness framework**:

1. **Parameter perturbation**: Do targets and pathways remain stable when SwissTarget probability thresholds, STRING confidence thresholds, and Top-K hub cutoffs change?
2. **Compound perturbation**: Do targets and pathways remain stable when the input compound set is perturbed by leave-one-compound-out or compound bootstrap analysis?

The goal is not to prove a wet-lab pharmacological mechanism. The goal is to make computational network pharmacology results more transparent, reproducible, and less dependent on arbitrary thresholds or single-compound dominance.

## Research Question

Traditional network pharmacology often uses fixed thresholds such as:

| Parameter | Common choice | Risk |
|---|---:|---|
| SwissTarget probability | 0, 0.1 | arbitrary target inclusion |
| STRING confidence | 0.4, 0.7 | network topology changes |
| Hub target cutoff | Top-20, Top-30 | subjective target set size |
| Enrichment p value | 0.01, 0.05 | reporting convention |

This repository asks:

> Which targets and pathways remain stable when these choices change, and which findings are sensitive to parameter or compound selection?

## Method Overview

### Current implemented analysis

The current repository implements the parameter-robustness analysis:

- compound target prediction and disease-gene intersection
- STRING PPI network construction
- Top-K hub target selection
- GO/KEGG enrichment
- parameter-grid scanning
- robust target ranking by frequency and network topology
- pathway stability analysis
- PubMed post-hoc pathway classification validation

Current parameter grid:

| Parameter | Values |
|---|---|
| SwissTarget probability | `0`, `0.09`, `0.1` |
| STRING confidence | `0.30`, `0.35`, `0.40`, `0.45`, `0.50`, `0.70` |
| Top-K hubs | `25`, `30`, `35`, `40` |
| Enrichment p value | `0.01` |
| Total combinations | `72` |

### Planned dual-perturbation extension

The next version will extend the analysis from parameter robustness to dual robustness:

- **leave-one-compound-out**: remove one compound at a time and rerun robust target/pathway analysis
- **compound bootstrap**: resample compounds with replacement and evaluate target/pathway retention
- **compound influence score**: quantify how much each compound changes the robust target set
- **dual-robust targets**: targets stable under both parameter and compound perturbation

Key planned metrics:

```text
parameter_retention_rate = target retained across parameter-grid runs / total parameter runs

compound_retention_rate = target retained across compound-perturbation runs / total perturbation runs

compound_influence_score = 1 - Jaccard(Q1_full, Q1_without_compound_i)
```

## Data Scope

The project is being upgraded from an initial small compound set to a 19-compound LC-MS/MS strawberry leaf compound list. Trehalose and alpha,alpha-trehalose should be collapsed as the same compound in the final analysis.

Current generated outputs in this repository reflect the implemented parameter-grid analysis. The 19-compound compound-level uncertainty extension is the next analysis stage.

## Repository Structure

```text
code/
├── scripts/
│   ├── run_pipeline.py
│   └── analysis/
│       ├── run_sensitivity_analysis.py
│       ├── compare_robust_vs_parametric.py
│       ├── plot_sensitivity_figures.py
│       └── validate_with_pubmed.py
├── src/
│   ├── connectors/
│   ├── pipeline/
│   ├── sensitivity/
│   ├── enrichment/
│   ├── visualization/
│   └── outputs/
│       ├── shared/
│       ├── traditional/
│       └── sensitivity/
├── tests/
├── docs/
└── reproduce.sh
```

## Installation

Python 3.10 or later is recommended.

```bash
pip install -e .
```

Core dependencies are listed in `pyproject.toml` and `requirements.txt`.

## Reproduce Main Analysis

Run the parameter-grid sensitivity analysis:

```bash
python scripts/analysis/run_sensitivity_analysis.py
```

Compare robust targets with parameter-specific target sets:

```bash
python scripts/analysis/compare_robust_vs_parametric.py
```

Generate figures:

```bash
python scripts/analysis/plot_sensitivity_figures.py
```

Or run the bundled reproduction script:

```bash
bash reproduce.sh
```

## Main Outputs

Key output files are written to `src/outputs/sensitivity/`:

| File | Description |
|---|---|
| `param_grid_results.csv` | parameter-combination summary |
| `target_robustness.csv` | robust target frequency and topology metrics |
| `pathway_rank_stability.csv` | pathway retention and ranking stability |
| `robust_vs_param_comparison.json` | robust targets vs parameter-specific sets |
| `pubmed_validation.json` | PubMed post-hoc validation of pathway classification |
| `figures/` | publication-oriented sensitivity figures |

## Interpretation Boundary

This repository provides computational robustness evidence. It can support statements such as:

- a target is not specific to one arbitrary parameter setting
- a pathway remains stable across the scanned parameter space
- a future dual-robust target is not driven by one compound alone

It does **not** prove:

- direct compound-target binding
- pathway activation or inhibition in cells
- in vivo antioxidant efficacy
- a complete pharmacological mechanism

Those claims require molecular docking, molecular dynamics, cell assays, qPCR/Western blot, or animal experiments.

## Manuscript Direction

Recommended manuscript framing:

> A dual-perturbation robustness framework for network pharmacology target prioritization: parameter sensitivity and compound-level uncertainty in strawberry leaf oxidative-stress-related signaling.

Current target positioning:

- parameter-grid analysis: completed
- compound-level uncertainty: planned next stage
- external validation: PubMed post-hoc validation completed; OMIM/DisGeNET/CTD target validation planned
- likely manuscript level after dual-perturbation extension: stable third-quartile target, second-quartile attempt possible depending on validation depth

## Documentation

Detailed project notes are under `docs/`, especially:

- `docs/V2_稳三区投稿计划.md`
- `docs/parameter_grid_design.md`
- `docs/pipeline_overview.md`
- `docs/传统方法.md`

## License

See `LICENSE`.
