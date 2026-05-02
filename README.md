# Network Pharmacology Parameter Robustness Analysis

**English** | [中文](./README_CN.md)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> **Dual-perturbation framework for robust target identification in network pharmacology:** parameter-grid sensitivity analysis combined with leave-one-compound-out validation, applied to strawberry leaf compounds against oxidative stress.

## Overview

This repository contains the complete reproducibility package for a systematic robustness analysis of network pharmacology pipelines. We address a critical but often overlooked problem: **parameter sensitivity**—the same dataset can yield vastly different targets and pathways depending on arbitrary threshold choices.

### Key Contributions

- **Dual-perturbation framework**: Simultaneously perturbs (1) analytical parameters and (2) input compounds to identify targets robust against both sources of variation.
- **112 parameter combinations**: Exhaustive grid over probability thresholds (7) × STRING confidence (4) × hub-gene cutoffs (4).
- **2,016 compound-sensitivity runs**: Leave-one-out analysis across 18 strawberry leaf compounds.
- **Three-tier target classification**: Core dual-robust (29 targets) vs. compound-specific vs. parameter-conditional.

### Key Findings

| Metric | Value |
|--------|-------|
| Total parameter combinations | **112** (7 × 4 × 4) |
| Total compound LOO combinations | **2,016** (18 × 112) |
| SwissTarget union targets | 684 |
| Disease-intersected targets | 639 |
| Q1 high-confidence candidates | 50 |
| **Dual-robust targets** | **29** (58% retention) |
| **Tier 1 target** | **STAT3** (freq = 100%, LOO = 100%) |
| Dual-robust pathways | 58 |
| Top compound influence | Cardiac glycoside (score = 0.306) |

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/belaViro/network-pharma-param-robustness.git
cd network-pharma-param-robustness
```

### 2. Install dependencies

```bash
# Option A: pip
pip install -r requirements.txt

# Option B: editable install (recommended)
pip install -e .
```

Requires Python ≥ 3.10.

### 3. Run the main analyses

```bash
# Parameter grid sensitivity (112 combinations)
python scripts/analysis/run_sensitivity_analysis.py

# Compound leave-one-out sensitivity (2,016 combinations)
python scripts/analysis/run_compound_sensitivity_analysis.py
```

### 4. One-click reproduction

```bash
bash reproduce.sh
```

This executes the full pipeline: parameter grid → robust-vs-parametric comparison → pathway quadrant augmentation → figure regeneration.

---

## Repository Structure

```
.
├── scripts/
│   ├── analysis/
│   │   ├── run_sensitivity_analysis.py           # Main: 112-combo parameter grid
│   │   ├── run_compound_sensitivity_analysis.py  # Main: 18-compound LOO
│   │   ├── compare_robust_vs_parametric.py       # Robust vs. parametric comparison
│   │   ├── augment_pathway_robustness.py         # Pathway quadrant enhancement
│   │   ├── plot_sensitivity_figures.py           # Publication figures
│   │   ├── validate_with_pubmed.py               # External validation (PubMed)
│   │   └── ...                                   # Additional diagnostic scripts
│   ├── run_full_pipeline.py                      # End-to-end pipeline
│   └── run_traditional_pipeline.py               # Baseline comparison
│
├── src/
│   ├── config.py                                 # Central configuration
│   ├── connectors/                               # API wrappers
│   │   ├── swiss_target.py                       # SwissTargetPrediction
│   │   ├── string_connector.py                   # STRING PPI
│   │   ├── pubchem.py                            # PubChem
│   │   └── pubmed_connector.py                   # PubMed validation
│   ├── enrichment/
│   │   ├── enrichment_analysis.py                # GO/KEGG enrichment (gseapy)
│   │   └── plotting.py                           # Enrichment visualization
│   ├── visualization/                            # Plotting utilities
│   ├── data/                                     # Reference data
│   │   ├── go-basic.obo                          # Gene Ontology ontology
│   │   └── kegg_pathway_classes.json             # KEGG pathway taxonomy
│   └── outputs/                                  # All analysis outputs
│       ├── sensitivity/                          # Parameter sensitivity results
│       │   ├── target_robustness.csv
│       │   ├── pathway_robustness.csv
│       │   ├── param_grid_results.csv
│       │   └── compound_sensitivity/             # LOO compound results
│       ├── shared/                               # Cross-stage shared data
│       └── traditional/                          # Baseline method outputs
│
├── tests/                                        # Unit & integration tests
├── docs/                                         # Project documentation
├── requirements.txt
├── pyproject.toml
└── reproduce.sh                                  # One-click reproduction
```

---

## Output Files

### Core Results (Parameter Sensitivity)

| File | Description |
|------|-------------|
| `src/outputs/sensitivity/target_robustness.csv` | Target frequency & degree statistics across 112 combos |
| `src/outputs/sensitivity/pathway_robustness.csv` | Pathway robustness metrics (signaling backbone, direct OX, cross-regulation) |
| `src/outputs/sensitivity/param_grid_results.csv` | Per-combo enrichment summary |
| `src/outputs/sensitivity/pathway_rank_stability.csv` | Pathway rank consistency across parameter space |

### Core Results (Compound Sensitivity)

| File | Description |
|------|-------------|
| `src/outputs/sensitivity/compound_sensitivity/compound_influence_scores.csv` | Influence score per compound (0–1) |
| `src/outputs/sensitivity/compound_sensitivity/dual_robust_targets.csv` | 29 dual-robust targets with LOO retention rates |
| `src/outputs/sensitivity/compound_sensitivity/loo_pathway_robustness.csv` | Pathway retention under compound dropout |

### Checkpoints

- `src/outputs/sensitivity/compound_sensitivity/checkpoints/combo_results/` — 2,016 individual JSON files (one per compound × combo), enabling full audit and resume.

---

## Parameter Grid

| Parameter | Values | Description |
|-----------|--------|-------------|
| `Prob` (SwissTarget) | 0.0, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20 | Probability threshold for target inclusion |
| `STRING` (PPI) | 0.4, 0.5, 0.7, 0.9 | STRING combined score threshold for edge retention |
| `TopK` (Hub) | 10, 20, 30, 50 | Number of top degree-ranked genes selected as hubs |
| `Pvalue` (fixed) | 0.01 | Enrichment significance cutoff |

**Total combinations: 7 × 4 × 4 = 112**

---

## Strawberry Leaf Compounds (n = 18)

LC-MS/MS identified 19 raw records; after deduplication (Trehalose ≡ alpha,alpha-Trehalose), 18 canonical compounds:

1. Miquelianin
2. Sucrose
3. Trehalose
4. Oleamide
5. Catechin
6. Indole-3-acrylic acid
7. DL-Tryptophan
8. (3β,9ξ)-3-(β-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide *(cardiac glycoside)*
9. Denin
10. 2-Methylglutaric acid
11. Procyanidin B1
12. Scutellarin
13. Choline
14. DL-Malic acid
15. alpha-Linolenic acid
16. Palmitic acid
17. Linoleic acid
18. Stearic acid

---

## Target Classification

| Tier | Criteria | Count | Examples |
|------|----------|-------|----------|
| **Tier 1** | param_freq = 100% AND LOO_freq = 100% | 1 | **STAT3** |
| **Tier 2** | Dual-robust (param-robust & compound-robust) | 28 | AKT1, MAPK3, EGFR, ... |
| **Compound-specific** | LOO_freq drops to 0 when cardiac glycoside removed | — | CASP3, MAPK1, PIK3CA |
| **Parameter-conditional** | Low-frequency dual-robust | — | MET, IL2, KIT |

---

## Dependencies

Core analysis:
- pandas ≥ 2.0
- numpy ≥ 1.24
- networkx ≥ 3.0
- scipy ≥ 1.10
- gseapy ≥ 1.0 (enrichment)
- goatools ≥ 1.3 (GO analysis)

Visualization:
- matplotlib ≥ 3.6
- seaborn ≥ 0.12

Data acquisition (optional):
- playwright ≥ 1.40
- selenium ≥ 4.0

See `requirements.txt` or `pyproject.toml` for full details.

---

## Data Availability

All source code, raw predictions, intermediate checkpoints, and final results are included in this repository.

- **Repository**: https://github.com/belaViro/network-pharma-param-robustness
- **Complete reproducibility package** (~33 MB compressed): Available via GitHub Release [v1.0.0](https://github.com/belaViro/network-pharma-param-robustness/releases)

External databases used:
- SwissTargetPrediction (http://www.swisstargetprediction.ch/)
- STRING v12 (https://string-db.org/)
- Enrichr (https://maayanlab.cloud/Enrichr/)
- Gene Ontology (http://purl.obolibrary.org/obo/go-basic.obo)

---

## Citation

If you use this code or framework, please cite:

> *Parameter robustness analysis for network pharmacology target identification: a dual-perturbation framework applied to strawberry leaf compounds against oxidative stress.*

---

## License

[MIT License](./LICENSE)
