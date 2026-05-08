# Network Pharmacology Dual-Perturbation Robustness Framework

**English** | [中文](./README_CN.md)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> **Dual-perturbation framework for robust target prioritization in network pharmacology:** parameter-grid perturbation improves fixed-threshold target selection, while leave-one-compound-out perturbation probes multi-compound mechanism dependence in strawberry leaf antioxidant analysis.

## Overview

This repository contains the reproducibility package for a dual-perturbation robustness analysis of network pharmacology pipelines. The project addresses two linked sources of uncertainty: **parameter sensitivity**, where arbitrary thresholds can change targets and GO BP/KEGG functional enrichment results, and **compound-input dependence**, where a candidate mechanism may be driven by one or a few input compounds rather than the whole multi-compound system.

### Key Contributions

- **Parameter perturbation as method improvement**: replaces a single fixed-threshold TopK list with parameter-space recurrence and network-position statistics.
- **Compound perturbation as mechanism analysis**: uses leave-one-compound-out runs to identify targets and GO BP/KEGG functional enrichment terms that do not depend on one specific input compound.
- **112 parameter combinations**: exhaustive grid over SwissTargetPrediction probability thresholds (7) × STRING confidence thresholds (4) × hub-gene cutoffs (4).
- **2,016 compound-perturbation runs**: 18 structure-normalized strawberry leaf compounds × 112 parameter combinations.
- **Auditable evidence chain**: PR-QI target/term prioritization, Top30 ablation, PubMed co-occurrence audit, degree-matched permutation, and docking seed-stability outputs.

### Key Findings

| Metric | Value |
|--------|-------|
| Total parameter combinations | **112** (7 × 4 × 4) |
| Total compound LOO combinations | **2,016** (18 × 112) |
| SwissTarget union targets | 684 |
| Positive-probability disease targets | 286 |
| Main reporting threshold targets (Probability > 0.05) | 219 |
| PR-QI parameter-robust targets | 50 |
| PR-QI GO BP/KEGG functional terms | 70 |
| **Dual-robust targets** | **29** (58% retention) |
| **Tier 1 target** | **STAT3** (freq = 100%, LOO = 100%) |
| Dual-robust GO BP/KEGG functional terms | 58 |
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

# Optional: fresh PubChem/SwissTargetPrediction acquisition
pip install -e ".[acquisition]"
```

Requires Python ≥ 3.10.

### 3. Run the main analyses manually

```bash
# Traditional fixed-threshold baseline
python scripts/run_traditional_pipeline.py

# Parameter grid sensitivity (112 combinations)
python scripts/analysis/run_sensitivity_analysis.py

# Compound leave-one-out perturbation (2,016 combinations)
python scripts/analysis/run_compound_sensitivity_analysis.py

# External evidence audit and Top30 ablation bundle
python scripts/analysis/robustness_validation_bundle.py

# Manuscript figure regeneration
python scripts/analysis/generate_paper_figures.py

# Supplementary table export
python scripts/analysis/export_supplementary_tables.py
```

### 4. One-click reproduction

```bash
# Full manuscript workflow from packaged inputs (default)
python reproduce.py

# Same workflow through the Bash wrapper
bash reproduce.sh

# Windows PowerShell wrapper
.\reproduce.ps1

# Fast parameter-robustness workflow only
python reproduce.py --core
```

By default, the reproduction entry point runs the manuscript-facing workflow from packaged repository inputs: traditional fixed-threshold baseline → 112-combo parameter grid → robust-vs-parametric comparison → GO BP/KEGG functional-term quadrant augmentation → compound LOO perturbation → external evidence and degree-matched validation bundle → final manuscript figures → supplementary tables.

Useful modes:

| Command | Use |
|---------|-----|
| `python reproduce.py` | Full manuscript workflow from included inputs and caches |
| `python reproduce.py --core` | Fast parameter-robustness workflow only |
| `python reproduce.py --skip-install` | Reuse the current environment without reinstalling the package |
| `python reproduce.py --pubmed-mode query` | Refresh PubMed co-occurrence counts through NCBI E-utilities |
| `python reproduce.py --full --include-acquisition` | Also rerun PubChem/SwissTargetPrediction acquisition |

PubMed validation uses cached counts by default and does not query the network. Fresh PubMed or acquisition runs require network access and may take substantially longer.

---

## Repository Structure

```
.
├── scripts/
│   ├── analysis/
│   │   ├── run_sensitivity_analysis.py           # Main: 112-combo parameter grid
│   │   ├── run_compound_sensitivity_analysis.py  # Main: 18-compound LOO
│   │   ├── compare_robust_vs_parametric.py       # Robust vs. parametric comparison
│   │   ├── augment_pathway_robustness.py         # GO BP/KEGG robustness quadrant enhancement
│   │   ├── robustness_validation_bundle.py        # Top30 ablation and external evidence audit
│   │   ├── generate_paper_figures.py             # Manuscript figure generation
│   │   ├── export_supplementary_tables.py        # Supplementary workbook export
│   │   ├── validate_with_pubmed.py               # Standalone PubMed diagnostic
│   │   └── ...                                   # Additional diagnostic scripts
│   ├── run_pipeline.py                           # Optional PubChem/SwissTarget acquisition
│   ├── run_full_pipeline.py                      # Legacy staged acquisition + baseline pipeline
│   └── run_traditional_pipeline.py               # Baseline comparison
│
├── src/
│   ├── config.py                                 # Central configuration
│   ├── connectors/                               # API wrappers
│   │   ├── swiss_target.py                       # SwissTargetPrediction
│   │   ├── string_connector.py                   # STRING PPI
│   │   ├── pubchem.py                            # PubChem
│   │   └── pubmed_connector.py                   # PubMed co-occurrence audit
│   ├── enrichment/
│   │   ├── enrichment_analysis.py                # GO/KEGG enrichment (gseapy)
│   │   └── plotting.py                           # Enrichment visualization
│   ├── visualization/                            # Plotting utilities
│   ├── data/                                     # Reference data
│   │   ├── oxidative_stress_go_terms.txt         # Curated GO oxidative-stress terms
│   │   └── kegg_pathway_classes.json             # KEGG pathway taxonomy
│   └── outputs/                                  # All analysis outputs
│       ├── paper_figures/                        # Manuscript-ready figures
│       ├── sensitivity/                          # Parameter sensitivity results
│       │   ├── target_robustness.csv
│       │   ├── pathway_robustness.csv
│       │   ├── param_grid_results.csv
│       │   └── compound_sensitivity/             # LOO compound results
│       ├── shared/                               # Cross-stage shared data
│       └── traditional/                          # Baseline method outputs
│
├── docs/                                         # Project documentation
├── 分子对接/                                      # Docking panel composition assets
├── requirements.txt
├── pyproject.toml
├── reproduce.py                                  # Cross-platform reproduction orchestrator
├── reproduce.ps1                                 # Windows PowerShell wrapper
└── reproduce.sh                                  # Bash wrapper
```

---

## Output Files

### Core Results (Parameter Perturbation)

| File | Description |
|------|-------------|
| `src/outputs/sensitivity/target_robustness.csv` | Target frequency & degree statistics across 112 combos |
| `src/outputs/sensitivity/pathway_robustness.csv` | GO BP/KEGG functional-term robustness metrics |
| `src/outputs/sensitivity/param_grid_results.csv` | Per-combo enrichment summary |
| `src/outputs/sensitivity/pathway_rank_stability.csv` | Functional-term rank consistency across parameter space |

### Core Results (Compound Perturbation)

| File | Description |
|------|-------------|
| `src/outputs/sensitivity/compound_sensitivity/compound_influence_scores.csv` | Influence score per compound (0–1) |
| `src/outputs/sensitivity/compound_sensitivity/dual_robust_targets.csv` | 29 dual-robust targets with LOO retention rates |
| `src/outputs/sensitivity/compound_sensitivity/dual_robust_pathways.csv` | 58 dual-robust GO BP/KEGG functional terms with LOO retention rates |
| `src/outputs/sensitivity/compound_sensitivity/leave_one_out_pathway_robustness.csv` | Functional-term robustness under compound dropout |

### External Evidence and Manuscript Figures

Compatibility note: historical filenames and some CSV columns use `pathway_*` for continuity. In this repository, those fields refer to mixed **KEGG pathways plus GO BP biological-process terms**; GO BP terms are not pathways.

| File | Description |
|------|-------------|
| `src/outputs/sensitivity/robustness_validation/` | Top30 ablation, PubMed co-occurrence audit, and degree-matched permutation outputs |
| `src/outputs/paper_figures/` | Manuscript-ready figures including workflow, parameter heatmap, compound LOO heatmap, PubMed threshold audit, and docking panels |
| `src/outputs/paper_figures/figure-catalog.md` | Figure purpose, source data, and interpretation notes |
| `src/outputs/reproduction_manifest.json` | Runtime manifest generated by `reproduce.py` and ignored by Git |

### Checkpoints

- `src/outputs/sensitivity/compound_sensitivity/checkpoints/combo_results/` — 2,016 individual JSON files (one per compound × combo), enabling full audit and resume.

### Supplementary Tables

| File | Description |
|------|-------------|
| `src/outputs/supplementary/supplementary_tables_v1.xlsx` | Combined workbook for Tables S1-S6 |
| `src/outputs/supplementary/csv/` | CSV copy for each supplementary table |

Generation command:

```bash
python scripts/analysis/export_supplementary_tables.py
```

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

LC-MS/MS identified 19 raw records; after structure-normalized deduplication (Trehalose ≡ alpha,alpha-Trehalose), 18 compounds were retained:

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
15. Citric acid
16. A-Ketoglutaric Acid
17. L-Allo-Isoleucine
18. Hexadecanamide

---

## Target Classification

| Tier | Criteria | Count | Examples |
|------|----------|-------|----------|
| **Tier 1** | param_freq = 100% AND LOO_freq = 100% | 1 | **STAT3** |
| **Tier 2** | param_freq >= 50% and dual-robust | 13 | EGFR, HSP90AA1, TNF, BCL2, SRC, ... |
| **Tier 3** | param_freq < 50% but dual-robust | 15 | MAPK14, MET, IL2, KIT, CA1, ... |
| **Compound-sensitive** | High parameter frequency but lost after dropping key compound | 3 | CASP3, MAPK1, PIK3CA |
| **Peripheral** | Low parameter frequency and non-dual-robust | 18 | APP, MAPK8, BCL2L1, MTOR, ... |

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
- **Reproducibility package**: the repository clone includes packaged inputs, cached intermediate files, checkpoints, and manuscript-facing outputs required by `reproduce.py`.

External databases used:
- SwissTargetPrediction (http://www.swisstargetprediction.ch/)
- STRING v12 (https://string-db.org/)
- Enrichr (https://maayanlab.cloud/Enrichr/)
- Gene Ontology (http://purl.obolibrary.org/obo/go-basic.obo)

---

## Citation

If you use this code or framework, please cite:

> *Parameter and compound dual-perturbation robustness framework for network pharmacology target prioritization: strawberry leaf antioxidant compounds as a case study.*

---

## License

[MIT License](./LICENSE)
