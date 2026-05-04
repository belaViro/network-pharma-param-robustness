# Analysis Report

## Main Questions

- Which targets and pathways remain robust under parameter and compound perturbation?
- Which compounds drive network-level instability?
- Are the selected docking pairs stable across fixed Vina seeds?

## Key Findings

- Top compound influence: Cardiac glycoside (score = 0.306).
- Q1 parameter-robust targets show strong PubMed co-occurrence support at >=50 papers: 32/50 targets, BH-FDR = 3.26e-07.
- Dual-robust targets are weak at PubMed >=10 but strong at stricter thresholds; at >=100 papers, 17/29 dual-robust targets overlap PubMed-supported genes, BH-FDR = 2.86e-06.
- Degree-matched PubMed threshold validation supports the dual-robust >=100 result: observed = 17, null mean = 8.21, empirical p>= = 0.001.
- Strongest docking pair: SRC x Miquelianin (mean affinity = -9.800 kcal/mol).
- Docking seed SD range: 0.008 to 0.024 kcal/mol across 9 receptor-ligand pairs.

## Caveats

- Docking replicates are computational seeds, not biological replicates.
- Vina affinity values are approximate scoring outputs and should not be interpreted as IC50/Kd.
- Functional enrichment visualizations summarize computational terms, not direct wet-lab pathway activation.
- PubMed co-occurrence is literature support, not experimental validation; traditional Top30 is also strong, so the PubMed result should be interpreted as support for a high-confidence core subset rather than proof that dual-robust targets dominate every external metric.
