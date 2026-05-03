# 分子对接 seed 稳定性验证报告

## 方法概览

- AutoDock Vina: v1.2.7
- exhaustiveness: 16
- seeds: 1001, 1002, 1003
- box size: 20 x 20 x 20 A
- core combinations: 3 receptors x 3 ligands = 9
- replicate dockings: 9 combinations x 3 seeds = 27

## 核心化合物重复对接汇总

| receptor | ligand | n | mean_affinity_kcal_mol | sd_affinity_kcal_mol | best_affinity_kcal_mol | stability |
| --- | --- | --- | --- | --- | --- | --- |
| SRC | Miquelianin | 3 | -9.8 | 0.024 | -9.825 | stable |
| EGFR | Miquelianin | 3 | -9.484 | 0.015 | -9.495 | stable |
| SRC | Catechin | 3 | -9.145 | 0.024 | -9.171 | stable |
| EGFR | Procyanidin_B1 | 3 | -8.497 | 0.016 | -8.514 | stable |
| SRC | Procyanidin_B1 | 3 | -8.32 | 0.016 | -8.331 | stable |
| EGFR | Catechin | 3 | -8.241 | 0.015 | -8.251 | stable |
| STAT3 | Miquelianin | 3 | -7.219 | 0.011 | -7.23 | stable |
| STAT3 | Procyanidin_B1 | 3 | -7.033 | 0.008 | -7.038 | stable |
| STAT3 | Catechin | 3 | -6.289 | 0.013 | -6.303 | stable |

## 判读标准

- sd < 0.3 kcal/mol: stable
- 0.3 <= sd <= 0.6 kcal/mol: acceptable
- sd > 0.6 kcal/mol: unstable

## 局限

- Vina scoring is approximate and cannot replace enzyme or cell experiments.
- Protein receptors are treated as rigid structures.
- Solvent and water-mediated interactions are not explicitly modeled.
- This report evaluates seed-level sampling stability only.
