# Stats Appendix

## PubMed Co-occurrence Validation

| target_set | threshold | overlap | expected | enrichment | p value | BH-FDR |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| q1_candidate_robust | 50 | 32 | 14.69 | 2.179 | 1.81e-08 | 3.26e-07 |
| dual_robust | 50 | 18 | 8.52 | 2.113 | 1.14e-04 | 2.74e-04 |
| dual_robust | 100 | 17 | 5.48 | 3.105 | 4.51e-07 | 2.86e-06 |
| dual_robust | 200 | 12 | 3.35 | 3.586 | 1.04e-05 | 3.73e-05 |

## PubMed Degree-Matched Threshold Validation

| target_set | threshold | observed | null mean | null 2.5%-97.5% | empirical p>= | z score |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| dual_robust | 50 | 18 | 13.75 | 10-18 | 0.036 | 2.05 |
| dual_robust | 100 | 17 | 8.21 | 5-12 | 0.001 | 4.70 |
| dual_robust | 200 | 12 | 5.75 | 3-9 | 0.001 | 3.84 |

## Docking Seed Stability

| receptor | ligand | n | mean_affinity_kcal_mol | sd_affinity_kcal_mol | stability |
| --- | --- | --- | --- | --- | --- |
| SRC | Miquelianin | 3 | -9.800 | 0.024 | stable |
| EGFR | Miquelianin | 3 | -9.484 | 0.015 | stable |
| SRC | Catechin | 3 | -9.145 | 0.024 | stable |
| EGFR | Procyanidin_B1 | 3 | -8.497 | 0.016 | stable |
| SRC | Procyanidin_B1 | 3 | -8.320 | 0.016 | stable |
| EGFR | Catechin | 3 | -8.241 | 0.015 | stable |
| STAT3 | Miquelianin | 3 | -7.219 | 0.011 | stable |
| STAT3 | Procyanidin_B1 | 3 | -7.033 | 0.008 | stable |
| STAT3 | Catechin | 3 | -6.289 | 0.013 | stable |

No inferential significance tests were run for docking, because seeds are repeated computational searches rather than independent biological samples.
