# Result Lock Robustness Validation Addendum

This addendum records validation artifacts generated after `result_lock_v1.md`.

## Validation Bundle

- Output directory: `src\outputs\sensitivity\robustness_validation`
- Degree-matched permutations: 1000
- Random seed: 20260503
- PubMed mode: `cache`

## Locked Robustness Validation Highlights

- Q1 targets OMIM overlap: 3/50, BH-FDR=0.1186.
- Dual-robust targets OMIM overlap: 2/29, BH-FDR=0.2186.
- Dual-robust targets PubMed overlap (PubMed_ge_10): 21/29, BH-FDR=0.2186.
- Dual-robust OMIM degree-matched empirical p>=: 0.02997.
- Dual-robust PubMed degree-matched empirical p>=: 0.6583.
- Q1 targets PubMed >=50 overlap: 32/50, BH-FDR=3.262e-07.
- Dual-robust targets PubMed >=100 overlap: 17/29, BH-FDR=2.861e-06.
- Dual-robust PubMed >=100 degree-matched empirical p>=: 0.000999.
- Dual-robust mean LOO Q1 retention: 0.977.
- Q1 robust Top30 mean parameter frequency: 0.5435; minimum frequency: 0.2679.
- Top compound perturbation driver: `(3beta,9xi)-3-(beta-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide`, influence score=0.3064.

## Files

- `external_target_validation.csv`
- `degree_matched_support_metrics.csv`
- `degree_matched_support_summary.csv`
- `degree_matched_random_sets.csv`
- `parameter_top30_ablation.csv`
- `ablation_comparison.csv`
- `compound_perturbation_ablation.csv`
- `target_set_context_comparison.csv`
- `target_pubmed_counts.csv`
- `pubmed_threshold_sensitivity.csv`
- `pubmed_cooccurrence_key_results.csv`
- `pubmed_dual_robust_target_counts.csv`
- `pubmed_threshold_degree_matched_summary.csv`
- `analysis-report.md`
- `stats-appendix.md`
- `figure-catalog.md`
