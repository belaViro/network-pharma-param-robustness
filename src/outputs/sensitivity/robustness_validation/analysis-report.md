# Robustness Validation Bundle

## Analysis Question

Do parameter perturbations identify a more stable Top30 target set, and which compounds most perturb the antioxidant network under leave-one-compound-out analysis?

## Key Findings

- Dual-robust targets overlap OMIM-supported oxidative-stress genes at 2/29 targets (hypergeometric p=0.0525, BH-FDR=0.219).
- Q1 candidate robust targets overlap OMIM-supported genes at 3/50 targets (hypergeometric p=0.0178, BH-FDR=0.119).
- Dual-robust targets overlap PubMed-supported genes (PubMed_ge_10) at 21/29 targets against the candidate target background set (hypergeometric p=0.0656, BH-FDR=0.219).
- Q1 candidate robust targets overlap PubMed-supported genes (PubMed_ge_10) at 39/50 targets (hypergeometric p=0.000894, BH-FDR=0.00894).
- PubMed threshold sensitivity is reported across multiple cutoffs; at >=50 papers, dual-robust targets overlap 18/29 high-literature genes (BH-FDR=0.000274).
- Q1 parameter-robust targets show strong PubMed co-occurrence support at >=50 papers: 32/50 targets, enrichment=2.18, BH-FDR=3.26e-07.
- Dual-robust targets are not strongest at the permissive >=10 PubMed threshold, but at >=100 papers they overlap 17/29 targets (enrichment=3.1, BH-FDR=2.86e-06).
- Multi-threshold degree-matched PubMed check confirms dual-robust support at >=100 papers: observed=17, null mean=8.21, empirical p>=0.000999.
- Degree-matched null test for dual-robust OMIM support: observed=2, null mean=0.745, empirical p>=0.03 (n=1000, seed=20260503).
- Degree-matched null test for dual-robust PubMed support: observed=21, null mean=21.4, empirical p>=0.658 (n=1000, seed=20260503).
- Parameter Top30 ablation controls target-set size at n=30: Q1 robust Top30 has mean parameter frequency 0.543 and minimum frequency 0.268, compared with fixed default Top30 mean frequency 0.426 and minimum frequency 0.0357.
- Compound perturbation ablation ranks `(3beta,9xi)-3-(beta-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide` as the largest network driver (influence score=0.306; mean target retention=0.757; mean pathway retention=0.655).

## Output Tables

- `external_target_validation.csv`: source-support hypergeometric tests with BH-FDR against the candidate target background set.
- `degree_matched_support_metrics.csv`: observed and degree-matched null support metrics.
- `degree_matched_support_summary.csv`: empirical p-values and null intervals.
- `parameter_top30_ablation.csv`: fixed default, single-best parameter, frequency Top30, and Q1 robust Top30 comparison with n=30 controlled.
- `ablation_comparison.csv`: compatibility alias for `parameter_top30_ablation.csv`.
- `compound_perturbation_ablation.csv`: leave-one-compound-out retention and influence scores.
- `target_set_context_comparison.csv`: contextual comparison of fixed Top30, all Q1 targets, dual-robust targets, and Tier 1 target; not the primary ablation table.
- `target_pubmed_counts.csv`: target-level PubMed co-occurrence counts when requested or cached.
- `pubmed_threshold_sensitivity.csv`: PubMed support enrichment over multiple count thresholds.
- `pubmed_cooccurrence_key_results.csv`: manuscript-facing PubMed threshold table for Top30, Q1, and dual-robust target sets.
- `pubmed_dual_robust_target_counts.csv`: dual-robust target-level PubMed counts sorted by literature support.
- `pubmed_threshold_degree_matched_summary.csv`: optional multi-threshold degree-matched PubMed support checks when generated.

## PubMed Status

- Mode: `cache`; status: cache_hit=286.

## Boundary Conditions

- GeneCards, OMIM, and TTD are source-support checks, not fully independent validation, because these sources also participate in disease-gene construction.
- PubMed target co-occurrence is independent literature support; short symbols use gene/protein phrases and descriptions instead of bare-symbol matching, and the raw query is stored for audit.
- PubMed co-occurrence should be interpreted as literature support, not experimental validation; high-degree and high-attention hub genes can have inflated literature counts.
- Traditional Top30 remains strong in PubMed co-occurrence, so the PubMed result should support the robustness evidence chain rather than claim that dual-robust targets universally beat hub-based selection.
- Degree matching uses nearest-degree non-observed genes from the available PPI target background set; exact matching is impossible for the highest-degree hubs.
