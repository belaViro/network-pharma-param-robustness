# Result Lock v1.0

Generated on 2026-05-03.

This file records the current locked result scope for manuscript writing and supplementary-table preparation.

## Locked Analysis Scope

| Item | Locked value |
|---|---:|
| LC-MS/MS raw compound records | 19 |
| Canonical compounds after Trehalose deduplication | 18 |
| Parameter combinations | 112 |
| Compound leave-one-out runs | 2016 |
| Completed LOO checkpoints | 2016 |
| Skipped LOO combinations | 0 |

## Locked Target Counts

| Target universe | Count |
|---|---:|
| SwissTarget union targets | 684 |
| Disease-filtered predicted targets | 639 |
| Positive-probability disease targets | 286 |
| Main reporting threshold targets, Probability > 0.05 | 219 |
| Q1 candidate robust targets | 50 |
| Dual-robust targets | 29 |
| Tier 1 dual-robust target | 1 |

Tier 1 target: STAT3.

## Locked Functional Enrichment Counts

| Functional enrichment term set | Count |
|---|---:|
| Functional terms observed across the parameter grid | 213 |
| Q1 candidate robust functional terms | 70 |
| Dual-robust functional terms | 58 |

Use "functional enrichment terms" for the overall set. Use "KEGG pathways" for KEGG records and "GO biological processes" for GO_BP records.

## Locked Compound Influence Summary

| Rank | Compound | Influence score |
|---:|---|---:|
| 1 | Cardiac glycoside | 0.306 |
| 2 | Catechin | 0.201 |
| 3 | DL-Tryptophan | 0.107 |

## Supplementary Tables

| Table | File/sheet | Rows | Description |
|---|---|---:|---|
| Table S1 | `S1_compounds` | 18 | Canonical compounds and SwissTargetPrediction summary |
| Table S2 | `S2_param_grid` | 112 | Parameter-grid result summary |
| Table S3 | `S3_Q1_targets` | 50 | Q1 candidate robust targets |
| Table S4 | `S4_dual_targets` | 29 | Dual-robust targets |
| Table S5 | `S5_compound_influence` | 18 | Compound influence scores |
| Table S6 | `S6_dual_terms` | 58 | Dual-robust functional enrichment terms |

Workbook:

```text
src/outputs/supplementary/supplementary_tables_v1.xlsx
```

CSV copies:

```text
src/outputs/supplementary/csv/
```

Regeneration command:

```bash
python scripts/analysis/export_supplementary_tables.py
```
