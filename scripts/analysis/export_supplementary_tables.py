"""Export publication supplementary tables from finalized CSV outputs.

Usage:
    cd code
    python scripts/analysis/export_supplementary_tables.py
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.config import (  # noqa: E402
    STRAWBERRY_LEAF_RAW_COMPOUNDS,
    canonicalize_compound_name,
    get_strawberry_leaf_molecules,
    safe_compound_file_stem,
)


OUTPUTS_DIR = ROOT / "src" / "outputs"
SHARED_DIR = OUTPUTS_DIR / "shared"
SENSITIVITY_DIR = OUTPUTS_DIR / "sensitivity"
COMPOUND_SENSITIVITY_DIR = SENSITIVITY_DIR / "compound_sensitivity"
SUPPLEMENTARY_DIR = OUTPUTS_DIR / "supplementary"
CSV_DIR = SUPPLEMENTARY_DIR / "csv"
WORKBOOK_PATH = SUPPLEMENTARY_DIR / "supplementary_tables_v1.xlsx"

PROBABILITY_THRESHOLDS = [0.0, 0.05, 0.10]


def _normalize_gene(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def _split_gene_cells(series: pd.Series) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for gene_cell, probability in zip(series["gene_cell"], series["probability"]):
        if not gene_cell or str(gene_cell).lower() == "nan":
            continue
        for gene in str(gene_cell).split():
            gene = _normalize_gene(gene)
            if gene:
                rows.append({"gene": gene, "probability": probability})
    return pd.DataFrame(rows, columns=["gene", "probability"])


def build_compound_overview() -> pd.DataFrame:
    disease_path = SHARED_DIR / "disease_genes" / "merged_disease_genes.csv"
    disease_df = pd.read_csv(disease_path)
    disease_genes = {_normalize_gene(gene) for gene in disease_df["Gene_Symbol"]}
    disease_genes.discard("")

    raw_records_by_canonical: dict[str, list[str]] = {}
    for raw_name in STRAWBERRY_LEAF_RAW_COMPOUNDS:
        canonical = canonicalize_compound_name(raw_name)
        raw_records_by_canonical.setdefault(canonical, []).append(raw_name)

    rows: list[dict[str, Any]] = []
    for compound in get_strawberry_leaf_molecules():
        stem = safe_compound_file_stem(compound)
        target_path = SHARED_DIR / "swiss_targets" / f"{stem}_targets.csv"
        row: dict[str, Any] = {
            "compound_name": compound,
            "compound_stem": stem,
            "raw_lcms_records_collapsed": "; ".join(raw_records_by_canonical.get(compound, [compound])),
            "target_file": str(target_path.relative_to(ROOT)),
            "file_exists": target_path.exists(),
        }
        if not target_path.exists():
            rows.append(row)
            continue

        df = pd.read_csv(target_path)
        df.columns = df.columns.str.strip()
        gene_col = "Common name" if "Common name" in df.columns else "Target"
        prob_col = next((col for col in df.columns if "Probability" in col), None)

        probability = (
            pd.to_numeric(df[prob_col], errors="coerce").fillna(0.0)
            if prob_col
            else pd.Series([1.0] * len(df))
        )
        gene_cells = df[gene_col].astype(str).str.strip()
        split_df = _split_gene_cells(
            pd.DataFrame({"gene_cell": gene_cells, "probability": probability})
        )
        split_unique = split_df.drop_duplicates("gene") if not split_df.empty else split_df
        disease_split = split_df[split_df["gene"].isin(disease_genes)] if not split_df.empty else split_df

        row.update(
            {
                "raw_swisstarget_rows": len(df),
                "unique_gene_cells": int(gene_cells.str.upper().nunique()),
                "unique_split_genes": int(split_unique["gene"].nunique()) if not split_unique.empty else 0,
                "probability_min": float(probability.min()) if len(probability) else None,
                "probability_max": float(probability.max()) if len(probability) else None,
            }
        )
        for threshold in PROBABILITY_THRESHOLDS:
            suffix = str(threshold).replace(".", "_")
            row[f"n_predictions_probability_gt_{suffix}"] = int((split_df["probability"] > threshold).sum())
            row[f"n_disease_predictions_probability_gt_{suffix}"] = int(
                (disease_split["probability"] > threshold).sum()
            )

        if split_df.empty:
            row["top_predicted_genes"] = ""
        else:
            top_genes = (
                split_df.sort_values(["probability", "gene"], ascending=[False, True])
                .drop_duplicates("gene")
                .head(10)["gene"]
                .tolist()
            )
            row["top_predicted_genes"] = "; ".join(top_genes)
        rows.append(row)

    return pd.DataFrame(rows)


def _boolean_true(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1"])


def build_tables() -> dict[str, pd.DataFrame]:
    table_s1 = build_compound_overview()

    table_s2 = pd.read_csv(SENSITIVITY_DIR / "param_grid_results.csv")
    table_s2 = table_s2.sort_values(["prob_threshold", "string_threshold", "top_k"]).reset_index(drop=True)

    target_robustness = pd.read_csv(SENSITIVITY_DIR / "target_robustness.csv")
    table_s3 = target_robustness[
        target_robustness["quadrant"].astype(str).str.contains("Q1", na=False)
    ].copy()
    table_s3 = table_s3.sort_values(
        ["frequency", "avg_degree_norm", "avg_degree_rank", "avg_betweenness"],
        ascending=[False, False, True, False],
    ).reset_index(drop=True)
    table_s3.insert(0, "q1_rank", range(1, len(table_s3) + 1))

    dual_targets = pd.read_csv(COMPOUND_SENSITIVITY_DIR / "dual_robust_targets.csv")
    table_s4 = dual_targets[_boolean_true(dual_targets["is_dual_robust"])].copy()

    def tier(row: pd.Series) -> str:
        if row["full_frequency"] == 1.0 and row["loo_q1_retention_rate"] == 1.0:
            return "Tier 1"
        if row["full_frequency"] >= 0.5:
            return "Tier 2"
        return "Tier 3"

    table_s4.insert(0, "tier", table_s4.apply(tier, axis=1))
    table_s4 = table_s4.sort_values(
        ["tier", "loo_q1_retention_rate", "full_frequency", "full_avg_degree_norm"],
        ascending=[True, False, False, False],
    ).reset_index(drop=True)
    table_s4.insert(0, "dual_robust_rank", range(1, len(table_s4) + 1))

    table_s5 = pd.read_csv(COMPOUND_SENSITIVITY_DIR / "compound_influence_scores.csv")
    table_s5 = table_s5.sort_values("compound_influence_score", ascending=False).reset_index(drop=True)
    table_s5.insert(0, "influence_rank", range(1, len(table_s5) + 1))

    dual_pathways = pd.read_csv(COMPOUND_SENSITIVITY_DIR / "dual_robust_pathways.csv")
    table_s6 = dual_pathways[_boolean_true(dual_pathways["is_dual_robust"])].copy()
    table_s6.insert(
        0,
        "term_type",
        table_s6["source"].map({"KEGG": "KEGG pathway", "GO_BP": "GO biological process"}).fillna(
            table_s6["source"]
        ),
    )
    table_s6 = table_s6.sort_values(
        ["loo_q1_retention_rate", "full_frequency", "full_avg_logp", "full_avg_ef"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    table_s6.insert(0, "dual_robust_term_rank", range(1, len(table_s6) + 1))

    readme = pd.DataFrame(
        [
            {
                "sheet": "S1_compounds",
                "table": "Table S1",
                "description": "Structure-normalized strawberry leaf compounds and SwissTargetPrediction summary.",
                "source": "src/outputs/shared/swiss_targets/*_targets.csv",
                "rows": len(table_s1),
            },
            {
                "sheet": "S2_param_grid",
                "table": "Table S2",
                "description": "Results for the 112 parameter combinations.",
                "source": "src/outputs/sensitivity/param_grid_results.csv",
                "rows": len(table_s2),
            },
            {
                "sheet": "S3_Q1_targets",
                "table": "Table S3",
                "description": "Q1 candidate robust targets from target robustness analysis.",
                "source": "src/outputs/sensitivity/target_robustness.csv",
                "rows": len(table_s3),
            },
            {
                "sheet": "S4_dual_targets",
                "table": "Table S4",
                "description": "Dual-robust targets retained under parameter and compound perturbation.",
                "source": "src/outputs/sensitivity/compound_sensitivity/dual_robust_targets.csv",
                "rows": len(table_s4),
            },
            {
                "sheet": "S5_compound_influence",
                "table": "Table S5",
                "description": "Compound influence scores from leave-one-compound-out analysis.",
                "source": "src/outputs/sensitivity/compound_sensitivity/compound_influence_scores.csv",
                "rows": len(table_s5),
            },
            {
                "sheet": "S6_dual_terms",
                "table": "Table S6",
                "description": "Dual-robust functional enrichment terms.",
                "source": "src/outputs/sensitivity/compound_sensitivity/dual_robust_pathways.csv",
                "rows": len(table_s6),
            },
        ]
    )
    readme.insert(0, "generated_at", datetime.now().isoformat(timespec="seconds"))

    return {
        "README": readme,
        "S1_compounds": table_s1,
        "S2_param_grid": table_s2,
        "S3_Q1_targets": table_s3,
        "S4_dual_targets": table_s4,
        "S5_compound_influence": table_s5,
        "S6_dual_terms": table_s6,
    }


def write_workbook(tables: dict[str, pd.DataFrame]) -> None:
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)
    CSV_DIR.mkdir(parents=True, exist_ok=True)

    for sheet_name, df in tables.items():
        df.to_csv(CSV_DIR / f"{sheet_name}.csv", index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(WORKBOOK_PATH, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

        workbook = writer.book
        header_fill = PatternFill("solid", fgColor="1F4E78")
        header_font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
        body_font = Font(name="Arial", size=10)

        for worksheet in workbook.worksheets:
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for row in worksheet.iter_rows():
                for cell in row:
                    cell.font = body_font
                    cell.alignment = Alignment(vertical="top", wrap_text=True)
            for cell in worksheet[1]:
                cell.fill = header_fill
                cell.font = header_font
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

            for column_cells in worksheet.columns:
                column_letter = get_column_letter(column_cells[0].column)
                max_length = max(len(str(cell.value)) if cell.value is not None else 0 for cell in column_cells)
                worksheet.column_dimensions[column_letter].width = min(max(max_length + 2, 10), 55)


def main() -> None:
    tables = build_tables()
    write_workbook(tables)
    print(f"Supplementary workbook written: {WORKBOOK_PATH}")
    print(f"CSV copies written: {CSV_DIR}")
    for sheet_name, df in tables.items():
        print(f"  {sheet_name}: {len(df)} rows")


if __name__ == "__main__":
    main()
