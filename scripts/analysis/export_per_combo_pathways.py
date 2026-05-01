"""
Export per-combination enrichment details from the cached sensitivity run.

This script does not call Enrichr or any network service. It reconstructs the
pathway table from:
  - src/outputs/sensitivity/param_grid_results.csv
  - src/outputs/sensitivity/combo_targets.json
  - src/outputs/sensitivity/enrichment_cache/*.json
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


CODE_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CODE_DIR))

from src.sensitivity.classification import classify_pathway  # noqa: E402


SENSITIVITY_DIR = CODE_DIR / "src" / "outputs" / "sensitivity"
CACHE_DIR = SENSITIVITY_DIR / "enrichment_cache"
PARAM_GRID_CSV = SENSITIVITY_DIR / "param_grid_results.csv"
COMBO_TARGETS_JSON = SENSITIVITY_DIR / "combo_targets.json"
OUTPUT_CSV = SENSITIVITY_DIR / "per_combo_pathways_20260501.csv"
LATEST_CSV = SENSITIVITY_DIR / "per_combo_pathways_latest.csv"

INCLUDED_SOURCES = ("GO_BP", "KEGG")


def cache_key(gene_list: Iterable[str], pvalue_cutoff: float) -> str:
    genes_str = ";".join(sorted(gene_list))
    return hashlib.md5(f"{genes_str}|{pvalue_cutoff}".encode()).hexdigest()[:12]


def to_float(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def load_cache_for_combo(genes: list[str], pvalue_cutoff: float) -> tuple[str, dict[str, list[dict[str, Any]]]]:
    key = cache_key(genes, pvalue_cutoff)
    cache_path = CACHE_DIR / f"{key}.json"
    if not cache_path.exists():
        raise FileNotFoundError(f"Missing enrichment cache for key={key}: {cache_path}")
    with cache_path.open("r", encoding="utf-8") as f:
        return key, json.load(f)


def build_records() -> tuple[pd.DataFrame, dict[str, Any]]:
    param_df = pd.read_csv(PARAM_GRID_CSV)
    with COMBO_TARGETS_JSON.open("r", encoding="utf-8") as f:
        combo_targets: dict[str, list[str]] = json.load(f)

    records: list[dict[str, Any]] = []
    missing_combos: list[str] = []
    count_mismatches: list[dict[str, Any]] = []

    for _, combo_row in param_df.iterrows():
        combo = str(combo_row["combo"])
        genes = combo_targets.get(combo)
        if not genes:
            missing_combos.append(combo)
            continue

        pvalue_cutoff = to_float(combo_row["pvalue_cutoff"])
        cache_id, cached = load_cache_for_combo(genes, pvalue_cutoff)

        combo_records: list[dict[str, Any]] = []
        for source in INCLUDED_SOURCES:
            for item in cached.get(source, []):
                pathway = str(item.get("Term", ""))
                adjusted_pvalue = to_float(item.get("Adjusted P-value"), 1.0)
                combo_records.append(
                    {
                        "combo": combo,
                        "prob_threshold": to_float(combo_row["prob_threshold"]),
                        "string_threshold": to_float(combo_row["string_threshold"]),
                        "top_k": int(combo_row["top_k"]),
                        "pvalue_cutoff": pvalue_cutoff,
                        "cache_id": cache_id,
                        "n_targets": int(combo_row["n_targets"]),
                        "n_ppi_nodes": int(combo_row["n_ppi_nodes"]),
                        "n_ppi_edges": int(combo_row["n_ppi_edges"]),
                        "total_significant": int(combo_row["significant_pathways"]),
                        "source": source,
                        "gene_set": item.get("Gene_set", ""),
                        "pathway": pathway,
                        "pathway_category": classify_pathway(pathway, source),
                        "adjusted_pvalue": adjusted_pvalue,
                        "raw_pvalue": to_float(item.get("P-value")),
                        "ef": to_float(item.get("EF")),
                        "overlap": item.get("Overlap", ""),
                        "k": item.get("k", ""),
                        "N": item.get("N", ""),
                        "odds_ratio": to_float(item.get("Odds Ratio")),
                        "combined_score": to_float(item.get("Combined Score")),
                        "genes": item.get("Genes", ""),
                    }
                )

        combo_records.sort(key=lambda row: row["adjusted_pvalue"])
        for rank, record in enumerate(combo_records, 1):
            record["rank_in_combo"] = rank
            records.append(record)

        expected = int(combo_row["significant_pathways"])
        observed = len(combo_records)
        if observed != expected:
            count_mismatches.append(
                {
                    "combo": combo,
                    "expected": expected,
                    "observed": observed,
                }
            )

    if missing_combos:
        raise RuntimeError(f"Missing target lists for combos: {missing_combos[:10]}")
    if count_mismatches:
        raise RuntimeError(f"Pathway count mismatches: {count_mismatches[:10]}")

    pathway_df = pd.DataFrame(records)
    if not pathway_df.empty:
        ordered_cols = [
            "combo",
            "prob_threshold",
            "string_threshold",
            "top_k",
            "pvalue_cutoff",
            "rank_in_combo",
            "source",
            "pathway",
            "pathway_category",
            "adjusted_pvalue",
            "raw_pvalue",
            "ef",
            "overlap",
            "k",
            "N",
            "odds_ratio",
            "combined_score",
            "genes",
            "gene_set",
            "cache_id",
            "n_targets",
            "n_ppi_nodes",
            "n_ppi_edges",
            "total_significant",
        ]
        pathway_df = pathway_df[ordered_cols]

    summary = {
        "combinations": int(param_df["combo"].nunique()),
        "rows": int(len(pathway_df)),
        "expected_rows": int(param_df["significant_pathways"].sum()),
        "sources": list(INCLUDED_SOURCES),
        "output": str(OUTPUT_CSV),
        "latest_output": str(LATEST_CSV),
    }
    return pathway_df, summary


def main() -> None:
    pathway_df, summary = build_records()
    pathway_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    pathway_df.to_csv(LATEST_CSV, index=False, encoding="utf-8-sig")

    print("Exported per-combo pathway details")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
