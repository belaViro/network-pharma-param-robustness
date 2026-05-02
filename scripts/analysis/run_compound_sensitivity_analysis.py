"""
Leave-one-compound-out sensitivity analysis with resumable checkpoints.

This script extends the parameter robustness workflow to compound perturbation:
for each canonical strawberry leaf compound, remove that compound, rerun the
full parameter grid, and compare targets/pathways against the full-set
reference results.

Checkpoint policy:
  - Each leave-one-out x parameter-combo result is written to one JSON file.
  - Re-running the script skips completed "ok" and "skipped" checkpoints.
  - Failed/error checkpoints are retried on the next run.
  - Aggregate CSV/JSON outputs are rebuilt from checkpoint files at the end.

Usage:
    cd code
    python scripts/analysis/run_compound_sensitivity_analysis.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Any

import gseapy as gp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.config import (  # noqa: E402
    BACKGROUND_GENES,
    MERGED_DISEASE_GENES_CSV,
    PPI_EDGES_CSV,
    PROBABILITY_THRESHOLDS,
    SENSITIVITY_DIR,
    SENSITIVITY_PPI_EDGES_CSV,
    STRAWBERRY_LEAF_MOLECULES,
    SWISS_TARGETS_DIR,
    safe_compound_file_stem,
)
from src.sensitivity.analysis import analyze_enrichment  # noqa: E402
from src.sensitivity.classification import classify_pathway  # noqa: E402
from src.sensitivity.enrichment_cache import CACHE_DIR as GLOBAL_ENRICHMENT_CACHE_DIR  # noqa: E402
from src.sensitivity.pathway_robustness import classify_pathway_robustness  # noqa: E402


STRING_THRESHOLDS = [0.40, 0.50, 0.70, 0.90]
PVALUE_CUTOFF = 0.01
TOP_K_VALUES = [10, 20, 30, 50]
PARAM_COMBINATIONS = list(product(PROBABILITY_THRESHOLDS, STRING_THRESHOLDS, TOP_K_VALUES))

OUT_DIR = SENSITIVITY_DIR / "compound_sensitivity"
CHECKPOINT_DIR = OUT_DIR / "checkpoints"
COMBO_CHECKPOINT_DIR = CHECKPOINT_DIR / "combo_results"
LOCAL_ENRICHMENT_CACHE_DIR = OUT_DIR / "enrichment_cache"
FIGURE_DIR = OUT_DIR / "figures"

LIBRARIES = {
    "GO_BP": "GO_Biological_Process_2021",
    "KEGG": "KEGG_2021_Human",
}


def json_safe(value: Any) -> Any:
    """Convert pandas/numpy values into JSON-serializable Python values."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, tuple):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
        return None if math.isnan(value) or math.isinf(value) else value
    if isinstance(value, float):
        return None if math.isnan(value) or math.isinf(value) else value
    if pd.isna(value):
        return None
    return value


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON via temp file + replace so interrupted writes do not corrupt checkpoints."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "w", encoding="utf-8") as fh:
        json.dump(json_safe(payload), fh, ensure_ascii=False, indent=2)
    os.replace(tmp_path, path)


def read_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def combo_id(prob_th: float, string_th: float, top_k: int) -> str:
    return f"P{prob_th}_S{string_th}_K{top_k}"


def checkpoint_path(drop_stem: str, combo: str) -> Path:
    raw = f"{drop_stem}__{combo}"
    digest = hashlib.md5(raw.encode("utf-8")).hexdigest()[:10]
    short_stem = drop_stem[:42]
    return COMBO_CHECKPOINT_DIR / f"{short_stem}__{combo}__{digest}.json"


def is_complete_checkpoint(payload: dict[str, Any]) -> bool:
    return payload.get("status") in {"ok", "skipped"}


def load_compound_target_records() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load all per-compound SwissTarget CSV files before deduplication."""
    rows: list[dict[str, Any]] = []
    compound_rows: list[dict[str, Any]] = []

    for compound_name in STRAWBERRY_LEAF_MOLECULES:
        stem = safe_compound_file_stem(compound_name)
        path = SWISS_TARGETS_DIR / f"{stem}_targets.csv"
        compound_rows.append(
            {
                "compound_name": compound_name,
                "compound_stem": stem,
                "target_file": str(path),
                "file_exists": path.exists(),
            }
        )
        if not path.exists():
            continue

        df = pd.read_csv(path)
        df.columns = df.columns.str.strip()

        if "Common name" in df.columns:
            gene_series = df["Common name"].astype(str).str.strip()
        elif "Target" in df.columns:
            gene_series = df["Target"].astype(str).str.strip()
        else:
            gene_series = df.iloc[:, 0].astype(str).str.strip()

        prob_col = next((col for col in df.columns if "Probability" in col), None)
        if prob_col is None:
            prob_series = pd.Series([1.0] * len(df))
        else:
            prob_series = pd.to_numeric(df[prob_col], errors="coerce").fillna(0.0)

        for raw_gene, probability in zip(gene_series.tolist(), prob_series.tolist()):
            if not raw_gene or raw_gene.lower() == "nan":
                continue
            # SwissTarget can emit multi-gene cells such as "FNTA FNTB".
            for gene in str(raw_gene).split():
                gene = gene.strip().upper()
                if not gene:
                    continue
                rows.append(
                    {
                        "Target_Gene": gene,
                        "Source_Compound": stem,
                        "Compound_Name": compound_name,
                        "Probability": float(probability),
                    }
                )

    raw_df = pd.DataFrame(rows)
    compounds_df = pd.DataFrame(compound_rows)
    if not raw_df.empty:
        counts = raw_df.groupby("Source_Compound")["Target_Gene"].count().to_dict()
        compounds_df["raw_target_records"] = compounds_df["compound_stem"].map(counts).fillna(0).astype(int)
    else:
        compounds_df["raw_target_records"] = 0
    return raw_df, compounds_df


def deduplicate_targets(
    raw_targets_df: pd.DataFrame,
    disease_genes: set[str],
    excluded_compound: str | None = None,
) -> pd.DataFrame:
    """Rebuild the deduplicated target table after optionally dropping one compound."""
    df = raw_targets_df.copy()
    if excluded_compound is not None:
        df = df[df["Source_Compound"] != excluded_compound].copy()
    df = df[df["Target_Gene"].isin(disease_genes)].copy()
    if df.empty:
        return pd.DataFrame(columns=["Target_Gene", "Source_Compound", "Compound_Name", "Probability"])

    df = df.sort_values(["Probability", "Target_Gene"], ascending=[False, True])
    idx = df.groupby("Target_Gene")["Probability"].idxmax()
    dedup = df.loc[idx].sort_values("Probability", ascending=False).reset_index(drop=True)
    return dedup[["Target_Gene", "Source_Compound", "Compound_Name", "Probability"]]


def _local_enrichment_cache_key(genes: list[str], pvalue_cutoff: float, top_n: int) -> str:
    genes_str = ";".join(sorted(genes))
    payload = f"{genes_str}|{pvalue_cutoff}|{top_n}|GO_BP_KEGG"
    return hashlib.md5(payload.encode("utf-8")).hexdigest()[:16]


def _global_enrichment_cache_key(genes: list[str], pvalue_cutoff: float) -> str:
    genes_str = ";".join(sorted(genes))
    return hashlib.md5(f"{genes_str}|{pvalue_cutoff}".encode()).hexdigest()[:12]


def _load_enrichment_json(path: Path) -> dict[str, pd.DataFrame] | None:
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as fh:
        cached = json.load(fh)
    return {source: pd.DataFrame(rows) for source, rows in cached.items()}


def _save_enrichment_json(path: Path, significant: dict[str, pd.DataFrame]) -> None:
    serializable = {
        source: df.to_dict(orient="records") if isinstance(df, pd.DataFrame) else []
        for source, df in significant.items()
    }
    atomic_write_json(path, serializable)


def cached_enrichment_go_bp_kegg(
    genes: list[str],
    pvalue_cutoff: float,
    background: int,
    top_n: int,
    sleep_seconds: float,
    retries: int,
    retry_sleep_seconds: float,
) -> dict[str, Any]:
    """Run GO_BP + KEGG enrichment with local cache and fallback to existing global cache."""
    genes = sorted(set(str(g).upper() for g in genes))
    local_key = _local_enrichment_cache_key(genes, pvalue_cutoff, top_n)
    local_path = LOCAL_ENRICHMENT_CACHE_DIR / f"{local_key}.json"
    cached = _load_enrichment_json(local_path)
    if cached is not None:
        return {"status": "成功", "significant": cached, "cache_id": local_key, "cache_source": "compound"}

    global_key = _global_enrichment_cache_key(genes, pvalue_cutoff)
    global_path = GLOBAL_ENRICHMENT_CACHE_DIR / f"{global_key}.json"
    global_cached = _load_enrichment_json(global_path)
    if global_cached is not None:
        significant = {
            source: global_cached.get(source, pd.DataFrame())
            for source in LIBRARIES
        }
        _save_enrichment_json(local_path, significant)
        return {"status": "成功", "significant": significant, "cache_id": local_key, "cache_source": "global"}

    LOCAL_ENRICHMENT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    significant: dict[str, pd.DataFrame] = {}
    n_genes = len(genes)

    for source, library in LIBRARIES.items():
        enr = None
        last_error: Exception | None = None
        for attempt in range(1, retries + 2):
            try:
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=library,
                    organism="hsapiens",
                    cutoff=pvalue_cutoff,
                    outdir=None,
                    background=background,
                )
                break
            except Exception as exc:
                last_error = exc
                print(
                    f"[WARN] {source} enrichment failed "
                    f"(attempt {attempt}/{retries + 1}): {exc}",
                    flush=True,
                )
                if attempt <= retries:
                    time.sleep(retry_sleep_seconds)

        if enr is None:
            raise RuntimeError(f"{source} enrichment failed after retries: {last_error}")

        if sleep_seconds > 0:
            time.sleep(sleep_seconds)

        if enr.res2d is None or len(enr.res2d) == 0:
            significant[source] = pd.DataFrame()
            continue

        df = enr.res2d[enr.res2d["Adjusted P-value"] < pvalue_cutoff].copy()
        if df.empty:
            significant[source] = pd.DataFrame()
            continue

        df["k"] = df["Overlap"].apply(lambda x: int(str(x).split("/")[0]))
        df["N"] = df["Overlap"].apply(lambda x: int(str(x).split("/")[1]))
        df["EF"] = (df["k"] / n_genes) / (df["N"] / background)
        df = df[(df["k"] >= 2) & (df["EF"] >= 3.0)]
        df = df.sort_values("Adjusted P-value").head(top_n)
        significant[source] = df

    _save_enrichment_json(local_path, significant)
    return {"status": "成功", "significant": significant, "cache_id": local_key, "cache_source": "new"}


def build_string_graphs(string_edges: pd.DataFrame) -> dict[float, nx.Graph]:
    graphs: dict[float, nx.Graph] = {}
    for string_th in STRING_THRESHOLDS:
        filtered_edges = string_edges[string_edges["combined_score"] >= string_th]
        graph = nx.Graph()
        for _, row in filtered_edges.iterrows():
            graph.add_edge(str(row["protein1"]).upper(), str(row["protein2"]).upper(), weight=row["combined_score"])
        graphs[string_th] = graph
    return graphs


def run_single_combo(
    filtered_targets: list[str],
    string_graph: nx.Graph,
    string_th: float,
    pvalue_cutoff: float,
    top_k: int,
    sleep_seconds: float,
    retries: int,
    retry_sleep_seconds: float,
) -> dict[str, Any]:
    n_targets = len(filtered_targets)
    if n_targets < 10:
        return {"status": "skipped", "skip_reason": "fewer_than_10_targets", "n_targets": n_targets}

    valid_nodes = [gene for gene in filtered_targets if gene in string_graph]
    graph = string_graph.subgraph(valid_nodes).copy()
    if graph.number_of_nodes() < 10:
        return {
            "status": "skipped",
            "skip_reason": "fewer_than_10_ppi_nodes",
            "n_targets": n_targets,
            "n_ppi_nodes": graph.number_of_nodes(),
            "n_ppi_edges": graph.number_of_edges(),
        }

    degree_sorted = sorted(graph.degree(), key=lambda x: x[1], reverse=True)
    topk_pairs = degree_sorted[:top_k]
    topk = [gene for gene, _ in topk_pairs]
    if len(topk) < 5:
        return {
            "status": "skipped",
            "skip_reason": "fewer_than_5_topk_genes",
            "n_targets": n_targets,
            "n_ppi_nodes": graph.number_of_nodes(),
            "n_ppi_edges": graph.number_of_edges(),
        }

    betweenness = nx.betweenness_centrality(graph)
    topk_info = [
        {
            "gene": gene,
            "degree": int(degree),
            "degree_rank": rank,
            "betweenness": float(betweenness.get(gene, 0.0)),
        }
        for rank, (gene, degree) in enumerate(topk_pairs, start=1)
    ]

    enrichment = cached_enrichment_go_bp_kegg(
        genes=topk,
        pvalue_cutoff=pvalue_cutoff,
        background=BACKGROUND_GENES,
        top_n=20,
        sleep_seconds=sleep_seconds,
        retries=retries,
        retry_sleep_seconds=retry_sleep_seconds,
    )
    metrics = analyze_enrichment(enrichment, include_sources=["KEGG", "GO_BP"])

    return {
        "status": "ok",
        "n_targets": n_targets,
        "n_ppi_nodes": graph.number_of_nodes(),
        "n_ppi_edges": graph.number_of_edges(),
        "top_k": len(topk),
        "topk_genes": topk,
        "topk_info": topk_info,
        "significant_pathways": metrics["total_significant"],
        "direct_ox_count": metrics["direct_ox"],
        "cross_talk_count": metrics["cross_talk"],
        "backbone_count": metrics["backbone"],
        "other_count": metrics["other"],
        "avg_ef": metrics["avg_ef"],
        "min_pvalue": metrics["min_pvalue"],
        "pathway_details": metrics["pathway_details"],
        "enrichment_cache_id": enrichment.get("cache_id"),
        "enrichment_cache_source": enrichment.get("cache_source"),
    }


def payload_to_result_row(payload: dict[str, Any], full_refs: dict[str, Any]) -> dict[str, Any]:
    combo = payload["combo"]
    topk_set = set(payload.get("topk_genes", []))
    pathway_set = {p["pathway"] for p in payload.get("pathway_details", [])}

    full_targets = set(full_refs["combo_targets"].get(combo, []))
    full_pathways = set(full_refs["combo_pathways"].get(combo, []))

    target_retention = len(topk_set & full_targets) / len(full_targets) if full_targets else np.nan
    pathway_retention = len(pathway_set & full_pathways) / len(full_pathways) if full_pathways else np.nan

    return {
        "drop_compound": payload["drop_compound"],
        "drop_compound_name": payload["drop_compound_name"],
        "combo": combo,
        "prob_threshold": payload["prob_threshold"],
        "string_threshold": payload["string_threshold"],
        "pvalue_cutoff": payload["pvalue_cutoff"],
        "top_k": payload["top_k_param"],
        "status": payload["status"],
        "skip_reason": payload.get("skip_reason", ""),
        "n_targets": payload.get("n_targets"),
        "n_ppi_nodes": payload.get("n_ppi_nodes"),
        "n_ppi_edges": payload.get("n_ppi_edges"),
        "significant_pathways": payload.get("significant_pathways"),
        "direct_ox_count": payload.get("direct_ox_count"),
        "cross_talk_count": payload.get("cross_talk_count"),
        "backbone_count": payload.get("backbone_count"),
        "other_count": payload.get("other_count"),
        "avg_ef": payload.get("avg_ef"),
        "min_pvalue": payload.get("min_pvalue"),
        "target_retention_vs_full_combo": target_retention,
        "pathway_retention_vs_full_combo": pathway_retention,
        "enrichment_cache_source": payload.get("enrichment_cache_source", ""),
        "checkpoint_file": payload.get("checkpoint_file", ""),
    }


def load_full_references() -> dict[str, Any]:
    with open(SENSITIVITY_DIR / "combo_targets.json", "r", encoding="utf-8") as fh:
        combo_targets = json.load(fh)

    combo_pathways_df = pd.read_csv(SENSITIVITY_DIR / "per_combo_pathways_latest.csv")
    combo_pathways = {
        combo: set(group["pathway"].astype(str).tolist())
        for combo, group in combo_pathways_df.groupby("combo")
    }

    target_robustness = pd.read_csv(SENSITIVITY_DIR / "target_robustness.csv")
    full_q1_targets = set(
        target_robustness[
            target_robustness["quadrant"].astype(str).str.contains("Q1", na=False)
        ]["gene"].astype(str)
    )

    pathway_robustness = pd.read_csv(SENSITIVITY_DIR / "pathway_robustness.csv")
    full_q1_pathways = set(
        pathway_robustness[
            pathway_robustness["is_candidate_robust"].astype(str).str.lower().isin(["true", "1"])
        ]["pathway"].astype(str)
    )

    return {
        "combo_targets": combo_targets,
        "combo_pathways": combo_pathways,
        "target_robustness": target_robustness,
        "pathway_robustness": pathway_robustness,
        "full_q1_targets": full_q1_targets,
        "full_q1_pathways": full_q1_pathways,
    }


def compute_target_robustness_by_drop(
    target_long_df: pd.DataFrame,
    result_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    ok_counts = result_df[result_df["status"] == "ok"].groupby("drop_compound")["combo"].nunique().to_dict()

    for drop_compound, group in target_long_df.groupby("drop_compound"):
        n_total = ok_counts.get(drop_compound, 0)
        if n_total == 0 or group.empty:
            continue

        group = group.copy()
        combo_max_degree = group.groupby("combo")["degree"].max().to_dict()
        group["degree_norm"] = group.apply(
            lambda row: row["degree"] / combo_max_degree[row["combo"]] if combo_max_degree[row["combo"]] else 0.0,
            axis=1,
        )

        stats: list[dict[str, Any]] = []
        for gene, gene_group in group.groupby("gene"):
            stats.append(
                {
                    "drop_compound": drop_compound,
                    "drop_compound_name": gene_group["drop_compound_name"].iloc[0],
                    "gene": gene,
                    "frequency": len(gene_group) / n_total,
                    "count": len(gene_group),
                    "avg_degree": gene_group["degree"].mean(),
                    "avg_degree_norm": gene_group["degree_norm"].mean(),
                    "avg_degree_rank": gene_group["degree_rank"].mean(),
                    "avg_betweenness": gene_group["betweenness"].mean(),
                    "min_rank": gene_group["degree_rank"].min(),
                    "max_rank": gene_group["degree_rank"].max(),
                }
            )

        stats_df = pd.DataFrame(stats)
        if stats_df.empty:
            continue
        freq_med = stats_df["frequency"].median()
        degree_med = stats_df["avg_degree_norm"].median()

        def classify(row: pd.Series) -> str:
            high_freq = row["frequency"] >= freq_med
            high_degree = row["avg_degree_norm"] >= degree_med
            if high_freq and high_degree:
                return "Q1: candidate_robust"
            if not high_freq and high_degree:
                return "Q2: low_frequency_high_importance"
            if high_freq and not high_degree:
                return "Q3: stable_but_peripheral"
            return "Q4: low_priority"

        stats_df["quadrant"] = stats_df.apply(classify, axis=1)
        rows.extend(stats_df.to_dict(orient="records"))

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows).sort_values(
        ["drop_compound", "frequency", "avg_degree_norm", "avg_degree_rank", "avg_betweenness"],
        ascending=[True, False, False, True, False],
    )


def compute_pathway_robustness_by_drop(
    pathway_long_df: pd.DataFrame,
    result_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    ok_counts = result_df[result_df["status"] == "ok"].groupby("drop_compound")["combo"].nunique().to_dict()

    for drop_compound, group in pathway_long_df.groupby("drop_compound"):
        n_total = ok_counts.get(drop_compound, 0)
        if n_total == 0 or group.empty:
            continue

        rank_rows: list[dict[str, Any]] = []
        for combo, combo_group in group.groupby("combo"):
            sorted_group = combo_group.sort_values("pvalue").reset_index(drop=True)
            for rank, (_, row) in enumerate(sorted_group.iterrows(), start=1):
                pvalue = float(row["pvalue"])
                rank_rows.append(
                    {
                        "pathway": row["pathway"],
                        "source": row["source"],
                        "combo": combo,
                        "rank": rank,
                        "pvalue": pvalue,
                        "ef": row["ef"],
                        "logp": -np.log10(max(pvalue, 1e-300)),
                        "pathway_category": row["pathway_category"],
                    }
                )

        rank_df = pd.DataFrame(rank_rows)
        if rank_df.empty:
            continue

        stats: list[dict[str, Any]] = []
        for pathway, pathway_group in rank_df.groupby("pathway"):
            category_counts = pathway_group["pathway_category"].value_counts()
            stats.append(
                {
                    "drop_compound": drop_compound,
                    "pathway": pathway,
                    "source": pathway_group["source"].iloc[0],
                    "pathway_category": category_counts.index[0],
                    "frequency": len(pathway_group) / n_total,
                    "count": len(pathway_group),
                    "avg_rank": pathway_group["rank"].mean(),
                    "min_rank": pathway_group["rank"].min(),
                    "max_rank": pathway_group["rank"].max(),
                    "rank_std": pathway_group["rank"].std(),
                    "avg_logp": pathway_group["logp"].mean(),
                    "avg_ef": pathway_group["ef"].mean(),
                }
            )

        stats_df = pd.DataFrame(stats)
        robust_df, _summary = classify_pathway_robustness(stats_df)
        robust_df["drop_compound"] = drop_compound
        rows.extend(robust_df.to_dict(orient="records"))

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def compute_dual_robust_targets(
    full_refs: dict[str, Any],
    target_robustness_df: pd.DataFrame,
    drop_compounds: list[str],
    threshold: float,
) -> pd.DataFrame:
    full_q1 = full_refs["target_robustness"][
        full_refs["target_robustness"]["quadrant"].astype(str).str.contains("Q1", na=False)
    ].copy()

    q1_by_drop = {
        drop: set(group[group["quadrant"].astype(str).str.contains("Q1", na=False)]["gene"].astype(str))
        for drop, group in target_robustness_df.groupby("drop_compound")
    }
    freq_by_drop_gene = {
        (row["drop_compound"], row["gene"]): row["frequency"]
        for _, row in target_robustness_df.iterrows()
    }

    rows: list[dict[str, Any]] = []
    for _, row in full_q1.iterrows():
        gene = str(row["gene"])
        retained = [drop for drop in drop_compounds if gene in q1_by_drop.get(drop, set())]
        missing = [drop for drop in drop_compounds if gene not in q1_by_drop.get(drop, set())]
        frequencies = [float(freq_by_drop_gene.get((drop, gene), 0.0)) for drop in drop_compounds]
        retention_rate = len(retained) / len(drop_compounds) if drop_compounds else 0.0
        rows.append(
            {
                "gene": gene,
                "is_dual_robust": retention_rate >= threshold and min(frequencies) > 0,
                "loo_q1_retention_rate": retention_rate,
                "loo_q1_retention_count": len(retained),
                "n_leave_one_out": len(drop_compounds),
                "min_loo_frequency": min(frequencies) if frequencies else 0.0,
                "mean_loo_frequency": float(np.mean(frequencies)) if frequencies else 0.0,
                "missing_q1_when_dropping": ";".join(missing),
                "full_frequency": row.get("frequency"),
                "full_avg_degree_norm": row.get("avg_degree_norm"),
                "full_avg_degree_rank": row.get("avg_degree_rank"),
                "full_avg_betweenness": row.get("avg_betweenness"),
                "full_quadrant": row.get("quadrant"),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["is_dual_robust", "loo_q1_retention_rate", "min_loo_frequency", "full_frequency", "full_avg_degree_norm"],
        ascending=[False, False, False, False, False],
    )


def compute_dual_robust_pathways(
    full_refs: dict[str, Any],
    pathway_robustness_df: pd.DataFrame,
    drop_compounds: list[str],
    threshold: float,
) -> pd.DataFrame:
    full_q1 = full_refs["pathway_robustness"][
        full_refs["pathway_robustness"]["is_candidate_robust"].astype(str).str.lower().isin(["true", "1"])
    ].copy()

    q1_by_drop = {
        drop: set(group[group["is_candidate_robust"].astype(str).str.lower().isin(["true", "1"])]["pathway"].astype(str))
        for drop, group in pathway_robustness_df.groupby("drop_compound")
    }
    freq_by_drop_pathway = {
        (row["drop_compound"], row["pathway"]): row["frequency"]
        for _, row in pathway_robustness_df.iterrows()
    }

    rows: list[dict[str, Any]] = []
    for _, row in full_q1.iterrows():
        pathway = str(row["pathway"])
        retained = [drop for drop in drop_compounds if pathway in q1_by_drop.get(drop, set())]
        missing = [drop for drop in drop_compounds if pathway not in q1_by_drop.get(drop, set())]
        frequencies = [float(freq_by_drop_pathway.get((drop, pathway), 0.0)) for drop in drop_compounds]
        retention_rate = len(retained) / len(drop_compounds) if drop_compounds else 0.0
        rows.append(
            {
                "pathway": pathway,
                "source": row.get("source"),
                "pathway_category": row.get("pathway_category"),
                "is_dual_robust": retention_rate >= threshold and min(frequencies) > 0,
                "loo_q1_retention_rate": retention_rate,
                "loo_q1_retention_count": len(retained),
                "n_leave_one_out": len(drop_compounds),
                "min_loo_frequency": min(frequencies) if frequencies else 0.0,
                "mean_loo_frequency": float(np.mean(frequencies)) if frequencies else 0.0,
                "missing_q1_when_dropping": ";".join(missing),
                "full_frequency": row.get("frequency"),
                "full_avg_logp": row.get("avg_logp"),
                "full_avg_rank": row.get("avg_rank"),
                "full_avg_ef": row.get("avg_ef"),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["is_dual_robust", "loo_q1_retention_rate", "min_loo_frequency", "full_frequency", "full_avg_logp"],
        ascending=[False, False, False, False, False],
    )


def compute_compound_influence_scores(
    result_df: pd.DataFrame,
    target_robustness_df: pd.DataFrame,
    pathway_robustness_df: pd.DataFrame,
    full_refs: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    full_q1_targets = full_refs["full_q1_targets"]
    full_q1_pathways = full_refs["full_q1_pathways"]

    q1_targets_by_drop = {
        drop: set(group[group["quadrant"].astype(str).str.contains("Q1", na=False)]["gene"].astype(str))
        for drop, group in target_robustness_df.groupby("drop_compound")
    }
    q1_pathways_by_drop = {
        drop: set(group[group["is_candidate_robust"].astype(str).str.lower().isin(["true", "1"])]["pathway"].astype(str))
        for drop, group in pathway_robustness_df.groupby("drop_compound")
    }

    for drop, group in result_df.groupby("drop_compound"):
        ok = group[group["status"] == "ok"]
        q1_target_retention = (
            len(q1_targets_by_drop.get(drop, set()) & full_q1_targets) / len(full_q1_targets)
            if full_q1_targets else np.nan
        )
        q1_pathway_retention = (
            len(q1_pathways_by_drop.get(drop, set()) & full_q1_pathways) / len(full_q1_pathways)
            if full_q1_pathways else np.nan
        )
        target_retention = float(ok["target_retention_vs_full_combo"].mean()) if not ok.empty else np.nan
        pathway_retention = float(ok["pathway_retention_vs_full_combo"].mean()) if not ok.empty else np.nan
        retention_values = [
            value for value in [target_retention, pathway_retention, q1_target_retention, q1_pathway_retention]
            if not pd.isna(value)
        ]
        mean_retention = float(np.mean(retention_values)) if retention_values else np.nan

        rows.append(
            {
                "drop_compound": drop,
                "drop_compound_name": group["drop_compound_name"].iloc[0],
                "ok_combinations": len(ok),
                "skipped_combinations": int((group["status"] == "skipped").sum()),
                "mean_target_retention_vs_full_combo": target_retention,
                "mean_pathway_retention_vs_full_combo": pathway_retention,
                "q1_target_retention_vs_full": q1_target_retention,
                "q1_pathway_retention_vs_full": q1_pathway_retention,
                "mean_retention": mean_retention,
                "compound_influence_score": 1.0 - mean_retention if not pd.isna(mean_retention) else np.nan,
            }
        )

    return pd.DataFrame(rows).sort_values("compound_influence_score", ascending=False)


def plot_outputs(influence_df: pd.DataFrame, dual_targets_df: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    if not influence_df.empty:
        plot_df = influence_df.sort_values("compound_influence_score", ascending=True)
        plt.figure(figsize=(9, max(4, len(plot_df) * 0.35)))
        plt.barh(plot_df["drop_compound"], plot_df["compound_influence_score"], color="#4C78A8")
        plt.xlabel("Compound influence score")
        plt.ylabel("Dropped compound")
        plt.tight_layout()
        plt.savefig(FIGURE_DIR / "fig_compound_influence.png", dpi=300)
        plt.savefig(FIGURE_DIR / "fig_compound_influence.svg")
        plt.close()

    if not dual_targets_df.empty:
        plot_df = dual_targets_df.head(30).sort_values("loo_q1_retention_rate", ascending=True)
        colors = ["#59A14F" if v else "#E15759" for v in plot_df["is_dual_robust"]]
        plt.figure(figsize=(8, max(4, len(plot_df) * 0.28)))
        plt.barh(plot_df["gene"], plot_df["loo_q1_retention_rate"], color=colors)
        plt.xlim(0, 1.02)
        plt.xlabel("Leave-one-compound-out Q1 retention rate")
        plt.ylabel("Full-set Q1 target")
        plt.tight_layout()
        plt.savefig(FIGURE_DIR / "fig_dual_robust_targets.png", dpi=300)
        plt.savefig(FIGURE_DIR / "fig_dual_robust_targets.svg")
        plt.close()


def aggregate_outputs(
    expected_payloads: list[dict[str, Any]],
    compounds_df: pd.DataFrame,
    full_refs: dict[str, Any],
    dual_retention_threshold: float,
) -> dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    result_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    pathway_rows: list[dict[str, Any]] = []

    for payload in expected_payloads:
        if payload.get("status") not in {"ok", "skipped"}:
            continue
        result_rows.append(payload_to_result_row(payload, full_refs))

        if payload.get("status") != "ok":
            continue

        for tinfo in payload.get("topk_info", []):
            target_rows.append(
                {
                    **tinfo,
                    "drop_compound": payload["drop_compound"],
                    "drop_compound_name": payload["drop_compound_name"],
                    "combo": payload["combo"],
                    "prob_threshold": payload["prob_threshold"],
                    "string_threshold": payload["string_threshold"],
                    "pvalue_cutoff": payload["pvalue_cutoff"],
                    "top_k": payload["top_k_param"],
                }
            )

        for pdet in payload.get("pathway_details", []):
            pathway_rows.append(
                {
                    **pdet,
                    "drop_compound": payload["drop_compound"],
                    "drop_compound_name": payload["drop_compound_name"],
                    "combo": payload["combo"],
                    "prob_threshold": payload["prob_threshold"],
                    "string_threshold": payload["string_threshold"],
                    "pvalue_cutoff": payload["pvalue_cutoff"],
                    "top_k": payload["top_k_param"],
                }
            )

    result_df = pd.DataFrame(result_rows)
    target_long_df = pd.DataFrame(target_rows)
    pathway_long_df = pd.DataFrame(pathway_rows)

    result_df.to_csv(OUT_DIR / "compound_sensitivity_results.csv", index=False)
    target_long_df.to_csv(OUT_DIR / "compound_targets_long.csv", index=False)
    pathway_long_df.to_csv(OUT_DIR / "compound_pathways_long.csv", index=False)
    compounds_df.to_csv(OUT_DIR / "compound_input_manifest.csv", index=False)

    target_robustness_df = compute_target_robustness_by_drop(target_long_df, result_df)
    pathway_robustness_df = compute_pathway_robustness_by_drop(pathway_long_df, result_df)
    target_robustness_df.to_csv(OUT_DIR / "leave_one_out_target_robustness.csv", index=False)
    pathway_robustness_df.to_csv(OUT_DIR / "leave_one_out_pathway_robustness.csv", index=False)

    drop_compounds = compounds_df["compound_stem"].tolist()
    dual_targets_df = compute_dual_robust_targets(
        full_refs,
        target_robustness_df,
        drop_compounds,
        dual_retention_threshold,
    )
    dual_pathways_df = compute_dual_robust_pathways(
        full_refs,
        pathway_robustness_df,
        drop_compounds,
        dual_retention_threshold,
    )
    dual_targets_df.to_csv(OUT_DIR / "dual_robust_targets.csv", index=False)
    dual_pathways_df.to_csv(OUT_DIR / "dual_robust_pathways.csv", index=False)

    influence_df = compute_compound_influence_scores(
        result_df,
        target_robustness_df,
        pathway_robustness_df,
        full_refs,
    )
    influence_df.to_csv(OUT_DIR / "compound_influence_scores.csv", index=False)
    plot_outputs(influence_df, dual_targets_df)

    n_expected = len(compounds_df) * len(PARAM_COMBINATIONS)
    n_complete = int(result_df["status"].isin(["ok", "skipped"]).sum()) if not result_df.empty else 0
    n_ok = int((result_df["status"] == "ok").sum()) if not result_df.empty else 0
    n_skipped = int((result_df["status"] == "skipped").sum()) if not result_df.empty else 0

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "analysis": "leave_one_compound_out_full_parameter_grid",
        "n_compounds": int(len(compounds_df)),
        "n_parameter_combinations": int(len(PARAM_COMBINATIONS)),
        "expected_combinations": int(n_expected),
        "completed_checkpoints": int(n_complete),
        "ok_combinations": int(n_ok),
        "skipped_combinations": int(n_skipped),
        "dual_retention_threshold": dual_retention_threshold,
        "dual_robust_target_count": int(dual_targets_df["is_dual_robust"].sum()) if not dual_targets_df.empty else 0,
        "dual_robust_pathway_count": int(dual_pathways_df["is_dual_robust"].sum()) if not dual_pathways_df.empty else 0,
        "top_influential_compounds": influence_df.head(10).to_dict(orient="records") if not influence_df.empty else [],
        "outputs": {
            "compound_sensitivity_results": str(OUT_DIR / "compound_sensitivity_results.csv"),
            "compound_influence_scores": str(OUT_DIR / "compound_influence_scores.csv"),
            "dual_robust_targets": str(OUT_DIR / "dual_robust_targets.csv"),
            "dual_robust_pathways": str(OUT_DIR / "dual_robust_pathways.csv"),
            "leave_one_out_target_robustness": str(OUT_DIR / "leave_one_out_target_robustness.csv"),
            "leave_one_out_pathway_robustness": str(OUT_DIR / "leave_one_out_pathway_robustness.csv"),
        },
    }
    atomic_write_json(OUT_DIR / "compound_sensitivity_summary.json", summary)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="rerun even completed checkpoints")
    parser.add_argument("--only-drop", action="append", default=None, help="run only selected compound stem(s)")
    parser.add_argument("--max-combos", type=int, default=None, help="debug limit for newly executed combos")
    parser.add_argument("--sleep-seconds", type=float, default=0.2, help="sleep between Enrichr library calls")
    parser.add_argument("--retries", type=int, default=3, help="Enrichr retry count per library")
    parser.add_argument("--retry-sleep-seconds", type=float, default=8.0, help="sleep before Enrichr retry")
    parser.add_argument("--dual-retention-threshold", type=float, default=0.80)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COMBO_CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_ENRICHMENT_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    raw_targets_df, compounds_df = load_compound_target_records()
    if args.only_drop:
        selected = set(args.only_drop)
        compounds_df = compounds_df[compounds_df["compound_stem"].isin(selected)].copy()
        missing = selected - set(compounds_df["compound_stem"])
        if missing:
            raise ValueError(f"Unknown --only-drop compound stem(s): {sorted(missing)}")

    disease_df = pd.read_csv(MERGED_DISEASE_GENES_CSV)
    disease_genes = set(disease_df["Gene_Symbol"].astype(str).str.upper())

    ppi_input = SENSITIVITY_PPI_EDGES_CSV if SENSITIVITY_PPI_EDGES_CSV.exists() else PPI_EDGES_CSV
    string_edges = pd.read_csv(ppi_input)
    string_graphs = build_string_graphs(string_edges)
    full_refs = load_full_references()

    print("=" * 78)
    print("Leave-one-compound-out full parameter-grid sensitivity analysis")
    print("=" * 78)
    print(f"Compounds: {len(compounds_df)}")
    print(f"Parameter combinations per compound: {len(PARAM_COMBINATIONS)}")
    print(f"Expected checkpoints: {len(compounds_df) * len(PARAM_COMBINATIONS)}")
    print(f"PPI input: {ppi_input}")
    print(f"Output dir: {OUT_DIR}")
    print(f"Checkpoint dir: {COMBO_CHECKPOINT_DIR}")
    print("")

    baseline_targets = deduplicate_targets(raw_targets_df, disease_genes, excluded_compound=None)
    print(f"Rebuilt full disease-filtered target universe: {len(baseline_targets)} targets")
    for _, compound in compounds_df.iterrows():
        drop_stem = compound["compound_stem"]
        loo_targets = deduplicate_targets(raw_targets_df, disease_genes, excluded_compound=drop_stem)
        compounds_df.loc[compounds_df["compound_stem"] == drop_stem, "loo_disease_filtered_targets"] = len(loo_targets)

    executed = 0
    expected_payloads: list[dict[str, Any]] = []

    for compound_index, compound in compounds_df.reset_index(drop=True).iterrows():
        drop_stem = str(compound["compound_stem"])
        drop_name = str(compound["compound_name"])
        pred_df = deduplicate_targets(raw_targets_df, disease_genes, excluded_compound=drop_stem)

        print(
            f"\n[{compound_index + 1}/{len(compounds_df)}] Drop {drop_stem} "
            f"-> {len(pred_df)} disease-filtered targets",
            flush=True,
        )

        for combo_index, (prob_th, string_th, top_k) in enumerate(PARAM_COMBINATIONS, start=1):
            combo = combo_id(prob_th, string_th, top_k)
            ckpt = checkpoint_path(drop_stem, combo)

            if ckpt.exists() and not args.force:
                payload = read_json(ckpt)
                if is_complete_checkpoint(payload):
                    payload["checkpoint_file"] = str(ckpt)
                    expected_payloads.append(payload)
                    continue

            if args.max_combos is not None and executed >= args.max_combos:
                continue

            filtered_targets = pred_df[pred_df["Probability"] > prob_th]["Target_Gene"].astype(str).tolist()
            start = time.time()
            try:
                result = run_single_combo(
                    filtered_targets=filtered_targets,
                    string_graph=string_graphs[string_th],
                    string_th=string_th,
                    pvalue_cutoff=PVALUE_CUTOFF,
                    top_k=top_k,
                    sleep_seconds=args.sleep_seconds,
                    retries=args.retries,
                    retry_sleep_seconds=args.retry_sleep_seconds,
                )
            except Exception as exc:
                result = {
                    "status": "error",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "n_targets": len(filtered_targets),
                }
            elapsed = time.time() - start
            payload = {
                **result,
                "drop_compound": drop_stem,
                "drop_compound_name": drop_name,
                "combo": combo,
                "prob_threshold": prob_th,
                "string_threshold": string_th,
                "pvalue_cutoff": PVALUE_CUTOFF,
                "top_k_param": top_k,
                "elapsed_seconds": round(elapsed, 3),
                "checkpoint_file": str(ckpt),
                "created_at": datetime.now().isoformat(timespec="seconds"),
            }
            atomic_write_json(ckpt, payload)
            expected_payloads.append(payload)
            executed += 1

            if result["status"] == "ok":
                print(
                    f"  [{combo_index:3d}/{len(PARAM_COMBINATIONS)}] OK   {combo} "
                    f"targets={result['n_targets']} ppi={result['n_ppi_nodes']}/{result['n_ppi_edges']} "
                    f"terms={result['significant_pathways']} cache={result.get('enrichment_cache_source')}",
                    flush=True,
                )
            elif result["status"] == "skipped":
                print(
                    f"  [{combo_index:3d}/{len(PARAM_COMBINATIONS)}] SKIP {combo} "
                    f"reason={result.get('skip_reason')}",
                    flush=True,
                )
            else:
                print(
                    f"  [{combo_index:3d}/{len(PARAM_COMBINATIONS)}] ERROR {combo} "
                    f"{result.get('error_type')}: {result.get('error_message')}",
                    flush=True,
                )

    # Load all complete checkpoints for selected compounds so aggregation is complete
    # even when this run only executed a subset of missing combos.
    all_payloads: list[dict[str, Any]] = []
    for _, compound in compounds_df.iterrows():
        drop_stem = str(compound["compound_stem"])
        for prob_th, string_th, top_k in PARAM_COMBINATIONS:
            combo = combo_id(prob_th, string_th, top_k)
            ckpt = checkpoint_path(drop_stem, combo)
            if ckpt.exists():
                payload = read_json(ckpt)
                if is_complete_checkpoint(payload):
                    payload["checkpoint_file"] = str(ckpt)
                    all_payloads.append(payload)

    summary = aggregate_outputs(
        expected_payloads=all_payloads,
        compounds_df=compounds_df,
        full_refs=full_refs,
        dual_retention_threshold=args.dual_retention_threshold,
    )

    print("\n" + "=" * 78)
    print("Compound sensitivity aggregation complete")
    print("=" * 78)
    print(f"Completed checkpoints: {summary['completed_checkpoints']}/{summary['expected_combinations']}")
    print(f"OK combinations: {summary['ok_combinations']}")
    print(f"Skipped combinations: {summary['skipped_combinations']}")
    print(f"Dual-robust targets: {summary['dual_robust_target_count']}")
    print(f"Dual-robust pathways: {summary['dual_robust_pathway_count']}")
    print(f"Summary: {OUT_DIR / 'compound_sensitivity_summary.json'}")


if __name__ == "__main__":
    main()
