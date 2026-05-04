"""Build the robustness validation bundle for manuscript support.

The bundle is intentionally reproducible offline.  It adds three pieces of
evidence on top of the locked robustness outputs:

1. Source-support validation with hypergeometric tests and BH-FDR correction.
2. Degree-matched permutation tests for target-source support metrics.
3. Separate ablation tables for parameter-robust Top30 target selection and
   leave-one-compound-out antioxidant-network perturbation.

PubMed target co-occurrence is optional because it requires network access.

Usage
-----
cd code
python scripts/analysis/robustness_validation_bundle.py --n-permutations 1000
python scripts/analysis/robustness_validation_bundle.py --pubmed-mode query
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.connectors.pubmed_connector import PubMedConfig, PubMedConnector  # noqa: E402
from src.pipeline.stage4_merge_disease_genes import (  # noqa: E402
    load_genecards,
    load_omim,
    load_ttd,
)


OUTPUTS_DIR = ROOT / "src" / "outputs"
SHARED_DIR = OUTPUTS_DIR / "shared"
SENSITIVITY_DIR = OUTPUTS_DIR / "sensitivity"
COMPOUND_SENSITIVITY_DIR = SENSITIVITY_DIR / "compound_sensitivity"
SUPPLEMENTARY_DIR = OUTPUTS_DIR / "supplementary"
DEFAULT_OUTPUT_DIR = SENSITIVITY_DIR / "robustness_validation"

TARGET_ROBUSTNESS = SENSITIVITY_DIR / "target_robustness.csv"
DUAL_ROBUST_TARGETS = COMPOUND_SENSITIVITY_DIR / "dual_robust_targets.csv"
POSITIVE_TARGET_UNIVERSE = SENSITIVITY_DIR / "intersected_targets_probability_gt0.csv"
DISEASE_GENE_DIR = SHARED_DIR / "disease_genes"
FULL_NODE_DEGREES = SENSITIVITY_DIR / "ppi_network_full" / "node_degrees.csv"
MAIN_HUB_NODES = SHARED_DIR / "ppi_network" / "hub_nodes.csv"
PARAM_GRID_RESULTS = SENSITIVITY_DIR / "param_grid_results.csv"
COMBO_TARGETS = SENSITIVITY_DIR / "combo_targets.json"
COMPOUND_INFLUENCE = COMPOUND_SENSITIVITY_DIR / "compound_influence_scores.csv"

OKABE_ITO = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "gray": "#7F7F7F",
}

PUBMED_KEY_THRESHOLDS = (10, 20, 50, 100, 200)


@dataclass(frozen=True)
class TargetSet:
    """Named target set used in validation and ablation tables."""

    set_name: str
    display_name: str
    genes: tuple[str, ...]
    set_type: str


def normalize_gene(value: Any) -> str:
    """Normalize a gene symbol for set comparisons."""
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().upper()


def boolean_true(series: pd.Series) -> pd.Series:
    """Parse common truthy values from CSV columns."""
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def save_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)


def frame_to_markdown(df: pd.DataFrame) -> str:
    """Render a compact markdown table without the optional tabulate dependency."""
    if df.empty:
        return "_No rows._"
    display_df = df.copy()
    display_df = display_df.where(pd.notna(display_df), "NA")
    columns = [str(column) for column in display_df.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in display_df.iterrows():
        values = [str(row[column]).replace("\n", " ") for column in display_df.columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def bh_adjust(p_values: Sequence[float | None]) -> list[float | None]:
    """Benjamini-Hochberg FDR adjustment preserving missing positions."""
    numeric: list[tuple[int, float]] = []
    for idx, value in enumerate(p_values):
        if value is None:
            continue
        if isinstance(value, float) and math.isnan(value):
            continue
        numeric.append((idx, float(value)))

    adjusted: list[float | None] = [None] * len(p_values)
    if not numeric:
        return adjusted

    ordered = sorted(numeric, key=lambda item: item[1])
    m = len(ordered)
    running_min = 1.0
    for rank_from_end, (idx, p_value) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        running_min = min(running_min, p_value * m / rank)
        adjusted[idx] = min(running_min, 1.0)
    return adjusted


def hypergeom_support_test(
    target_genes: Iterable[str],
    validation_genes: Iterable[str],
    background_genes: Iterable[str],
) -> dict[str, Any]:
    """One-sided enrichment test for overlap with a validation gene set."""
    background = {normalize_gene(gene) for gene in background_genes}
    background.discard("")
    target_set = {normalize_gene(gene) for gene in target_genes} & background
    validation_set = {normalize_gene(gene) for gene in validation_genes} & background

    m = len(background)
    n_success = len(validation_set)
    n_draws = len(target_set)
    overlap_genes = sorted(target_set & validation_set)
    k_overlap = len(overlap_genes)

    if m == 0 or n_draws == 0 or n_success == 0:
        p_value = 1.0
    else:
        p_value = float(hypergeom.sf(k_overlap - 1, m, n_success, n_draws))

    expected = n_draws * n_success / m if m else 0.0
    target_rate = k_overlap / n_draws if n_draws else 0.0
    background_rate = n_success / m if m else 0.0
    enrichment_ratio = target_rate / background_rate if background_rate else math.inf

    return {
        "background_n": m,
        "validation_n_in_background": n_success,
        "target_n_in_background": n_draws,
        "overlap_n": k_overlap,
        "expected_overlap": expected,
        "observed_minus_expected": k_overlap - expected,
        "target_support_rate": target_rate,
        "background_support_rate": background_rate,
        "enrichment_ratio": enrichment_ratio,
        "p_value": p_value,
        "overlap_genes": "; ".join(overlap_genes),
    }


def empirical_pvalue_greater_equal(observed: float, null_values: Sequence[float]) -> float:
    if not null_values:
        return math.nan
    return (sum(value >= observed for value in null_values) + 1) / (len(null_values) + 1)


def percentile_less(observed: float, null_values: Sequence[float]) -> float:
    if not null_values:
        return math.nan
    return 100.0 * sum(value < observed for value in null_values) / len(null_values)


def load_target_universe() -> set[str]:
    df = safe_read_csv(POSITIVE_TARGET_UNIVERSE)
    gene_col = "Target_Gene" if "Target_Gene" in df.columns else df.columns[0]
    genes = {normalize_gene(gene) for gene in df[gene_col]}
    genes.discard("")
    return genes


def load_candidate_target_background() -> set[str]:
    """Alias with manuscript-facing terminology for the target background set."""
    return load_target_universe()


def load_source_sets() -> dict[str, set[str]]:
    """Load source-specific support sets used for target validation."""
    loaders = {
        "GeneCards": (load_genecards, DISEASE_GENE_DIR / "GeneCards-SearchResults.csv"),
        "OMIM": (load_omim, DISEASE_GENE_DIR / "OMIM-Entry-Search.csv"),
        "TTD": (load_ttd, DISEASE_GENE_DIR / "TTD_Search_Results.xlsx"),
    }
    source_sets: dict[str, set[str]] = {}
    for source, (loader, path) in loaders.items():
        if not path.exists():
            source_sets[source] = set()
            continue
        df = loader(str(path))
        genes = {normalize_gene(gene) for gene in df["Gene_Symbol"]}
        genes.discard("")
        source_sets[source] = genes
    any_source = set().union(*source_sets.values()) if source_sets else set()
    source_sets["AnySource"] = any_source
    return source_sets


def load_gene_descriptions() -> dict[str, str]:
    path = DISEASE_GENE_DIR / "merged_disease_genes.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    descriptions: dict[str, str] = {}
    for _, row in df.iterrows():
        gene = normalize_gene(row.get("Gene_Symbol"))
        desc = str(row.get("Description", "")).strip()
        if gene and desc and desc.lower() != "nan":
            descriptions[gene] = desc
    return descriptions


def load_target_sets() -> list[TargetSet]:
    target_df = safe_read_csv(TARGET_ROBUSTNESS)
    target_df["gene"] = target_df["gene"].map(normalize_gene)
    q1 = target_df[target_df["quadrant"].astype(str).str.contains("Q1", na=False)].copy()
    q1 = q1.sort_values(
        ["frequency", "avg_degree_norm", "avg_degree_rank", "avg_betweenness"],
        ascending=[False, False, True, False],
    )

    dual_df = safe_read_csv(DUAL_ROBUST_TARGETS)
    dual_df["gene"] = dual_df["gene"].map(normalize_gene)
    dual = dual_df[boolean_true(dual_df["is_dual_robust"])].copy()
    tier1 = dual[
        (dual["full_frequency"] == 1.0) & (dual["loo_q1_retention_rate"] == 1.0)
    ].copy()

    hub_df = safe_read_csv(MAIN_HUB_NODES)
    hub_df["Gene"] = hub_df["Gene"].map(normalize_gene)
    traditional_top30 = tuple(hub_df.sort_values("Rank").head(30)["Gene"].tolist())

    return [
        TargetSet(
            set_name="traditional_top30",
            display_name="Traditional fixed Top30 hubs",
            genes=traditional_top30,
            set_type="baseline",
        ),
        TargetSet(
            set_name="q1_candidate_robust",
            display_name="Parameter-Q1 robust targets",
            genes=tuple(q1["gene"].tolist()),
            set_type="parameter_robust",
        ),
        TargetSet(
            set_name="dual_robust",
            display_name="Dual-robust targets",
            genes=tuple(dual["gene"].tolist()),
            set_type="dual_robust",
        ),
        TargetSet(
            set_name="tier1",
            display_name="Tier 1 target",
            genes=tuple(tier1["gene"].tolist()),
            set_type="tier1",
        ),
    ]


def build_gene_pubmed_query(gene: str, description: str | None, context_term: str) -> str:
    """Build a conservative Title/Abstract co-occurrence query for one gene."""
    symbol = normalize_gene(gene)
    safe_context = context_term.replace('"', "")
    terms = []
    if len(symbol) > 3:
        terms.append(f'"{symbol}"[Title/Abstract]')
    terms.extend(
        [
            f'"{symbol} gene"[Title/Abstract]',
            f'"{symbol} protein"[Title/Abstract]',
        ]
    )
    if description:
        cleaned = description.replace('"', "").strip()
        if 4 <= len(cleaned) <= 90 and ";" not in cleaned:
            terms.append(f'"{cleaned}"[Title/Abstract]')
    return f"({' OR '.join(dict.fromkeys(terms))}) AND (\"{safe_context}\"[Title/Abstract])"


def read_pubmed_cache(cache_dir: Path) -> dict[str, int]:
    cache_path = cache_dir / "pubmed_cache.json"
    if not cache_path.exists():
        return {}
    with open(cache_path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    counts: dict[str, int] = {}
    for key, value in payload.items():
        try:
            counts[key] = int(value["count"])
        except (KeyError, TypeError, ValueError):
            continue
    return counts


def collect_pubmed_counts(
    genes: Sequence[str],
    descriptions: dict[str, str],
    output_dir: Path,
    mode: str,
    context_term: str,
    delay_seconds: float,
) -> pd.DataFrame:
    """Collect target-level PubMed co-occurrence counts.

    mode="cache" never touches the network; mode="query" uses NCBI E-utilities.
    """
    rows: list[dict[str, Any]] = []
    unique_genes = sorted({normalize_gene(gene) for gene in genes if normalize_gene(gene)})
    output_dir.mkdir(parents=True, exist_ok=True)

    cached_counts = read_pubmed_cache(output_dir)
    connector: PubMedConnector | None = None
    if mode == "query":
        connector = PubMedConnector(
            config=PubMedConfig(delay_seconds=delay_seconds),
            cache_dir=output_dir,
        )

    for idx, gene in enumerate(unique_genes, start=1):
        query = build_gene_pubmed_query(gene, descriptions.get(gene), context_term)
        count: int | None = None
        status = "not_queried"

        if mode == "skip":
            status = "skipped"
        elif mode == "cache":
            if query in cached_counts:
                count = cached_counts[query]
                status = "cache_hit"
            else:
                status = "not_cached"
        elif connector is not None:
            count = connector.search_count(query, use_cache=True)
            status = "queried" if count >= 0 else "query_failed"
            print(f"  PubMed [{idx:3d}/{len(unique_genes)}] {gene:<10} n={count}")
        else:
            raise ValueError(f"Unsupported pubmed mode: {mode}")

        rows.append(
            {
                "gene": gene,
                "pubmed_count": count,
                "query_status": status,
                "context_term": context_term,
                "query": query,
            }
        )

    return pd.DataFrame(rows)


def run_external_validation(
    target_sets: Sequence[TargetSet],
    source_sets: dict[str, set[str]],
    background: set[str],
    pubmed_df: pd.DataFrame,
    pubmed_threshold: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    validation_sets: dict[str, set[str]] = dict(source_sets)

    if not pubmed_df.empty and "pubmed_count" in pubmed_df.columns:
        pubmed_valid = pubmed_df.dropna(subset=["pubmed_count"]).copy()
        pubmed_valid["pubmed_count"] = pd.to_numeric(pubmed_valid["pubmed_count"], errors="coerce")
        pubmed_valid = pubmed_valid[pubmed_valid["pubmed_count"] >= 0]
        if not pubmed_valid.empty:
            pubmed_genes = set(
                pubmed_valid.loc[
                    pubmed_valid["pubmed_count"] >= pubmed_threshold,
                    "gene",
                ].map(normalize_gene)
            )
            validation_sets[f"PubMed_ge_{pubmed_threshold}"] = pubmed_genes

    for target_set in target_sets:
        for source_name, validation_genes in validation_sets.items():
            result = hypergeom_support_test(target_set.genes, validation_genes, background)
            rows.append(
                {
                    "target_set": target_set.set_name,
                    "display_name": target_set.display_name,
                    "set_type": target_set.set_type,
                    "validation_source": source_name,
                    **result,
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df["fdr_bh"] = bh_adjust(df["p_value"].tolist())
        df["minus_log10_fdr"] = df["fdr_bh"].apply(
            lambda value: -math.log10(max(float(value), 1e-300)) if pd.notna(value) else math.nan
        )
    return df


def build_pubmed_threshold_sensitivity(
    target_sets: Sequence[TargetSet],
    pubmed_df: pd.DataFrame,
    background: set[str],
    thresholds: Sequence[int],
) -> pd.DataFrame:
    """Test PubMed support enrichment over multiple count thresholds."""
    if pubmed_df.empty or "pubmed_count" not in pubmed_df.columns:
        return pd.DataFrame()

    counts_df = pubmed_df.copy()
    counts_df["gene"] = counts_df["gene"].map(normalize_gene)
    counts_df["pubmed_count"] = pd.to_numeric(counts_df["pubmed_count"], errors="coerce")
    counts_df = counts_df[counts_df["pubmed_count"].notna() & (counts_df["pubmed_count"] >= 0)]
    if counts_df.empty:
        return pd.DataFrame()

    count_map = dict(zip(counts_df["gene"], counts_df["pubmed_count"].astype(int)))
    rows: list[dict[str, Any]] = []
    for threshold in thresholds:
        validation_genes = {gene for gene, count in count_map.items() if count >= threshold}
        for target_set in target_sets:
            result = hypergeom_support_test(target_set.genes, validation_genes, background)
            rows.append(
                {
                    "target_set": target_set.set_name,
                    "display_name": target_set.display_name,
                    "pubmed_threshold": threshold,
                    **result,
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df["fdr_bh"] = bh_adjust(df["p_value"].tolist())
    return df


def build_pubmed_key_results(
    pubmed_threshold_df: pd.DataFrame,
    thresholds: Sequence[int] = PUBMED_KEY_THRESHOLDS,
) -> pd.DataFrame:
    """Return the manuscript-facing PubMed co-occurrence threshold table."""
    if pubmed_threshold_df.empty:
        return pd.DataFrame()
    target_order = ["traditional_top30", "q1_candidate_robust", "dual_robust"]
    key_df = pubmed_threshold_df[
        pubmed_threshold_df["target_set"].isin(target_order)
        & pubmed_threshold_df["pubmed_threshold"].isin(thresholds)
    ].copy()
    if key_df.empty:
        return pd.DataFrame()
    key_df["target_set"] = pd.Categorical(key_df["target_set"], categories=target_order, ordered=True)
    key_df["pubmed_threshold"] = pd.Categorical(
        key_df["pubmed_threshold"],
        categories=list(thresholds),
        ordered=True,
    )
    key_df = key_df.sort_values(["pubmed_threshold", "target_set"]).reset_index(drop=True)
    key_df["pubmed_threshold"] = key_df["pubmed_threshold"].astype(int)
    return key_df[
        [
            "target_set",
            "display_name",
            "pubmed_threshold",
            "validation_n_in_background",
            "target_n_in_background",
            "overlap_n",
            "expected_overlap",
            "target_support_rate",
            "background_support_rate",
            "enrichment_ratio",
            "p_value",
            "fdr_bh",
            "overlap_genes",
        ]
    ]


def build_dual_pubmed_counts(
    pubmed_df: pd.DataFrame,
    target_sets: Sequence[TargetSet],
) -> pd.DataFrame:
    """List dual-robust targets with their PubMed co-occurrence counts."""
    if pubmed_df.empty or "pubmed_count" not in pubmed_df.columns:
        return pd.DataFrame()
    dual_sets = [target_set for target_set in target_sets if target_set.set_name == "dual_robust"]
    if not dual_sets:
        return pd.DataFrame()
    dual_genes = {normalize_gene(gene) for gene in dual_sets[0].genes}
    counts_df = pubmed_df.copy()
    counts_df["gene"] = counts_df["gene"].map(normalize_gene)
    counts_df["pubmed_count"] = pd.to_numeric(counts_df["pubmed_count"], errors="coerce")
    counts_df = counts_df[counts_df["gene"].isin(dual_genes)].copy()
    if counts_df.empty:
        return pd.DataFrame()
    return counts_df.sort_values(["pubmed_count", "gene"], ascending=[False, True]).reset_index(drop=True)


def load_pubmed_threshold_degree_summary(output_dir: Path) -> pd.DataFrame:
    """Load optional multi-threshold degree-matched PubMed results if present."""
    path = output_dir / "pubmed_threshold_degree_matched_summary.csv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    if df.empty or "metric" not in df.columns:
        return pd.DataFrame()
    focus = df[df["metric"].astype(str).str.match(r"PubMed_ge_\d+_count")].copy()
    if focus.empty:
        return pd.DataFrame()
    return focus.sort_values(["pubmed_threshold", "target_set"]).reset_index(drop=True)


def load_degree_table(background: set[str]) -> pd.DataFrame:
    df = safe_read_csv(FULL_NODE_DEGREES)
    df = df.rename(columns={"Gene": "gene", "Degree": "degree"})
    df["gene"] = df["gene"].map(normalize_gene)
    df["degree"] = pd.to_numeric(df["degree"], errors="coerce").fillna(0.0).astype(float)
    df = df[df["gene"].isin(background)].copy()
    return df[["gene", "degree"]].drop_duplicates("gene").reset_index(drop=True)


def weighted_degree_choice(pool: pd.DataFrame, target_degree: float, rng: random.Random) -> str:
    distances = (pool["degree"] - target_degree).abs().to_numpy(dtype=float)
    weights = 1.0 / np.power(distances + 1.0, 1.25)
    pick = rng.random() * float(weights.sum())
    cumulative = 0.0
    for gene, weight in zip(pool["gene"].tolist(), weights):
        cumulative += float(weight)
        if cumulative >= pick:
            return str(gene)
    return str(pool["gene"].iloc[-1])


def sample_degree_matched_genes(
    observed_genes: Sequence[str],
    degree_table: pd.DataFrame,
    rng: random.Random,
    match_pool_size: int,
) -> tuple[list[str], list[dict[str, Any]]]:
    degree_map = dict(zip(degree_table["gene"], degree_table["degree"]))
    target_rows = [
        {"gene": gene, "degree": float(degree_map[gene])}
        for gene in observed_genes
        if gene in degree_map
    ]
    if not target_rows:
        raise ValueError("No observed genes have degree values.")

    selected: set[str] = set()
    excluded = {row["gene"] for row in target_rows}
    available = degree_table[~degree_table["gene"].isin(excluded)].copy()
    sample_rows: list[dict[str, Any]] = []

    shuffled = target_rows[:]
    rng.shuffle(shuffled)
    for target in shuffled:
        pool = available[~available["gene"].isin(selected)].copy()
        pool["degree_distance"] = (pool["degree"] - target["degree"]).abs()
        pool = pool.nsmallest(min(match_pool_size, len(pool)), "degree_distance")
        if pool.empty:
            raise RuntimeError("No candidates left for degree-matched sampling.")

        picked = weighted_degree_choice(pool, float(target["degree"]), rng)
        picked_degree = float(pool.loc[pool["gene"] == picked, "degree"].iloc[0])
        selected.add(picked)
        sample_rows.append(
            {
                "matched_target": target["gene"],
                "matched_target_degree": float(target["degree"]),
                "gene": picked,
                "degree": picked_degree,
                "degree_distance": abs(picked_degree - float(target["degree"])),
            }
        )

    return [row["gene"] for row in sample_rows], sample_rows


def support_metrics_for_genes(
    genes: Sequence[str],
    source_sets: dict[str, set[str]],
    pubmed_counts: dict[str, int] | None,
    pubmed_threshold: int,
    degree_map: dict[str, float],
) -> dict[str, Any]:
    gene_set = {normalize_gene(gene) for gene in genes}
    gene_set.discard("")
    degrees = [float(degree_map[gene]) for gene in gene_set if gene in degree_map]

    row: dict[str, Any] = {
        "n_targets": len(gene_set),
        "mean_degree": float(np.mean(degrees)) if degrees else math.nan,
        "median_degree": float(np.median(degrees)) if degrees else math.nan,
    }
    for source_name, source_genes in source_sets.items():
        overlap = gene_set & source_genes
        row[f"{source_name}_count"] = len(overlap)
        row[f"{source_name}_rate"] = len(overlap) / len(gene_set) if gene_set else 0.0

    if pubmed_counts:
        known = {gene: count for gene, count in pubmed_counts.items() if gene in gene_set and count >= 0}
        positives = {gene for gene, count in known.items() if count >= pubmed_threshold}
        row["PubMed_known_count"] = len(known)
        row[f"PubMed_ge_{pubmed_threshold}_count"] = len(positives)
        row[f"PubMed_ge_{pubmed_threshold}_rate"] = len(positives) / len(known) if known else math.nan
        row["PubMed_count_sum"] = int(sum(known.values())) if known else 0
        row["PubMed_count_median"] = float(np.median(list(known.values()))) if known else math.nan
    return row


def run_degree_matched_permutation(
    target_sets: Sequence[TargetSet],
    source_sets: dict[str, set[str]],
    degree_table: pd.DataFrame,
    n_permutations: int,
    seed: int,
    match_pool_size: int,
    pubmed_df: pd.DataFrame,
    pubmed_threshold: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = random.Random(seed)
    degree_map = dict(zip(degree_table["gene"], degree_table["degree"]))
    pubmed_counts: dict[str, int] = {}
    if not pubmed_df.empty and "pubmed_count" in pubmed_df.columns:
        for _, row in pubmed_df.dropna(subset=["pubmed_count"]).iterrows():
            pubmed_counts[normalize_gene(row["gene"])] = int(row["pubmed_count"])

    metric_rows: list[dict[str, Any]] = []
    random_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    permutation_sets = [item for item in target_sets if item.set_name != "tier1"]
    for target_set in permutation_sets:
        observed_genes = [gene for gene in target_set.genes if gene in degree_map]
        observed_metrics = support_metrics_for_genes(
            observed_genes,
            source_sets,
            pubmed_counts or None,
            pubmed_threshold,
            degree_map,
        )
        metric_rows.append(
            {
                "target_set": target_set.set_name,
                "sample_id": "observed",
                "is_observed": True,
                **observed_metrics,
            }
        )

        for i in range(1, n_permutations + 1):
            sampled_genes, sample_rows = sample_degree_matched_genes(
                observed_genes,
                degree_table,
                rng,
                match_pool_size,
            )
            sample_id = f"{target_set.set_name}_random_{i:04d}"
            for row in sample_rows:
                row["target_set"] = target_set.set_name
                row["sample_id"] = sample_id
            random_rows.extend(sample_rows)
            metric_rows.append(
                {
                    "target_set": target_set.set_name,
                    "sample_id": sample_id,
                    "is_observed": False,
                    **support_metrics_for_genes(
                        sampled_genes,
                        source_sets,
                        pubmed_counts or None,
                        pubmed_threshold,
                        degree_map,
                    ),
                }
            )

    metrics_df = pd.DataFrame(metric_rows)
    random_sets_df = pd.DataFrame(random_rows)

    skip_columns = {"target_set", "sample_id", "is_observed"}
    metric_columns = [
        column
        for column in metrics_df.columns
        if column not in skip_columns and pd.api.types.is_numeric_dtype(metrics_df[column])
    ]
    for target_set in sorted(metrics_df["target_set"].unique()):
        observed = metrics_df[(metrics_df["target_set"] == target_set) & metrics_df["is_observed"]].iloc[0]
        null = metrics_df[(metrics_df["target_set"] == target_set) & ~metrics_df["is_observed"]]
        for metric in metric_columns:
            null_values = [float(value) for value in null[metric].dropna().tolist()]
            observed_value = float(observed[metric])
            if not null_values:
                continue
            null_mean = float(np.mean(null_values))
            null_sd = float(np.std(null_values, ddof=1)) if len(null_values) > 1 else 0.0
            summary_rows.append(
                {
                    "target_set": target_set,
                    "metric": metric,
                    "observed": observed_value,
                    "null_mean": null_mean,
                    "null_sd": null_sd,
                    "null_q025": float(np.quantile(null_values, 0.025)),
                    "null_q975": float(np.quantile(null_values, 0.975)),
                    "observed_percentile": percentile_less(observed_value, null_values),
                    "empirical_p_greater_equal": empirical_pvalue_greater_equal(observed_value, null_values),
                    "z_score": (observed_value - null_mean) / null_sd if null_sd > 0 else math.nan,
                    "n_permutations": n_permutations,
                    "seed": seed,
                    "match_pool_size": match_pool_size,
                }
            )

    return metrics_df, random_sets_df, pd.DataFrame(summary_rows)


def build_ablation_table(
    target_sets: Sequence[TargetSet],
    source_sets: dict[str, set[str]],
    degree_table: pd.DataFrame,
    pubmed_df: pd.DataFrame,
    pubmed_threshold: int,
) -> pd.DataFrame:
    target_robustness = safe_read_csv(TARGET_ROBUSTNESS)
    target_robustness["gene"] = target_robustness["gene"].map(normalize_gene)
    robust_by_gene = target_robustness.drop_duplicates("gene").set_index("gene")

    dual_df = safe_read_csv(DUAL_ROBUST_TARGETS)
    dual_df["gene"] = dual_df["gene"].map(normalize_gene)
    dual_by_gene = dual_df.drop_duplicates("gene").set_index("gene")
    dual_genes = set(dual_df.loc[boolean_true(dual_df["is_dual_robust"]), "gene"])
    q1_genes = set(target_sets[1].genes)
    degree_map = dict(zip(degree_table["gene"], degree_table["degree"]))

    pubmed_counts: dict[str, int] = {}
    if not pubmed_df.empty and "pubmed_count" in pubmed_df.columns:
        for _, row in pubmed_df.dropna(subset=["pubmed_count"]).iterrows():
            pubmed_counts[normalize_gene(row["gene"])] = int(row["pubmed_count"])

    rows: list[dict[str, Any]] = []
    for target_set in target_sets:
        genes = {normalize_gene(gene) for gene in target_set.genes}
        genes.discard("")
        robust_subset = robust_by_gene.loc[robust_by_gene.index.intersection(genes)]
        dual_subset = dual_by_gene.loc[dual_by_gene.index.intersection(genes)]
        support = support_metrics_for_genes(
            sorted(genes),
            source_sets,
            pubmed_counts or None,
            pubmed_threshold,
            degree_map,
        )

        full_frequency = pd.to_numeric(robust_subset.get("frequency"), errors="coerce")
        avg_degree_norm = pd.to_numeric(robust_subset.get("avg_degree_norm"), errors="coerce")
        loo_retention = pd.to_numeric(
            dual_subset.get("loo_q1_retention_rate"),
            errors="coerce",
        )
        row: dict[str, Any] = {
            "target_set": target_set.set_name,
            "display_name": target_set.display_name,
            "set_type": target_set.set_type,
            "n_targets": len(genes),
            "q1_overlap_n": len(genes & q1_genes),
            "q1_overlap_rate": len(genes & q1_genes) / len(genes) if genes else 0.0,
            "dual_overlap_n": len(genes & dual_genes),
            "dual_overlap_rate": len(genes & dual_genes) / len(genes) if genes else 0.0,
            "mean_parameter_frequency": float(full_frequency.mean()) if len(full_frequency) else math.nan,
            "median_parameter_frequency": float(full_frequency.median()) if len(full_frequency) else math.nan,
            "mean_degree_norm": float(avg_degree_norm.mean()) if len(avg_degree_norm) else math.nan,
            "mean_loo_q1_retention_rate": float(loo_retention.mean()) if len(loo_retention) else math.nan,
            **support,
            "genes": "; ".join(sorted(genes)),
        }
        rows.append(row)

    return pd.DataFrame(rows)


def load_combo_targets() -> dict[str, list[str]]:
    if not COMBO_TARGETS.exists():
        raise FileNotFoundError(f"Missing required file: {COMBO_TARGETS}")
    with open(COMBO_TARGETS, "r", encoding="utf-8") as fh:
        raw = json.load(fh)
    return {
        str(combo): [normalize_gene(gene) for gene in genes if normalize_gene(gene)]
        for combo, genes in raw.items()
    }


def select_best_top30_combo(param_grid: pd.DataFrame) -> pd.Series:
    top30 = param_grid[pd.to_numeric(param_grid["top_k"], errors="coerce") == 30].copy()
    if top30.empty:
        raise ValueError("No TopK=30 rows are available in the parameter grid.")

    for column in [
        "direct_ox_count",
        "cross_talk_count",
        "backbone_count",
        "significant_pathways",
        "avg_ef",
        "min_pvalue",
    ]:
        top30[column] = pd.to_numeric(top30[column], errors="coerce")
    top30["oxidative_relevant_terms"] = (
        top30["direct_ox_count"].fillna(0)
        + top30["cross_talk_count"].fillna(0)
        + top30["backbone_count"].fillna(0)
    )
    top30 = top30.sort_values(
        [
            "oxidative_relevant_terms",
            "direct_ox_count",
            "significant_pathways",
            "avg_ef",
            "min_pvalue",
        ],
        ascending=[False, False, False, False, True],
    )
    return top30.iloc[0]


def summarize_parameter_target_set(
    row_id: str,
    display_name: str,
    selection_rule: str,
    genes: Sequence[str],
    robust_by_gene: pd.DataFrame,
    q1_top30_genes: set[str],
    q1_all_genes: set[str],
    combo_row: pd.Series | None = None,
) -> dict[str, Any]:
    normalized_genes = [normalize_gene(gene) for gene in genes if normalize_gene(gene)]
    normalized_genes = list(dict.fromkeys(normalized_genes))[:30]
    subset = robust_by_gene.loc[robust_by_gene.index.intersection(normalized_genes)]

    frequencies = pd.to_numeric(subset.get("frequency"), errors="coerce")
    degree_norm = pd.to_numeric(subset.get("avg_degree_norm"), errors="coerce")
    degree_rank = pd.to_numeric(subset.get("avg_degree_rank"), errors="coerce")
    betweenness = pd.to_numeric(subset.get("avg_betweenness"), errors="coerce")
    observed_gene_set = set(normalized_genes)

    output: dict[str, Any] = {
        "method": row_id,
        "display_name": display_name,
        "selection_rule": selection_rule,
        "n_targets": len(normalized_genes),
        "q1_top30_overlap_n": len(observed_gene_set & q1_top30_genes),
        "q1_top30_overlap_rate": len(observed_gene_set & q1_top30_genes) / len(normalized_genes)
        if normalized_genes
        else 0.0,
        "q1_all_overlap_n": len(observed_gene_set & q1_all_genes),
        "q1_all_overlap_rate": len(observed_gene_set & q1_all_genes) / len(normalized_genes)
        if normalized_genes
        else 0.0,
        "mean_parameter_frequency": float(frequencies.mean()) if len(frequencies) else math.nan,
        "median_parameter_frequency": float(frequencies.median()) if len(frequencies) else math.nan,
        "min_parameter_frequency": float(frequencies.min()) if len(frequencies) else math.nan,
        "mean_degree_norm": float(degree_norm.mean()) if len(degree_norm) else math.nan,
        "median_avg_degree_rank": float(degree_rank.median()) if len(degree_rank) else math.nan,
        "mean_betweenness": float(betweenness.mean()) if len(betweenness) else math.nan,
        "genes": "; ".join(normalized_genes),
    }
    if combo_row is not None:
        output.update(
            {
                "representative_combo": str(combo_row.get("combo", "")),
                "prob_threshold": float(combo_row.get("prob_threshold"))
                if pd.notna(combo_row.get("prob_threshold"))
                else math.nan,
                "string_threshold": float(combo_row.get("string_threshold"))
                if pd.notna(combo_row.get("string_threshold"))
                else math.nan,
                "top_k": int(combo_row.get("top_k")) if pd.notna(combo_row.get("top_k")) else 30,
                "significant_pathways": int(combo_row.get("significant_pathways"))
                if pd.notna(combo_row.get("significant_pathways"))
                else math.nan,
                "oxidative_relevant_terms": int(combo_row.get("oxidative_relevant_terms"))
                if pd.notna(combo_row.get("oxidative_relevant_terms"))
                else math.nan,
                "direct_ox_count": int(combo_row.get("direct_ox_count"))
                if pd.notna(combo_row.get("direct_ox_count"))
                else math.nan,
                "avg_ef": float(combo_row.get("avg_ef")) if pd.notna(combo_row.get("avg_ef")) else math.nan,
                "min_pvalue": float(combo_row.get("min_pvalue"))
                if pd.notna(combo_row.get("min_pvalue"))
                else math.nan,
            }
        )
    else:
        output.update(
            {
                "representative_combo": "across_grid",
                "prob_threshold": math.nan,
                "string_threshold": math.nan,
                "top_k": 30,
                "significant_pathways": math.nan,
                "oxidative_relevant_terms": math.nan,
                "direct_ox_count": math.nan,
                "avg_ef": math.nan,
                "min_pvalue": math.nan,
            }
        )
    return output


def build_parameter_top30_ablation() -> pd.DataFrame:
    target_robustness = safe_read_csv(TARGET_ROBUSTNESS)
    target_robustness["gene"] = target_robustness["gene"].map(normalize_gene)
    target_robustness = target_robustness.drop_duplicates("gene")
    robust_by_gene = target_robustness.set_index("gene")

    q1 = target_robustness[
        target_robustness["quadrant"].astype(str).str.contains("Q1", na=False)
    ].copy()
    q1_sorted = q1.sort_values(
        ["frequency", "avg_degree_norm", "avg_degree_rank", "avg_betweenness"],
        ascending=[False, False, True, False],
    )
    q1_top30 = q1_sorted.head(30)["gene"].tolist()
    q1_top30_genes = set(q1_top30)
    q1_all_genes = set(q1["gene"])

    frequency_top30 = target_robustness.sort_values(
        ["frequency", "avg_degree_norm", "avg_degree_rank", "avg_betweenness"],
        ascending=[False, False, True, False],
    ).head(30)["gene"].tolist()

    hub_df = safe_read_csv(MAIN_HUB_NODES)
    hub_df["Gene"] = hub_df["Gene"].map(normalize_gene)
    default_top30 = hub_df.sort_values("Rank").head(30)["Gene"].tolist()

    param_grid = safe_read_csv(PARAM_GRID_RESULTS)
    for column in ["prob_threshold", "string_threshold", "top_k"]:
        param_grid[column] = pd.to_numeric(param_grid[column], errors="coerce")
    param_grid["oxidative_relevant_terms"] = (
        pd.to_numeric(param_grid["direct_ox_count"], errors="coerce").fillna(0)
        + pd.to_numeric(param_grid["cross_talk_count"], errors="coerce").fillna(0)
        + pd.to_numeric(param_grid["backbone_count"], errors="coerce").fillna(0)
    )
    combo_targets = load_combo_targets()

    default_combo_match = param_grid[
        (param_grid["prob_threshold"] == 0.05)
        & (param_grid["string_threshold"] == 0.4)
        & (param_grid["top_k"] == 30)
    ]
    default_combo = default_combo_match.iloc[0] if not default_combo_match.empty else None

    best_combo = select_best_top30_combo(param_grid)
    best_combo_name = str(best_combo["combo"])
    best_combo_genes = combo_targets.get(best_combo_name, [])[:30]
    if len(best_combo_genes) < 30:
        raise ValueError(f"Representative combo {best_combo_name} has fewer than 30 targets.")

    rows = [
        summarize_parameter_target_set(
            "fixed_default_top30",
            "Fixed default Top30",
            "One conventional parameter setting: Prob>0.05, STRING>=0.4, degree Top30.",
            default_top30,
            robust_by_gene,
            q1_top30_genes,
            q1_all_genes,
            default_combo,
        ),
        summarize_parameter_target_set(
            "single_best_parameter_top30",
            "Single-best parameter Top30",
            "Best TopK=30 parameter setting selected by oxidative-relevant enrichment terms, then degree Top30.",
            best_combo_genes,
            robust_by_gene,
            q1_top30_genes,
            q1_all_genes,
            best_combo,
        ),
        summarize_parameter_target_set(
            "frequency_top30",
            "Frequency Top30",
            "Top30 by target recurrence frequency across the full 112-combination parameter grid.",
            frequency_top30,
            robust_by_gene,
            q1_top30_genes,
            q1_all_genes,
        ),
        summarize_parameter_target_set(
            "q1_composite_top30",
            "Q1 robust Top30",
            "Top30 Q1 targets ranked by frequency, normalized degree, average rank, and betweenness.",
            q1_top30,
            robust_by_gene,
            q1_top30_genes,
            q1_all_genes,
        ),
    ]
    return pd.DataFrame(rows)


def build_compound_perturbation_ablation() -> pd.DataFrame:
    influence_df = safe_read_csv(COMPOUND_INFLUENCE)
    required = [
        "mean_target_retention_vs_full_combo",
        "mean_pathway_retention_vs_full_combo",
        "q1_target_retention_vs_full",
        "q1_pathway_retention_vs_full",
        "mean_retention",
        "compound_influence_score",
    ]
    for column in required:
        influence_df[column] = pd.to_numeric(influence_df[column], errors="coerce")
    output = influence_df.sort_values("compound_influence_score", ascending=False).copy()
    output.insert(0, "influence_rank", range(1, len(output) + 1))
    output["mean_target_loss_vs_full_combo"] = 1.0 - output["mean_target_retention_vs_full_combo"]
    output["mean_pathway_loss_vs_full_combo"] = 1.0 - output["mean_pathway_retention_vs_full_combo"]
    output["q1_target_loss_vs_full"] = 1.0 - output["q1_target_retention_vs_full"]
    output["q1_pathway_loss_vs_full"] = 1.0 - output["q1_pathway_retention_vs_full"]
    return output[
        [
            "influence_rank",
            "drop_compound",
            "drop_compound_name",
            "ok_combinations",
            "skipped_combinations",
            "mean_target_retention_vs_full_combo",
            "mean_target_loss_vs_full_combo",
            "mean_pathway_retention_vs_full_combo",
            "mean_pathway_loss_vs_full_combo",
            "q1_target_retention_vs_full",
            "q1_target_loss_vs_full",
            "q1_pathway_retention_vs_full",
            "q1_pathway_loss_vs_full",
            "mean_retention",
            "compound_influence_score",
        ]
    ]


def plot_parameter_top30_ablation(parameter_df: pd.DataFrame, output_dir: Path) -> None:
    if parameter_df.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    ordered = parameter_df.copy()
    labels = ["Default", "Best", "Freq", "Q1"]
    panels = [
        ("mean_parameter_frequency", "Mean parameter frequency"),
        ("min_parameter_frequency", "Minimum parameter frequency"),
        ("q1_top30_overlap_rate", "Overlap with Q1 Top30"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.5), dpi=160)
    colors = [OKABE_ITO["gray"], OKABE_ITO["orange"], OKABE_ITO["blue"], OKABE_ITO["green"]]
    for ax, (metric, title) in zip(axes, panels):
        values = ordered[metric].astype(float).fillna(0.0).tolist()
        ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.8)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.25, linestyle="--")
    fig.tight_layout()
    fig.savefig(figure_dir / "fig3_parameter_top30_ablation.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig3_parameter_top30_ablation.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "fig3_parameter_top30_ablation.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_compound_perturbation_ablation(compound_df: pd.DataFrame, output_dir: Path) -> None:
    if compound_df.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    plot_df = compound_df.head(10).iloc[::-1].copy()
    labels = [
        str(name)
        .replace("(3beta,9xi)-3-(beta-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide", "Cardiac glycoside")
        for name in plot_df["drop_compound_name"]
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=160)
    ax.barh(labels, plot_df["compound_influence_score"].astype(float), color=OKABE_ITO["red"])
    ax.set_xlabel("Influence score (1 - mean retention)")
    ax.grid(axis="x", alpha=0.25, linestyle="--")
    fig.tight_layout()
    fig.savefig(figure_dir / "fig4_compound_perturbation_ablation.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig4_compound_perturbation_ablation.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "fig4_compound_perturbation_ablation.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_external_validation(df: pd.DataFrame, output_dir: Path) -> None:
    if df.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    available_sources = set(df["validation_source"])
    pubmed_sources = sorted(source for source in available_sources if source.startswith("PubMed_ge_"))
    preferred_sources = ["OMIM", *pubmed_sources, "TTD"]
    sources = [source for source in preferred_sources if source in available_sources]
    if not sources:
        sources = [source for source in ["OMIM", "TTD", "AnySource"] if source in available_sources]
    plot_df = df[df["validation_source"].isin(sources)].copy()
    if plot_df.empty:
        return

    target_order = ["traditional_top30", "q1_candidate_robust", "dual_robust", "tier1"]
    source_order = sources
    x = np.arange(len(target_order))
    width = 0.22

    fig, ax = plt.subplots(figsize=(8.2, 4.5), dpi=160)
    palette = [
        OKABE_ITO["blue"],
        OKABE_ITO["orange"],
        OKABE_ITO["green"],
        OKABE_ITO["purple"],
        OKABE_ITO["sky"],
    ]
    for idx, source in enumerate(source_order):
        values = []
        for target in target_order:
            row = plot_df[(plot_df["target_set"] == target) & (plot_df["validation_source"] == source)]
            values.append(float(row["target_support_rate"].iloc[0]) if not row.empty else 0.0)
        ax.bar(
            x + (idx - (len(source_order) - 1) / 2) * width,
            values,
            width,
            label=source,
            color=palette[idx % len(palette)],
        )

    ax.set_ylabel("Overlap rate within target set")
    ax.set_xticks(x)
    ax.set_xticklabels(["Top30", "Q1", "Dual", "Tier1"], rotation=0)
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.25, linestyle="--")
    ax.legend(frameon=False, ncol=len(source_order), loc="upper center", bbox_to_anchor=(0.5, 1.12))
    fig.tight_layout()
    fig.savefig(figure_dir / "fig1_external_source_support.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig1_external_source_support.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "fig1_external_source_support.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_permutation(summary_df: pd.DataFrame, metrics_df: pd.DataFrame, output_dir: Path) -> None:
    if summary_df.empty or metrics_df.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    candidate_metrics = ["OMIM_count", "TTD_count", "AnySource_count"]
    pubmed_metrics = [col for col in metrics_df.columns if col.startswith("PubMed_ge_") and col.endswith("_count")]
    candidate_metrics.extend(pubmed_metrics[:1])

    target_set = "dual_robust" if "dual_robust" in set(metrics_df["target_set"]) else str(metrics_df["target_set"].iloc[0])
    available = [metric for metric in candidate_metrics if metric in metrics_df.columns]
    if not available:
        return

    n_panels = len(available)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.2 * n_panels, 3.4), dpi=160)
    if n_panels == 1:
        axes = [axes]

    subset = metrics_df[metrics_df["target_set"] == target_set].copy()
    observed = subset[subset["is_observed"]].iloc[0]
    null = subset[~subset["is_observed"]]
    for ax, metric in zip(axes, available):
        null_values = null[metric].dropna().astype(float)
        ax.hist(null_values, bins=min(20, max(5, int(null_values.nunique()))), color="#B8BEC8", edgecolor="white")
        ax.axvline(float(observed[metric]), color=OKABE_ITO["red"], linewidth=2.2, label="Observed")
        row = summary_df[(summary_df["target_set"] == target_set) & (summary_df["metric"] == metric)]
        p_text = ""
        if not row.empty:
            p_text = f"p={float(row['empirical_p_greater_equal'].iloc[0]):.3f}"
        ax.set_title(f"{metric}\n{p_text}", fontsize=9)
        ax.set_ylabel("Matched null sets")
        ax.grid(alpha=0.25, linestyle="--")
        ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(figure_dir / "fig2_degree_matched_support_permutation.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig2_degree_matched_support_permutation.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "fig2_degree_matched_support_permutation.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_pubmed_threshold_sensitivity(pubmed_key_df: pd.DataFrame, output_dir: Path) -> None:
    if pubmed_key_df.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    plot_df = pubmed_key_df.copy()
    for column in ["pubmed_threshold", "enrichment_ratio", "fdr_bh"]:
        plot_df[column] = pd.to_numeric(plot_df[column], errors="coerce")
    plot_df["minus_log10_fdr"] = plot_df["fdr_bh"].clip(lower=1e-300).map(lambda value: -math.log10(value))

    target_order = ["traditional_top30", "q1_candidate_robust", "dual_robust"]
    labels = {
        "traditional_top30": "Traditional Top30",
        "q1_candidate_robust": "Parameter-Q1 robust",
        "dual_robust": "Dual-robust",
    }
    colors = {
        "traditional_top30": OKABE_ITO["gray"],
        "q1_candidate_robust": OKABE_ITO["blue"],
        "dual_robust": OKABE_ITO["red"],
    }

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), dpi=160)
    for target_set in target_order:
        current = plot_df[plot_df["target_set"] == target_set].sort_values("pubmed_threshold")
        if current.empty:
            continue
        axes[0].plot(
            current["pubmed_threshold"],
            current["minus_log10_fdr"],
            marker="o",
            linewidth=2,
            color=colors[target_set],
            label=labels[target_set],
        )
        axes[1].plot(
            current["pubmed_threshold"],
            current["enrichment_ratio"],
            marker="o",
            linewidth=2,
            color=colors[target_set],
            label=labels[target_set],
        )
    axes[0].axhline(-math.log10(0.05), color="#333333", linestyle="--", linewidth=1, label="FDR 0.05")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(list(PUBMED_KEY_THRESHOLDS))
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.grid(axis="y", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("-log10(BH-FDR)")
    axes[0].set_xlabel("PubMed co-occurrence threshold")
    axes[1].set_ylabel("Enrichment ratio")
    axes[1].set_xlabel("PubMed co-occurrence threshold")
    axes[0].set_title("Significance across literature thresholds")
    axes[1].set_title("Effect size across literature thresholds")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "fig5_pubmed_cooccurrence_threshold_sensitivity.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig5_pubmed_cooccurrence_threshold_sensitivity.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "fig5_pubmed_cooccurrence_threshold_sensitivity.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_ablation(ablation_df: pd.DataFrame, output_dir: Path) -> None:
    if ablation_df.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    labels = ["Top30", "Q1", "Dual", "Tier1"]
    ordered = ablation_df.set_index("target_set").loc[
        ["traditional_top30", "q1_candidate_robust", "dual_robust", "tier1"]
    ].reset_index()

    panels = [
        ("dual_overlap_rate", "Dual-robust overlap rate"),
        ("mean_parameter_frequency", "Mean parameter frequency"),
        ("mean_loo_q1_retention_rate", "Mean LOO Q1 retention"),
        ("OMIM_count", "OMIM-supported targets"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.5), dpi=160)
    colors = [OKABE_ITO["gray"], OKABE_ITO["blue"], OKABE_ITO["green"], OKABE_ITO["orange"]]
    for ax, (metric, title) in zip(axes, panels):
        values = ordered[metric].astype(float).fillna(0.0).tolist()
        ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.8)
        ax.set_title(title, fontsize=9)
        ax.grid(axis="y", alpha=0.25, linestyle="--")
        ax.tick_params(axis="x", labelrotation=0)
        if metric.endswith("rate") or metric == "mean_parameter_frequency":
            ax.set_ylim(0, 1.05)

    fig.tight_layout()
    fig.savefig(figure_dir / "fig3_ablation_comparison.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig3_ablation_comparison.svg", bbox_inches="tight")
    fig.savefig(figure_dir / "fig3_ablation_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_analysis_report(
    output_dir: Path,
    external_df: pd.DataFrame,
    pubmed_threshold_df: pd.DataFrame,
    pubmed_key_df: pd.DataFrame,
    pubmed_threshold_degree_df: pd.DataFrame,
    permutation_summary: pd.DataFrame,
    context_comparison_df: pd.DataFrame,
    parameter_ablation_df: pd.DataFrame,
    compound_ablation_df: pd.DataFrame,
    pubmed_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    def fmt(value: Any, digits: int = 3) -> str:
        if value is None or pd.isna(value):
            return "NA"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{digits}g}"
        return str(value)

    dual_omim = external_df[
        (external_df["target_set"] == "dual_robust") & (external_df["validation_source"] == "OMIM")
    ]
    q1_omim = external_df[
        (external_df["target_set"] == "q1_candidate_robust") & (external_df["validation_source"] == "OMIM")
    ]
    pubmed_sources = sorted(
        source for source in external_df["validation_source"].dropna().unique() if str(source).startswith("PubMed_ge_")
    )
    pubmed_source = pubmed_sources[0] if pubmed_sources else None
    dual_pubmed = pd.DataFrame()
    q1_pubmed = pd.DataFrame()
    if pubmed_source is not None:
        dual_pubmed = external_df[
            (external_df["target_set"] == "dual_robust")
            & (external_df["validation_source"] == pubmed_source)
        ]
        q1_pubmed = external_df[
            (external_df["target_set"] == "q1_candidate_robust")
            & (external_df["validation_source"] == pubmed_source)
        ]
    pubmed_status = "not requested"
    if not pubmed_df.empty:
        status_counts = pubmed_df["query_status"].value_counts().to_dict()
        pubmed_status = ", ".join(f"{key}={value}" for key, value in status_counts.items())

    report_lines = [
        "# Robustness Validation Bundle",
        "",
        "## Analysis Question",
        "",
        "Do parameter perturbations identify a more stable Top30 target set, and which compounds most perturb the antioxidant network under leave-one-compound-out analysis?",
        "",
        "## Key Findings",
        "",
    ]
    if not dual_omim.empty:
        row = dual_omim.iloc[0]
        report_lines.append(
            f"- Dual-robust targets overlap OMIM-supported oxidative-stress genes at "
            f"{int(row['overlap_n'])}/{int(row['target_n_in_background'])} targets "
            f"(hypergeometric p={fmt(row['p_value'])}, BH-FDR={fmt(row['fdr_bh'])})."
        )
    if not q1_omim.empty:
        row = q1_omim.iloc[0]
        report_lines.append(
            f"- Q1 candidate robust targets overlap OMIM-supported genes at "
            f"{int(row['overlap_n'])}/{int(row['target_n_in_background'])} targets "
            f"(hypergeometric p={fmt(row['p_value'])}, BH-FDR={fmt(row['fdr_bh'])})."
        )
    if not dual_pubmed.empty:
        row = dual_pubmed.iloc[0]
        report_lines.append(
            f"- Dual-robust targets overlap PubMed-supported genes ({pubmed_source}) at "
            f"{int(row['overlap_n'])}/{int(row['target_n_in_background'])} targets "
            f"against the candidate target background set "
            f"(hypergeometric p={fmt(row['p_value'])}, BH-FDR={fmt(row['fdr_bh'])})."
        )
    if not q1_pubmed.empty:
        row = q1_pubmed.iloc[0]
        report_lines.append(
            f"- Q1 candidate robust targets overlap PubMed-supported genes ({pubmed_source}) at "
            f"{int(row['overlap_n'])}/{int(row['target_n_in_background'])} targets "
            f"(hypergeometric p={fmt(row['p_value'])}, BH-FDR={fmt(row['fdr_bh'])})."
        )
    if not pubmed_threshold_df.empty:
        dual_high = pubmed_threshold_df[
            (pubmed_threshold_df["target_set"] == "dual_robust")
            & (pubmed_threshold_df["pubmed_threshold"] == 50)
        ]
        if not dual_high.empty:
            row = dual_high.iloc[0]
            report_lines.append(
                f"- PubMed threshold sensitivity is reported across multiple cutoffs; at >=50 papers, "
                f"dual-robust targets overlap {int(row['overlap_n'])}/{int(row['target_n_in_background'])} "
                f"high-literature genes (BH-FDR={fmt(row['fdr_bh'])})."
            )
    if not pubmed_key_df.empty:
        q1_high = pubmed_key_df[
            (pubmed_key_df["target_set"] == "q1_candidate_robust")
            & (pubmed_key_df["pubmed_threshold"] == 50)
        ]
        dual_strict = pubmed_key_df[
            (pubmed_key_df["target_set"] == "dual_robust")
            & (pubmed_key_df["pubmed_threshold"] == 100)
        ]
        if not q1_high.empty:
            row = q1_high.iloc[0]
            report_lines.append(
                f"- Q1 parameter-robust targets show strong PubMed co-occurrence support at >=50 papers: "
                f"{int(row['overlap_n'])}/{int(row['target_n_in_background'])} targets, "
                f"enrichment={fmt(row['enrichment_ratio'])}, BH-FDR={fmt(row['fdr_bh'])}."
            )
        if not dual_strict.empty:
            row = dual_strict.iloc[0]
            report_lines.append(
                f"- Dual-robust targets are not strongest at the permissive >=10 PubMed threshold, but at >=100 papers "
                f"they overlap {int(row['overlap_n'])}/{int(row['target_n_in_background'])} targets "
                f"(enrichment={fmt(row['enrichment_ratio'])}, BH-FDR={fmt(row['fdr_bh'])})."
            )
    if not pubmed_threshold_degree_df.empty:
        dual_degree = pubmed_threshold_degree_df[
            (pubmed_threshold_degree_df["target_set"] == "dual_robust")
            & (pubmed_threshold_degree_df["pubmed_threshold"] == 100)
        ]
        if not dual_degree.empty:
            row = dual_degree.iloc[0]
            report_lines.append(
                f"- Multi-threshold degree-matched PubMed check confirms dual-robust support at >=100 papers: "
                f"observed={fmt(row['observed'])}, null mean={fmt(row['null_mean'])}, "
                f"empirical p>={fmt(row['empirical_p_greater_equal'])}."
            )

    if not permutation_summary.empty:
        focus = permutation_summary[
            (permutation_summary["target_set"] == "dual_robust")
            & (permutation_summary["metric"] == "OMIM_count")
        ]
        if not focus.empty:
            row = focus.iloc[0]
            report_lines.append(
                f"- Degree-matched null test for dual-robust OMIM support: observed={fmt(row['observed'])}, "
                f"null mean={fmt(row['null_mean'])}, empirical p>={fmt(row['empirical_p_greater_equal'])} "
                f"(n={args.n_permutations}, seed={args.seed})."
            )
        if pubmed_source is not None:
            pubmed_metric = f"{pubmed_source}_count"
            pubmed_focus = permutation_summary[
                (permutation_summary["target_set"] == "dual_robust")
                & (permutation_summary["metric"] == pubmed_metric)
            ]
            if not pubmed_focus.empty:
                row = pubmed_focus.iloc[0]
                report_lines.append(
                    f"- Degree-matched null test for dual-robust PubMed support: observed={fmt(row['observed'])}, "
                    f"null mean={fmt(row['null_mean'])}, empirical p>={fmt(row['empirical_p_greater_equal'])} "
                    f"(n={args.n_permutations}, seed={args.seed})."
                )

    default_top30 = parameter_ablation_df[parameter_ablation_df["method"] == "fixed_default_top30"]
    q1_top30 = parameter_ablation_df[parameter_ablation_df["method"] == "q1_composite_top30"]
    if not default_top30.empty and not q1_top30.empty:
        default_row = default_top30.iloc[0]
        q1_row = q1_top30.iloc[0]
        report_lines.append(
            f"- Parameter Top30 ablation controls target-set size at n=30: Q1 robust Top30 has "
            f"mean parameter frequency {fmt(q1_row['mean_parameter_frequency'])} and minimum frequency "
            f"{fmt(q1_row['min_parameter_frequency'])}, compared with fixed default Top30 mean frequency "
            f"{fmt(default_row['mean_parameter_frequency'])} and minimum frequency {fmt(default_row['min_parameter_frequency'])}."
        )
    if not compound_ablation_df.empty:
        top_compound = compound_ablation_df.iloc[0]
        report_lines.append(
            f"- Compound perturbation ablation ranks `{top_compound['drop_compound_name']}` as the largest "
            f"network driver (influence score={fmt(top_compound['compound_influence_score'])}; "
            f"mean target retention={fmt(top_compound['mean_target_retention_vs_full_combo'])}; "
            f"mean pathway retention={fmt(top_compound['mean_pathway_retention_vs_full_combo'])})."
        )

    report_lines.extend(
        [
            "",
            "## Output Tables",
            "",
            "- `external_target_validation.csv`: source-support hypergeometric tests with BH-FDR against the candidate target background set.",
            "- `degree_matched_support_metrics.csv`: observed and degree-matched null support metrics.",
            "- `degree_matched_support_summary.csv`: empirical p-values and null intervals.",
            "- `parameter_top30_ablation.csv`: fixed default, single-best parameter, frequency Top30, and Q1 robust Top30 comparison with n=30 controlled.",
            "- `ablation_comparison.csv`: compatibility alias for `parameter_top30_ablation.csv`.",
            "- `compound_perturbation_ablation.csv`: leave-one-compound-out retention and influence scores.",
            "- `target_set_context_comparison.csv`: contextual comparison of fixed Top30, all Q1 targets, dual-robust targets, and Tier 1 target; not the primary ablation table.",
            "- `target_pubmed_counts.csv`: target-level PubMed co-occurrence counts when requested or cached.",
            "- `pubmed_threshold_sensitivity.csv`: PubMed support enrichment over multiple count thresholds.",
            "- `pubmed_cooccurrence_key_results.csv`: manuscript-facing PubMed threshold table for Top30, Q1, and dual-robust target sets.",
            "- `pubmed_dual_robust_target_counts.csv`: dual-robust target-level PubMed counts sorted by literature support.",
            "- `pubmed_threshold_degree_matched_summary.csv`: optional multi-threshold degree-matched PubMed support checks when generated.",
            "",
            "## PubMed Status",
            "",
            f"- Mode: `{args.pubmed_mode}`; status: {pubmed_status}.",
            "",
            "## Boundary Conditions",
            "",
            "- GeneCards, OMIM, and TTD are source-support checks, not fully independent validation, because these sources also participate in disease-gene construction.",
            "- PubMed target co-occurrence is independent literature support; short symbols use gene/protein phrases and descriptions instead of bare-symbol matching, and the raw query is stored for audit.",
            "- PubMed co-occurrence should be interpreted as literature support, not experimental validation; high-degree and high-attention hub genes can have inflated literature counts.",
            "- Traditional Top30 remains strong in PubMed co-occurrence, so the PubMed result should support the robustness evidence chain rather than claim that dual-robust targets universally beat hub-based selection.",
            "- Degree matching uses nearest-degree non-observed genes from the available PPI target background set; exact matching is impossible for the highest-degree hubs.",
        ]
    )
    (output_dir / "analysis-report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")


def write_stats_appendix(
    output_dir: Path,
    external_df: pd.DataFrame,
    pubmed_key_df: pd.DataFrame,
    pubmed_threshold_degree_df: pd.DataFrame,
    permutation_summary: pd.DataFrame,
    parameter_ablation_df: pd.DataFrame,
    compound_ablation_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    lines = [
        "# Stats Appendix",
        "",
        "## Tests",
        "",
        "- Source-support enrichment uses a one-sided hypergeometric survival test.",
        "- Multiple source/set tests are adjusted with Benjamini-Hochberg FDR.",
        "- Degree-matched support tests use empirical one-sided p-values with +1 correction.",
        "",
        "## Units",
        "",
        "- Hypergeometric unit: gene within the positive-probability disease-filtered predicted target background set.",
        "- Permutation unit: degree-matched target set sampled from the PPI target background set.",
        "",
        "## Parameters",
        "",
        f"- n_permutations: {args.n_permutations}",
        f"- seed: {args.seed}",
        f"- match_pool_size: {args.match_pool_size}",
        f"- pubmed_threshold: {args.pubmed_threshold}",
        "",
        "## External Validation Summary",
        "",
        frame_to_markdown(
            external_df[
                [
                    "target_set",
                    "validation_source",
                    "target_n_in_background",
                    "validation_n_in_background",
                    "overlap_n",
                    "expected_overlap",
                    "enrichment_ratio",
                    "p_value",
                    "fdr_bh",
                ]
            ]
        ),
        "",
        "## PubMed Co-occurrence Threshold Sensitivity",
        "",
        frame_to_markdown(
            pubmed_key_df[
                [
                    "target_set",
                    "pubmed_threshold",
                    "validation_n_in_background",
                    "target_n_in_background",
                    "overlap_n",
                    "expected_overlap",
                    "target_support_rate",
                    "enrichment_ratio",
                    "p_value",
                    "fdr_bh",
                ]
            ]
        )
        if not pubmed_key_df.empty
        else "No PubMed threshold sensitivity table generated.",
        "",
        "## PubMed Multi-Threshold Degree-Matched Summary",
        "",
        frame_to_markdown(
            pubmed_threshold_degree_df[
                [
                    "target_set",
                    "pubmed_threshold",
                    "observed",
                    "null_mean",
                    "null_q025",
                    "null_q975",
                    "empirical_p_greater_equal",
                    "z_score",
                ]
            ]
        )
        if not pubmed_threshold_degree_df.empty
        else "No multi-threshold degree-matched PubMed summary found.",
        "",
        "## Degree-Matched Summary",
        "",
        frame_to_markdown(permutation_summary) if not permutation_summary.empty else "No permutation summary generated.",
        "",
        "## Parameter Top30 Ablation",
        "",
        frame_to_markdown(
            parameter_ablation_df[
                [
                    "method",
                    "n_targets",
                    "representative_combo",
                    "q1_top30_overlap_n",
                    "mean_parameter_frequency",
                    "min_parameter_frequency",
                    "mean_degree_norm",
                ]
            ]
        )
        if not parameter_ablation_df.empty
        else "No parameter Top30 ablation generated.",
        "",
        "## Compound Perturbation Ablation",
        "",
        frame_to_markdown(
            compound_ablation_df[
                [
                    "influence_rank",
                    "drop_compound_name",
                    "mean_target_retention_vs_full_combo",
                    "mean_pathway_retention_vs_full_combo",
                    "compound_influence_score",
                ]
            ].head(10)
        )
        if not compound_ablation_df.empty
        else "No compound perturbation ablation generated.",
    ]
    (output_dir / "stats-appendix.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_pubmed_snapshot(
    output_dir: Path,
    pubmed_key_df: pd.DataFrame,
    pubmed_threshold_degree_df: pd.DataFrame,
    dual_pubmed_counts_df: pd.DataFrame,
    pubmed_df: pd.DataFrame,
) -> None:
    if pubmed_key_df.empty:
        return
    status = "not available"
    if not pubmed_df.empty and "query_status" in pubmed_df.columns:
        counts = pubmed_df["query_status"].value_counts().to_dict()
        status = ", ".join(f"{key}={value}" for key, value in counts.items())
    lines = [
        "# PubMed Co-occurrence Validation Snapshot",
        "",
        f"Input: `target_pubmed_counts.csv`; query status: {status}.",
        "",
        "## Main Reading",
        "",
        "- Parameter-Q1 robust targets are consistently supported by PubMed co-occurrence.",
        "- Dual-robust targets are weak at the permissive threshold >=10, but become clearly enriched at stricter thresholds >=50/100/200.",
        "- Traditional Top30 remains very strong, so PubMed co-occurrence alone should not be used to claim that the dual framework beats hub-based selection.",
        "",
        "## Key Thresholds",
        "",
        frame_to_markdown(
            pubmed_key_df[
                [
                    "target_set",
                    "pubmed_threshold",
                    "validation_n_in_background",
                    "overlap_n",
                    "expected_overlap",
                    "target_support_rate",
                    "enrichment_ratio",
                    "p_value",
                    "fdr_bh",
                ]
            ]
        ),
        "",
        "## Degree-Matched PubMed Thresholds",
        "",
        frame_to_markdown(
            pubmed_threshold_degree_df[
                [
                    "target_set",
                    "pubmed_threshold",
                    "observed",
                    "null_mean",
                    "null_q025",
                    "null_q975",
                    "empirical_p_greater_equal",
                    "z_score",
                ]
            ]
        )
        if not pubmed_threshold_degree_df.empty
        else "No multi-threshold degree-matched PubMed summary found.",
        "",
        "## Dual-Robust High-Count Targets",
        "",
        frame_to_markdown(dual_pubmed_counts_df[["gene", "pubmed_count"]].head(15))
        if not dual_pubmed_counts_df.empty
        else "No dual-robust PubMed target counts generated.",
    ]
    (output_dir / "pubmed_cooccurrence_validation_snapshot.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def write_figure_catalog(output_dir: Path) -> None:
    lines = [
        "# Figure Catalog",
        "",
        "| Figure | Purpose | Data source | Key observation to check | Caveat |",
        "|---|---|---|---|---|",
        "| `figures/fig1_external_source_support.*` | Compare source-support overlap rates across target-selection strategies. | `external_target_validation.csv` | Robust sets should retain source support after filtering. | GeneCards/OMIM/TTD are support checks, not independent gold standards. |",
        "| `figures/fig2_degree_matched_support_permutation.*` | Compare observed dual-robust support to matched null target sets. | `degree_matched_support_metrics.csv` and `degree_matched_support_summary.csv` | Observed line relative to null distribution indicates whether support exceeds degree effects. | Exact degree matching is limited for high-degree hubs. |",
        "| `figures/fig3_parameter_top30_ablation.*` | Compare Top30 target-selection strategies under a controlled target-set size. | `parameter_top30_ablation.csv` | Q1 robust Top30 should improve parameter-frequency metrics relative to fixed or single-best settings. | This table evaluates parameter robustness, not compound contribution. |",
        "| `figures/fig4_compound_perturbation_ablation.*` | Rank compounds by leave-one-compound-out impact on the antioxidant network. | `compound_perturbation_ablation.csv` | Higher influence scores indicate larger network/pathway retention loss after deleting the compound. | This is a computational network impact score, not a wet-lab activity assay. |",
        "| `figures/fig5_pubmed_cooccurrence_threshold_sensitivity.*` | Show PubMed co-occurrence support across literature thresholds. | `pubmed_cooccurrence_key_results.csv` | Q1 should be broadly supported; dual-robust targets should be interpreted at stricter high-literature thresholds. | PubMed co-occurrence is literature support and remains sensitive to hub-gene research attention. |",
        "| `figures/fig3_ablation_comparison.*` | Legacy contextual comparison of endpoint target sets. | `target_set_context_comparison.csv` | Use only as context; it mixes different target-set sizes and therefore is not the primary ablation. | Traditional Top30 can look better on literature coverage due to hub bias and equal-size mismatch. |",
    ]
    (output_dir / "figure-catalog.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_result_lock_addendum(
    output_dir: Path,
    external_df: pd.DataFrame,
    pubmed_key_df: pd.DataFrame,
    pubmed_threshold_degree_df: pd.DataFrame,
    permutation_summary: pd.DataFrame,
    context_comparison_df: pd.DataFrame,
    parameter_ablation_df: pd.DataFrame,
    compound_ablation_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    addendum = SUPPLEMENTARY_DIR / "result_lock_robustness_validation_addendum.md"
    dual_omim = external_df[
        (external_df["target_set"] == "dual_robust") & (external_df["validation_source"] == "OMIM")
    ]
    q1_omim = external_df[
        (external_df["target_set"] == "q1_candidate_robust") & (external_df["validation_source"] == "OMIM")
    ]
    pubmed_sources = sorted(
        source for source in external_df["validation_source"].dropna().unique() if str(source).startswith("PubMed_ge_")
    )
    pubmed_source = pubmed_sources[0] if pubmed_sources else None

    lines = [
        "# Result Lock Robustness Validation Addendum",
        "",
        "This addendum records validation artifacts generated after `result_lock_v1.md`.",
        "",
        "## Validation Bundle",
        "",
        f"- Output directory: `{output_dir.relative_to(ROOT)}`",
        f"- Degree-matched permutations: {args.n_permutations}",
        f"- Random seed: {args.seed}",
        f"- PubMed mode: `{args.pubmed_mode}`",
        "",
        "## Locked Robustness Validation Highlights",
        "",
    ]
    if not q1_omim.empty:
        row = q1_omim.iloc[0]
        lines.append(
            f"- Q1 targets OMIM overlap: {int(row['overlap_n'])}/{int(row['target_n_in_background'])}, "
            f"BH-FDR={float(row['fdr_bh']):.4g}."
        )
    if not dual_omim.empty:
        row = dual_omim.iloc[0]
        lines.append(
            f"- Dual-robust targets OMIM overlap: {int(row['overlap_n'])}/{int(row['target_n_in_background'])}, "
            f"BH-FDR={float(row['fdr_bh']):.4g}."
        )
    if pubmed_source is not None:
        dual_pubmed = external_df[
            (external_df["target_set"] == "dual_robust")
            & (external_df["validation_source"] == pubmed_source)
        ]
        if not dual_pubmed.empty:
            row = dual_pubmed.iloc[0]
            lines.append(
                f"- Dual-robust targets PubMed overlap ({pubmed_source}): "
                f"{int(row['overlap_n'])}/{int(row['target_n_in_background'])}, "
                f"BH-FDR={float(row['fdr_bh']):.4g}."
            )
    focus = permutation_summary[
        (permutation_summary["target_set"] == "dual_robust")
        & (permutation_summary["metric"] == "OMIM_count")
    ]
    if not focus.empty:
        row = focus.iloc[0]
        lines.append(
            f"- Dual-robust OMIM degree-matched empirical p>=: "
            f"{float(row['empirical_p_greater_equal']):.4g}."
        )
    if pubmed_source is not None:
        pubmed_focus = permutation_summary[
            (permutation_summary["target_set"] == "dual_robust")
            & (permutation_summary["metric"] == f"{pubmed_source}_count")
        ]
        if not pubmed_focus.empty:
            row = pubmed_focus.iloc[0]
            lines.append(
                f"- Dual-robust PubMed degree-matched empirical p>=: "
                f"{float(row['empirical_p_greater_equal']):.4g}."
            )
    if not pubmed_key_df.empty:
        dual_high = pubmed_key_df[
            (pubmed_key_df["target_set"] == "dual_robust")
            & (pubmed_key_df["pubmed_threshold"] == 100)
        ]
        q1_high = pubmed_key_df[
            (pubmed_key_df["target_set"] == "q1_candidate_robust")
            & (pubmed_key_df["pubmed_threshold"] == 50)
        ]
        if not q1_high.empty:
            row = q1_high.iloc[0]
            lines.append(
                f"- Q1 targets PubMed >=50 overlap: {int(row['overlap_n'])}/{int(row['target_n_in_background'])}, "
                f"BH-FDR={float(row['fdr_bh']):.4g}."
            )
        if not dual_high.empty:
            row = dual_high.iloc[0]
            lines.append(
                f"- Dual-robust targets PubMed >=100 overlap: {int(row['overlap_n'])}/{int(row['target_n_in_background'])}, "
                f"BH-FDR={float(row['fdr_bh']):.4g}."
            )
    if not pubmed_threshold_degree_df.empty:
        dual_degree = pubmed_threshold_degree_df[
            (pubmed_threshold_degree_df["target_set"] == "dual_robust")
            & (pubmed_threshold_degree_df["pubmed_threshold"] == 100)
        ]
        if not dual_degree.empty:
            row = dual_degree.iloc[0]
            lines.append(
                f"- Dual-robust PubMed >=100 degree-matched empirical p>=: "
                f"{float(row['empirical_p_greater_equal']):.4g}."
            )
    dual = context_comparison_df[context_comparison_df["target_set"] == "dual_robust"]
    if not dual.empty:
        row = dual.iloc[0]
        lines.append(
            f"- Dual-robust mean LOO Q1 retention: {float(row['mean_loo_q1_retention_rate']):.4g}."
        )
    q1_top30 = parameter_ablation_df[parameter_ablation_df["method"] == "q1_composite_top30"]
    if not q1_top30.empty:
        row = q1_top30.iloc[0]
        lines.append(
            f"- Q1 robust Top30 mean parameter frequency: "
            f"{float(row['mean_parameter_frequency']):.4g}; minimum frequency: "
            f"{float(row['min_parameter_frequency']):.4g}."
        )
    if not compound_ablation_df.empty:
        row = compound_ablation_df.iloc[0]
        lines.append(
            f"- Top compound perturbation driver: `{row['drop_compound_name']}`, "
            f"influence score={float(row['compound_influence_score']):.4g}."
        )

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `external_target_validation.csv`",
            "- `degree_matched_support_metrics.csv`",
            "- `degree_matched_support_summary.csv`",
            "- `degree_matched_random_sets.csv`",
            "- `parameter_top30_ablation.csv`",
            "- `ablation_comparison.csv`",
            "- `compound_perturbation_ablation.csv`",
            "- `target_set_context_comparison.csv`",
            "- `target_pubmed_counts.csv`",
            "- `pubmed_threshold_sensitivity.csv`",
            "- `pubmed_cooccurrence_key_results.csv`",
            "- `pubmed_dual_robust_target_counts.csv`",
            "- `pubmed_threshold_degree_matched_summary.csv`",
            "- `analysis-report.md`",
            "- `stats-appendix.md`",
            "- `figure-catalog.md`",
        ]
    )
    addendum.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260503)
    parser.add_argument("--match-pool-size", type=int, default=80)
    parser.add_argument(
        "--pubmed-mode",
        choices=["skip", "cache", "query"],
        default="cache",
        help="cache never uses network; query calls NCBI E-utilities.",
    )
    parser.add_argument(
        "--pubmed-scope",
        choices=["target_sets", "background", "universe"],
        default="target_sets",
        help="Genes to query when pubmed-mode=query. Use background for the 286-gene candidate target background set.",
    )
    parser.add_argument("--pubmed-threshold", type=int, default=10)
    parser.add_argument("--pubmed-delay-seconds", type=float, default=0.35)
    parser.add_argument("--context-term", default="oxidative stress")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading locked target sets and validation sources...")
    background = load_candidate_target_background()
    source_sets = load_source_sets()
    target_sets = load_target_sets()
    degree_table = load_degree_table(background)
    descriptions = load_gene_descriptions()

    if args.pubmed_scope in {"background", "universe"}:
        pubmed_genes = sorted(background)
    else:
        pubmed_genes = sorted(set().union(*(set(item.genes) for item in target_sets)))

    print(f"Candidate target background genes: {len(background)}")
    print(f"PPI degree genes in background: {len(degree_table)}")
    print(f"PubMed mode: {args.pubmed_mode}; genes considered: {len(pubmed_genes)}")
    pubmed_df = collect_pubmed_counts(
        pubmed_genes,
        descriptions,
        output_dir,
        args.pubmed_mode,
        args.context_term,
        args.pubmed_delay_seconds,
    )
    pubmed_df.to_csv(output_dir / "target_pubmed_counts.csv", index=False, encoding="utf-8-sig")

    print("Running source-support hypergeometric validation...")
    external_df = run_external_validation(
        target_sets,
        source_sets,
        background,
        pubmed_df,
        args.pubmed_threshold,
    )
    external_df.to_csv(output_dir / "external_target_validation.csv", index=False, encoding="utf-8-sig")
    pubmed_threshold_df = build_pubmed_threshold_sensitivity(
        target_sets,
        pubmed_df,
        background,
        thresholds=[1, 5, 10, 20, 50, 100, 200, 500, 1000],
    )
    pubmed_threshold_df.to_csv(
        output_dir / "pubmed_threshold_sensitivity.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pubmed_key_df = build_pubmed_key_results(pubmed_threshold_df)
    pubmed_key_df.to_csv(
        output_dir / "pubmed_cooccurrence_key_results.csv",
        index=False,
        encoding="utf-8-sig",
    )
    dual_pubmed_counts_df = build_dual_pubmed_counts(pubmed_df, target_sets)
    dual_pubmed_counts_df.to_csv(
        output_dir / "pubmed_dual_robust_target_counts.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(f"Running degree-matched support permutations: n={args.n_permutations}")
    metrics_df, random_sets_df, permutation_summary = run_degree_matched_permutation(
        target_sets,
        source_sets,
        degree_table,
        args.n_permutations,
        args.seed,
        args.match_pool_size,
        pubmed_df,
        args.pubmed_threshold,
    )
    metrics_df.to_csv(output_dir / "degree_matched_support_metrics.csv", index=False, encoding="utf-8-sig")
    random_sets_df.to_csv(output_dir / "degree_matched_random_sets.csv", index=False, encoding="utf-8-sig")
    permutation_summary.to_csv(output_dir / "degree_matched_support_summary.csv", index=False, encoding="utf-8-sig")
    save_json(output_dir / "degree_matched_support_summary.json", {"rows": permutation_summary.to_dict("records")})
    pubmed_threshold_degree_df = load_pubmed_threshold_degree_summary(output_dir)

    print("Building context comparison and separated ablation tables...")
    context_comparison_df = build_ablation_table(
        target_sets,
        source_sets,
        degree_table,
        pubmed_df,
        args.pubmed_threshold,
    )
    parameter_ablation_df = build_parameter_top30_ablation()
    compound_ablation_df = build_compound_perturbation_ablation()
    context_comparison_df.to_csv(
        output_dir / "target_set_context_comparison.csv",
        index=False,
        encoding="utf-8-sig",
    )
    parameter_ablation_df.to_csv(
        output_dir / "parameter_top30_ablation.csv",
        index=False,
        encoding="utf-8-sig",
    )
    parameter_ablation_df.to_csv(
        output_dir / "ablation_comparison.csv",
        index=False,
        encoding="utf-8-sig",
    )
    compound_ablation_df.to_csv(
        output_dir / "compound_perturbation_ablation.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("Writing figures and reports...")
    plot_external_validation(external_df, output_dir)
    plot_permutation(permutation_summary, metrics_df, output_dir)
    plot_parameter_top30_ablation(parameter_ablation_df, output_dir)
    plot_compound_perturbation_ablation(compound_ablation_df, output_dir)
    plot_pubmed_threshold_sensitivity(pubmed_key_df, output_dir)
    plot_ablation(context_comparison_df, output_dir)
    write_pubmed_snapshot(
        output_dir,
        pubmed_key_df,
        pubmed_threshold_degree_df,
        dual_pubmed_counts_df,
        pubmed_df,
    )
    write_analysis_report(
        output_dir,
        external_df,
        pubmed_threshold_df,
        pubmed_key_df,
        pubmed_threshold_degree_df,
        permutation_summary,
        context_comparison_df,
        parameter_ablation_df,
        compound_ablation_df,
        pubmed_df,
        args,
    )
    write_stats_appendix(
        output_dir,
        external_df,
        pubmed_key_df,
        pubmed_threshold_degree_df,
        permutation_summary,
        parameter_ablation_df,
        compound_ablation_df,
        args,
    )
    write_figure_catalog(output_dir)
    write_result_lock_addendum(
        output_dir,
        external_df,
        pubmed_key_df,
        pubmed_threshold_degree_df,
        permutation_summary,
        context_comparison_df,
        parameter_ablation_df,
        compound_ablation_df,
        args,
    )

    print("\nRobustness validation bundle written:")
    print(f"  {output_dir}")
    print(f"  {SUPPLEMENTARY_DIR / 'result_lock_robustness_validation_addendum.md'}")


if __name__ == "__main__":
    main()
