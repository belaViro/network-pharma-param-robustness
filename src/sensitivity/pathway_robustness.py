"""
Functional enrichment-term robustness classification utilities.

The target-level robustness analysis uses occurrence frequency and relative
network importance. For GO BP/KEGG functional enrichment terms, we use the same structure:
occurrence frequency for robustness and average -log10(adjusted p-value) for
enrichment strength.

Historical column and file names use "pathway" for compatibility. In this
project that name means KEGG pathways plus GO BP biological-process terms.
"""

from __future__ import annotations

from typing import Any

import pandas as pd


QUADRANT_CANDIDATE = "Q1: 候选稳健功能条目"
QUADRANT_LOW_FREQ_HIGH_SIG = "Q2: 低频高显著"
QUADRANT_STABLE_WEAK_EFFECT = "Q3: 稳健但效应弱"
QUADRANT_LOW_PRIORITY = "Q4: 可忽略"


def classify_pathway_robustness(rank_stats_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add robustness labels to GO BP/KEGG functional-term rank-stability statistics.

    Parameters
    ----------
    rank_stats_df:
        DataFrame produced by the functional-term rank stability step. Required
        columns: frequency, avg_logp, avg_rank, avg_ef.

    Returns
    -------
    tuple[pd.DataFrame, dict]
        Augmented DataFrame and a JSON-serializable summary.
    """
    required_cols = {"pathway", "frequency", "avg_logp", "avg_rank", "avg_ef"}
    missing = required_cols - set(rank_stats_df.columns)
    if missing:
        raise ValueError(f"Missing required pathway robustness columns: {sorted(missing)}")

    df = rank_stats_df.copy()
    frequency_median = float(df["frequency"].median())
    avg_logp_median = float(df["avg_logp"].median())

    def classify_quadrant(row: pd.Series) -> str:
        high_frequency = row["frequency"] >= frequency_median
        high_significance = row["avg_logp"] >= avg_logp_median
        if high_frequency and high_significance:
            return QUADRANT_CANDIDATE
        if not high_frequency and high_significance:
            return QUADRANT_LOW_FREQ_HIGH_SIG
        if high_frequency and not high_significance:
            return QUADRANT_STABLE_WEAK_EFFECT
        return QUADRANT_LOW_PRIORITY

    df["is_absolute_robust"] = df["frequency"] == 1.0
    df["robustness_quadrant"] = df.apply(classify_quadrant, axis=1)
    df["is_candidate_robust"] = df["robustness_quadrant"] == QUADRANT_CANDIDATE

    candidate_order = df[df["is_candidate_robust"]].sort_values(
        ["frequency", "avg_logp", "avg_rank", "avg_ef"],
        ascending=[False, False, True, False],
    )
    candidate_ranks = {
        pathway: rank
        for rank, pathway in enumerate(candidate_order["pathway"].tolist(), start=1)
    }
    df["candidate_rank"] = df["pathway"].map(candidate_ranks).astype("Int64")

    output_order = [
        "is_absolute_robust",
        "is_candidate_robust",
        "candidate_rank",
        "robustness_quadrant",
        "pathway",
        "source",
        "pathway_category",
        "frequency",
        "count",
        "avg_logp",
        "avg_rank",
        "min_rank",
        "max_rank",
        "rank_std",
        "avg_ef",
    ]
    existing_order = [col for col in output_order if col in df.columns]
    remaining = [col for col in df.columns if col not in existing_order]
    df = df[existing_order + remaining]
    df = df.sort_values(
        ["is_candidate_robust", "candidate_rank", "frequency", "avg_logp", "avg_rank"],
        ascending=[False, True, False, False, True],
        na_position="last",
    )

    quadrant_summary: dict[str, dict[str, Any]] = {}
    for quadrant in [
        QUADRANT_CANDIDATE,
        QUADRANT_LOW_FREQ_HIGH_SIG,
        QUADRANT_STABLE_WEAK_EFFECT,
        QUADRANT_LOW_PRIORITY,
    ]:
        sub = df[df["robustness_quadrant"] == quadrant]
        quadrant_summary[quadrant] = {
            "count": int(len(sub)),
            "percentage": round(float(len(sub) / len(df) * 100), 2) if len(df) else 0.0,
            "avg_frequency": round(float(sub["frequency"].mean()), 4) if len(sub) else 0.0,
            "avg_logp": round(float(sub["avg_logp"].mean()), 4) if len(sub) else 0.0,
        }

    absolute_df = df[df["is_absolute_robust"]].sort_values(
        ["avg_logp", "avg_rank"], ascending=[False, True]
    )
    candidate_df = df[df["is_candidate_robust"]].sort_values(
        ["candidate_rank"], ascending=True
    )

    summary = {
        "method": "Functional-term robustness = occurrence frequency x average -log10(adjusted p-value)",
        "compatibility_note": (
            "Historical pathway_* names refer to mixed KEGG pathways and GO BP biological-process terms; "
            "GO BP terms are not pathways."
        ),
        "thresholds": {
            "frequency_median": round(frequency_median, 4),
            "avg_logp_median": round(avg_logp_median, 4),
        },
        "total_pathways": int(len(df)),
        "absolute_robust_count": int(len(absolute_df)),
        "candidate_robust_count": int(len(candidate_df)),
        "quadrant_summary": quadrant_summary,
        "absolute_robust_pathways": absolute_df[
            ["pathway", "source", "pathway_category", "frequency", "count", "avg_logp", "avg_rank", "avg_ef"]
        ].to_dict(orient="records"),
        "candidate_robust_pathways_top20": candidate_df.head(20)[
            [
                "candidate_rank",
                "pathway",
                "source",
                "pathway_category",
                "frequency",
                "count",
                "avg_logp",
                "avg_rank",
                "avg_ef",
            ]
        ].to_dict(orient="records"),
    }
    return df, summary
