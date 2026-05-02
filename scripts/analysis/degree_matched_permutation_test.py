"""
Degree-matched permutation pilot for Top30 candidate robust targets.

This script tests whether enrichment from the observed Top30 candidate robust
targets is stronger than enrichment from degree-matched random target sets.

Outputs:
  src/outputs/sensitivity/degree_matched_permutation/
    - permutation_metrics.csv
    - random_sets.csv
    - permutation_summary.json
    - figures/fig_degree_matched_permutation_pilot.png/.svg
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.config import BACKGROUND_GENES
from src.sensitivity.classification import classify_pathway


SENSITIVITY_DIR = ROOT / "src" / "outputs" / "sensitivity"
DEFAULT_OUTPUT_DIR = SENSITIVITY_DIR / "degree_matched_permutation"

TARGET_ROBUSTNESS = SENSITIVITY_DIR / "target_robustness.csv"
NODE_DEGREES = SENSITIVITY_DIR / "ppi_network_full" / "node_degrees.csv"
REAL_ENRICHMENT = (
    SENSITIVITY_DIR
    / "top30_candidate_targets_enrichment"
    / "top30_candidate_targets_enrichment_GO_BP_KEGG.csv"
)

LIBRARIES = {
    "GO_BP": "GO_Biological_Process_2021",
    "KEGG": "KEGG_2021_Human",
}

MECHANISM_TERMS = {
    "AGE-RAGE signaling pathway in diabetic complications",
    "TNF signaling pathway",
    "HIF-1 signaling pathway",
    "PI3K-Akt signaling pathway",
    "MAPK cascade (GO:0000165)",
    "response to reactive oxygen species (GO:0000302)",
    "cellular response to oxidative stress (GO:0034599)",
    "cellular response to reactive oxygen species (GO:0034614)",
    "Lipid and atherosclerosis",
    "Fluid shear stress and atherosclerosis",
    "cytokine-mediated signaling pathway (GO:0019221)",
    "cellular response to cytokine stimulus (GO:0071345)",
    "transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169)",
}


def load_top_targets(k: int) -> List[str]:
    df = pd.read_csv(TARGET_ROBUSTNESS)
    q1 = df[df["quadrant"].str.contains("Q1", na=False)].copy()
    q1 = q1.sort_values(
        ["frequency", "avg_degree_norm", "avg_degree_rank", "avg_betweenness"],
        ascending=[False, False, True, False],
    )
    return q1.head(k)["gene"].tolist()


def load_degree_table() -> pd.DataFrame:
    df = pd.read_csv(NODE_DEGREES)
    df = df.rename(columns={"Gene": "gene", "Degree": "degree"})
    df["gene"] = df["gene"].astype(str)
    df["degree"] = pd.to_numeric(df["degree"], errors="coerce").fillna(0).astype(float)
    return df[["gene", "degree"]].drop_duplicates("gene")


def weighted_choice(candidates: pd.DataFrame, target_degree: float, rng: random.Random) -> str:
    distances = (candidates["degree"] - target_degree).abs().to_numpy(dtype=float)
    weights = 1.0 / np.power(distances + 1.0, 1.25)
    total = float(weights.sum())
    pick = rng.random() * total
    cumulative = 0.0
    for gene, weight in zip(candidates["gene"].tolist(), weights):
        cumulative += float(weight)
        if cumulative >= pick:
            return gene
    return str(candidates["gene"].iloc[-1])


def sample_degree_matched_set(
    target_degrees: pd.DataFrame,
    universe: pd.DataFrame,
    excluded_genes: set[str],
    rng: random.Random,
    match_pool_size: int,
) -> List[Dict]:
    """Sample one independent null set with nearest-degree matching.

    Top30 observed targets are excluded from the null set. Because the observed
    set contains very high-degree hubs, exact matching is impossible; this
    nearest-neighbor design records degree diagnostics so the limitation is
    explicit.
    """
    selected: set[str] = set()
    rows: List[Dict] = []

    shuffled = target_degrees.sample(
        frac=1.0,
        random_state=rng.randrange(1_000_000_000),
    )
    available = universe[~universe["gene"].isin(excluded_genes)].copy()

    for _, target in shuffled.iterrows():
        target_gene = str(target["gene"])
        target_degree = float(target["degree"])
        pool = available[~available["gene"].isin(selected)].copy()
        pool["distance"] = (pool["degree"] - target_degree).abs()
        pool = pool.nsmallest(match_pool_size, "distance")
        if pool.empty:
            raise RuntimeError("No candidates left for degree-matched sampling.")

        gene = weighted_choice(pool, target_degree, rng)
        gene_degree = float(pool.loc[pool["gene"] == gene, "degree"].iloc[0])
        selected.add(gene)
        rows.append(
            {
                "matched_target": target_gene,
                "matched_target_degree": target_degree,
                "gene": gene,
                "degree": gene_degree,
                "degree_distance": abs(gene_degree - target_degree),
            }
        )

    return rows


def cache_key(genes: Sequence[str], pvalue_cutoff: float, top_n: int) -> str:
    payload = ";".join(sorted(genes)) + f"|{pvalue_cutoff}|{top_n}|GO_BP_KEGG"
    return hashlib.md5(payload.encode("utf-8")).hexdigest()[:16]


def save_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)


def run_enrichment_go_bp_kegg(
    genes: Sequence[str],
    cache_dir: Path,
    pvalue_cutoff: float,
    top_n: int,
    sleep_seconds: float,
    retries: int,
    retry_sleep_seconds: float,
) -> Dict[str, pd.DataFrame]:
    key = cache_key(genes, pvalue_cutoff, top_n)
    cache_file = cache_dir / f"{key}.json"
    if cache_file.exists():
        with open(cache_file, "r", encoding="utf-8") as fh:
            cached = json.load(fh)
        return {source: pd.DataFrame(rows) for source, rows in cached.items()}

    cache_dir.mkdir(parents=True, exist_ok=True)
    significant: Dict[str, pd.DataFrame] = {}
    n_genes = len(genes)

    for source, library in LIBRARIES.items():
        last_error: Exception | None = None
        enr = None
        for attempt in range(1, retries + 2):
            try:
                enr = gp.enrichr(
                    gene_list=list(genes),
                    gene_sets=library,
                    organism="hsapiens",
                    cutoff=pvalue_cutoff,
                    outdir=None,
                    background=BACKGROUND_GENES,
                )
                break
            except Exception as exc:  # Enrichr can intermittently reject requests.
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
        if len(df) == 0:
            significant[source] = pd.DataFrame()
            continue

        df["k"] = df["Overlap"].apply(lambda x: int(str(x).split("/")[0]))
        df["N"] = df["Overlap"].apply(lambda x: int(str(x).split("/")[1]))
        df["EF"] = (df["k"] / n_genes) / (df["N"] / BACKGROUND_GENES)
        df = df[(df["k"] >= 2) & (df["EF"] >= 3.0)]
        df = df.sort_values("Adjusted P-value").head(top_n)
        significant[source] = df

    save_json(
        cache_file,
        {source: df.to_dict(orient="records") for source, df in significant.items()},
    )
    return significant


def load_real_enrichment() -> pd.DataFrame:
    if not REAL_ENRICHMENT.exists():
        raise FileNotFoundError(f"Missing real enrichment file: {REAL_ENRICHMENT}")
    df = pd.read_csv(REAL_ENRICHMENT)
    return df[df["source"].isin(["GO_BP", "KEGG"])].copy()


def summarize_enrichment(label: str, genes: Sequence[str], enrichment_df: pd.DataFrame) -> Dict:
    rows = []
    for _, row in enrichment_df.iterrows():
        pvalue = float(row.get("adjusted_pvalue", row.get("Adjusted P-value", 1.0)))
        ef = float(row.get("EF", 0.0))
        source = str(row.get("source", ""))
        pathway = str(row.get("pathway", row.get("Term", "")))
        rows.append(
            {
                "source": source,
                "pathway": pathway,
                "adjusted_pvalue": pvalue,
                "logp": -math.log10(max(pvalue, 1e-300)),
                "ef": ef,
                "category": classify_pathway(pathway, source),
            }
        )

    details = pd.DataFrame(rows)
    if details.empty:
        return {
            "sample_id": label,
            "n_genes": len(genes),
            "total_significant": 0,
            "kegg_significant": 0,
            "go_bp_significant": 0,
            "direct_ox": 0,
            "cross_talk": 0,
            "backbone": 0,
            "other": 0,
            "mechanism_hits": 0,
            "mechanism_logp_sum": 0.0,
            "avg_logp": 0.0,
            "max_logp": 0.0,
            "min_adjusted_pvalue": 1.0,
            "avg_ef": 0.0,
        }

    return {
        "sample_id": label,
        "n_genes": len(genes),
        "total_significant": int(len(details)),
        "kegg_significant": int((details["source"] == "KEGG").sum()),
        "go_bp_significant": int((details["source"] == "GO_BP").sum()),
        "direct_ox": int((details["category"] == "direct_ox").sum()),
        "cross_talk": int((details["category"] == "cross_talk").sum()),
        "backbone": int((details["category"] == "signaling_backbone").sum()),
        "other": int((details["category"] == "other").sum()),
        "mechanism_hits": int(details["pathway"].isin(MECHANISM_TERMS).sum()),
        "mechanism_logp_sum": round(
            float(details.loc[details["pathway"].isin(MECHANISM_TERMS), "logp"].sum()),
            6,
        ),
        "avg_logp": round(float(details["logp"].mean()), 6),
        "max_logp": round(float(details["logp"].max()), 6),
        "min_adjusted_pvalue": float(details["adjusted_pvalue"].min()),
        "avg_ef": round(float(details["ef"].mean()), 6),
    }


def enrichment_dict_to_frame(enrichment: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for source, df in enrichment.items():
        if df is None or len(df) == 0:
            continue
        tmp = df.copy()
        tmp["source"] = source
        tmp["pathway"] = tmp["Term"]
        tmp["adjusted_pvalue"] = tmp["Adjusted P-value"]
        frames.append(tmp)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def empirical_pvalue(real_value: float, null_values: Iterable[float]) -> float:
    values = list(null_values)
    return (sum(v >= real_value for v in values) + 1) / (len(values) + 1)


def percentile(real_value: float, null_values: Iterable[float]) -> float:
    values = list(null_values)
    if not values:
        return 0.0
    return 100.0 * sum(v < real_value for v in values) / len(values)


def write_summary(metrics: pd.DataFrame, degree_diag: Dict, args: argparse.Namespace, output_dir: Path) -> Dict:
    real = metrics[metrics["sample_id"] == "real_top30"].iloc[0]
    null = metrics[metrics["sample_id"] != "real_top30"].copy()
    metric_names = [
        "mechanism_hits",
        "mechanism_logp_sum",
        "total_significant",
        "kegg_significant",
        "go_bp_significant",
        "avg_logp",
        "max_logp",
        "avg_ef",
        "direct_ox",
        "cross_talk",
        "backbone",
    ]

    comparisons = {}
    for metric in metric_names:
        comparisons[metric] = {
            "real": float(real[metric]),
            "null_mean": round(float(null[metric].mean()), 6),
            "null_sd": round(float(null[metric].std(ddof=1)), 6),
            "null_median": round(float(null[metric].median()), 6),
            "real_percentile": round(percentile(float(real[metric]), null[metric].tolist()), 2),
            "empirical_p_greater_equal": round(
                empirical_pvalue(float(real[metric]), null[metric].tolist()),
                6,
            ),
        }

    summary = {
        "analysis": "degree-matched permutation pilot",
        "n_permutations": int(args.n),
        "seed": int(args.seed),
        "top_k_targets": int(args.top_k),
        "pvalue_cutoff": float(args.pvalue_cutoff),
        "top_n_per_source": int(args.top_n),
        "match_pool_size": int(args.match_pool_size),
        "null_design": (
            "Nearest-degree random target sets from PPI universe, excluding observed Top30 targets; "
            "one matched null gene per observed target, unique genes within each random set."
        ),
        "mechanism_terms": sorted(MECHANISM_TERMS),
        "degree_diagnostics": degree_diag,
        "comparisons": comparisons,
    }
    save_json(output_dir / "permutation_summary.json", summary)
    return summary


def plot_metrics(metrics: pd.DataFrame, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    real = metrics[metrics["sample_id"] == "real_top30"].iloc[0]
    null = metrics[metrics["sample_id"] != "real_top30"].copy()

    panels = [
        ("mechanism_hits", "Mechanism term hits"),
        ("mechanism_logp_sum", "Mechanism logP sum"),
        ("total_significant", "Significant terms"),
        ("kegg_significant", "KEGG pathways"),
        ("avg_logp", "Average -log10(adj. p)"),
        ("avg_ef", "Average enrichment factor"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), dpi=160)
    axes = axes.ravel()
    for ax, (metric, title) in zip(axes, panels):
        ax.hist(null[metric], bins=18, color="#B8BEC8", edgecolor="white")
        ax.axvline(real[metric], color="#C0392B", linewidth=2.2, label="Observed Top30")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("Null sets")
        ax.grid(alpha=0.25, linestyle="--")
        ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Degree-matched permutation pilot (n=100)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(figure_dir / "fig_degree_matched_permutation_pilot.png", bbox_inches="tight")
    fig.savefig(figure_dir / "fig_degree_matched_permutation_pilot.svg", bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100, help="Number of random target sets.")
    parser.add_argument("--top-k", type=int, default=30, help="Observed candidate target count.")
    parser.add_argument("--seed", type=int, default=20260501)
    parser.add_argument("--pvalue-cutoff", type=float, default=0.01)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--match-pool-size", type=int, default=80)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-sleep-seconds", type=float, default=8.0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "enrichment_cache"

    rng = random.Random(args.seed)
    top_targets = load_top_targets(args.top_k)
    degree_table = load_degree_table()
    degree_map = dict(zip(degree_table["gene"], degree_table["degree"]))
    missing = [gene for gene in top_targets if gene not in degree_map]
    if missing:
        raise ValueError(f"Top targets missing degree values: {missing}")

    target_degrees = pd.DataFrame(
        [{"gene": gene, "degree": float(degree_map[gene])} for gene in top_targets]
    )
    excluded = set(top_targets)

    real_enrichment = load_real_enrichment()
    metrics = [summarize_enrichment("real_top30", top_targets, real_enrichment)]
    random_rows: List[Dict] = []

    print(f"Observed Top{args.top_k}: {', '.join(top_targets)}")
    print(f"Universe genes: {len(degree_table)}")
    print(f"Running degree-matched randomization: n={args.n}")

    for i in range(1, args.n + 1):
        sample_rows = sample_degree_matched_set(
            target_degrees=target_degrees,
            universe=degree_table,
            excluded_genes=excluded,
            rng=rng,
            match_pool_size=args.match_pool_size,
        )
        genes = [row["gene"] for row in sample_rows]
        for row in sample_rows:
            row["sample_id"] = f"random_{i:04d}"
        random_rows.extend(sample_rows)

        enrichment = run_enrichment_go_bp_kegg(
            genes=genes,
            cache_dir=cache_dir,
            pvalue_cutoff=args.pvalue_cutoff,
            top_n=args.top_n,
            sleep_seconds=args.sleep_seconds,
            retries=args.retries,
            retry_sleep_seconds=args.retry_sleep_seconds,
        )
        enrichment_df = enrichment_dict_to_frame(enrichment)
        record = summarize_enrichment(f"random_{i:04d}", genes, enrichment_df)
        metrics.append(record)

        if i == 1 or i % 10 == 0 or i == args.n:
            print(
                f"[{i:>3}/{args.n}] "
                f"terms={record['total_significant']}, "
                f"KEGG={record['kegg_significant']}, "
                f"mechanism_hits={record['mechanism_hits']}, "
                f"avg_logp={record['avg_logp']:.2f}"
            )

    metrics_df = pd.DataFrame(metrics)
    random_sets_df = pd.DataFrame(random_rows)
    metrics_df.to_csv(output_dir / "permutation_metrics.csv", index=False)
    random_sets_df.to_csv(output_dir / "random_sets.csv", index=False)

    degree_diag = {
        "observed_mean_degree": round(float(target_degrees["degree"].mean()), 6),
        "observed_median_degree": round(float(target_degrees["degree"].median()), 6),
        "observed_min_degree": round(float(target_degrees["degree"].min()), 6),
        "observed_max_degree": round(float(target_degrees["degree"].max()), 6),
        "null_mean_degree_mean": round(
            float(random_sets_df.groupby("sample_id")["degree"].mean().mean()),
            6,
        ),
        "null_mean_degree_sd": round(
            float(random_sets_df.groupby("sample_id")["degree"].mean().std(ddof=1)),
            6,
        ),
        "null_degree_distance_mean": round(float(random_sets_df["degree_distance"].mean()), 6),
        "null_degree_distance_median": round(float(random_sets_df["degree_distance"].median()), 6),
        "note": (
            "Observed Top30 are extremely high-degree hubs. Exact matching while excluding "
            "Top30 is impossible in this universe, so the null uses nearest-degree non-Top30 genes."
        ),
    }
    summary = write_summary(metrics_df, degree_diag, args, output_dir)
    plot_metrics(metrics_df, output_dir)

    print("\nSummary:")
    for metric, item in summary["comparisons"].items():
        if metric in {"mechanism_hits", "mechanism_logp_sum", "total_significant", "kegg_significant", "avg_logp"}:
            print(
                f"  {metric}: real={item['real']}, "
                f"null_mean={item['null_mean']}, "
                f"percentile={item['real_percentile']}%, "
                f"p>={item['empirical_p_greater_equal']}"
            )
    print(f"\nOutputs written to: {output_dir}")


if __name__ == "__main__":
    main()
