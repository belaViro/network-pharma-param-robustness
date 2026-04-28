"""
Generate publication-quality figures for network pharmacology sensitivity analysis.

Figures:
1. Heatmap: STRING × TopK matrix by significant pathway count, faceted by Prob
2. Target robustness quadrant scatter plot (frequency vs avg_degree)
3. Pathway rank stability (top pathways with rank ranges)
4. Pvalue threshold comparison bar chart
5. Robust vs parametric comparison (Q1 / Tier1 vs 72 combos)

Output: PNG (300 DPI) + SVG to src/outputs/sensitivity/figures/
"""

import json
import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, BoundaryNorm
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore", category=UserWarning)

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
FIGURES_DIR = OUTPUTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

DATA = {
    "grid": OUTPUTS_DIR / "param_grid_results.csv",
    "robustness": OUTPUTS_DIR / "target_robustness.csv",
    "pathway_stability": OUTPUTS_DIR / "pathway_rank_stability.csv",
    "comparison": OUTPUTS_DIR / "robust_vs_param_comparison.json",
}

# ── Style ──────────────────────────────────────────────────────────────
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["SimHei", "Microsoft YaHei", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    }
)

# ── Color palettes ─────────────────────────────────────────────────────
COLOR_PROB = {0.0: "#2166AC", 0.09: "#92C5DE", 0.1: "#F4A582"}
COLOR_QUADRANT = {
    "Q1: 候选稳健": "#2166AC",
    "Q2: 低频高重要": "#E69F00",  # not present in current data
    "Q3: 稳健但边缘": "#B2182B",
    "Q4: 可忽略": "#D3D3D3",
}
COLOR_CATEGORY = {
    "signaling_backbone": "#2166AC",
    "cross_talk": "#E69F00",
    "direct_ox": "#D6604D",
    "other": "#D3D3D3",
}
COLOR_COMPARISON = {"Q1 (n=33)": "#2166AC", "Tier1 (n=13)": "#B2182B", "72 Combos (mean)": "#D3D3D3"}

# ── Helpers ────────────────────────────────────────────────────────────


def save_fig(fig: plt.Figure, name: str):
    """Save figure as PNG + SVG."""
    for ext in ["png", "svg"]:
        path = FIGURES_DIR / f"{name}.{ext}"
        fig.savefig(
            path,
            dpi=300 if ext == "png" else 150,
            bbox_inches="tight",
            pad_inches=0.12,
        )
        print(f"  Saved: {path}")


def clean_pathway_name(name: str, max_len: int = 46) -> str:
    """Make pathway names readable in compact figures."""
    import re

    name = re.sub(r"\s*\(GO:\d+\)", "", str(name))
    replacements = {
        "transmembrane receptor protein tyrosine kinase": "TM receptor tyrosine kinase",
        "positive regulation of": "+ regulation of",
        "regulation of": "reg. of",
        "cellular response to": "cell response to",
        "phosphatidylinositol 3-kinase": "PI3K",
        "protein serine/threonine kinase": "ser/thr kinase",
        "protein kinase B": "AKT",
        "signal transduction": "signaling",
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    if len(name) > max_len:
        name = name[: max_len - 3].rstrip() + "..."
    return name


# ══════════════════════════════════════════════════════════════════════
# FIGURE 1 — Heatmap: STRING × TopK faceted by Prob
# ══════════════════════════════════════════════════════════════════════


def fig1_heatmap(df: pd.DataFrame):
    """Faceted heatmap: rows=STRING, cols=TopK, color=significant_pathways, panels=Prob.

    Each Prob panel gets its own color scale (vmin/vmax per panel) so that
    the gradient is visible within each panel — P_prob=0.2 was excluded
    from the grid as it falls into the enrichment failure zone.
    """
    prob_vals = sorted(df["prob_threshold"].unique())
    string_vals = sorted(df["string_threshold"].unique())
    topk_vals = sorted(df["top_k"].unique())

    n_probs = len(prob_vals)
    ncols = 2
    nrows = int(np.ceil(n_probs / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5.5 * ncols, 4 * nrows)
    )
    fig.subplots_adjust(top=0.92)
    axes_flat = axes.flatten() if n_probs > 1 else [axes]

    for ax, prob in zip(axes_flat, prob_vals):
        sub = df[df["prob_threshold"] == prob].pivot_table(
            index="string_threshold", columns="top_k", values="significant_pathways"
        )
        # Ensure consistent ordering
        sub = sub.reindex(index=string_vals, columns=topk_vals)

        # Per-panel color scale (avoid outlier flattening)
        vmin_panel = sub.values.min()
        vmax_panel = sub.values.max()

        im = ax.imshow(sub.values, cmap="RdYlBu_r", aspect="auto",
                       vmin=vmin_panel, vmax=vmax_panel)
        ax.set_title(f"Prob = {prob}", fontsize=10, fontweight="bold")

        # Tick labels
        ax.set_xticks(range(len(topk_vals)))
        ax.set_xticklabels([str(k) for k in topk_vals])
        ax.set_yticks(range(len(string_vals)))
        ax.set_yticklabels([f"{s:.2f}" for s in string_vals])

        # Value annotations (use per-panel midpoint for text color)
        mid_panel = vmin_panel + (vmax_panel - vmin_panel) * 0.5
        for i in range(len(string_vals)):
            for j in range(len(topk_vals)):
                val = sub.values[i, j]
                if np.isnan(val):
                    continue
                color = "white" if val < mid_panel else "black"
                ax.text(j, i, str(int(val)), ha="center", va="center", fontsize=7, color=color)

        # Show the per-panel range as annotation
        ax.text(0.98, 0.02,
                f"范围: [{int(vmin_panel)}, {int(vmax_panel)}]",
                transform=ax.transAxes, fontsize=6.5, color="#555555",
                ha="right", va="bottom", style="italic")

        if prob == prob_vals[0]:
            ax.set_ylabel("STRING阈值", fontsize=9)
        if prob in [prob_vals[-2], prob_vals[-1]]:
            ax.set_xlabel("TopK", fontsize=9)

    # Hide unused subplots
    for ax in axes_flat[n_probs:]:
        ax.set_visible(False)

    # Single colorbar (use last panel's im)
    cbar = fig.colorbar(im, ax=axes_flat[:n_probs].tolist(), shrink=0.6, pad=0.02)
    cbar.set_label("显著通路数（每面板独立刻度）", fontsize=9)

    fig.suptitle("参数敏感性：显著通路数", fontsize=12, y=0.98)
    save_fig(fig, "fig1_heatmap_significant_pathways")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# FIGURE 2 — Target Robustness Quadrant Scatter
# ══════════════════════════════════════════════════════════════════════


def fig2_target_robustness(df: pd.DataFrame):
    """Scatter: frequency (x) vs avg_degree_norm (y), colored by quadrant, labeled for Tier1 genes.

    FIXED: Now uses avg_degree_norm (consistent with run_sensitivity_analysis.py classification)
    so that quadrant colors, background shading, and divider lines are fully aligned.
    Tier1 genes are explicitly labeled.
    """
    from matplotlib.font_manager import FontProperties

    # Compute data-driven medians using NORMALIZED degree (consistent with classification)
    freq_median = df["frequency"].median()
    degree_median = df["avg_degree_norm"].median()

    # Font setup
    font_en = FontProperties(family="Arial", size=9)
    font_en_bold = FontProperties(family="Arial", size=9, weight="bold")
    font_en_small = FontProperties(family="Arial", size=7)
    font_cn = FontProperties(family="SimHei", size=9)
    font_cn_bold = FontProperties(family="SimHei", size=9, weight="bold")
    font_cn_small = FontProperties(family="SimHei", size=7.5)
    font_label = FontProperties(family="Arial", size=6.5, weight="bold")

    fig, ax = plt.subplots(figsize=(8.2, 5.8), constrained_layout=True)

    # Quadrant definitions
    quadrants = {
        "Q1: 候选稳健": {"marker": "o", "size": 55, "zorder": 3},
        "Q2: 低频高重要": {"marker": "s", "size": 55, "zorder": 4},
        "Q3: 稳健但边缘": {"marker": "^", "size": 50, "zorder": 2},
        "Q4: 可忽略": {"marker": "o", "size": 25, "zorder": 1},
    }

    # Background shading — aligned to normalized degree median
    x_max = 1.12
    y_min = -0.05
    y_max = df["avg_degree_norm"].max() * 1.15

    # Compute correct axis-ratio for the degree median line
    y_med_ratio = (degree_median - y_min) / (y_max - y_min)
    y_med_ratio = max(0.0, min(1.0, y_med_ratio))

    ax.axvspan(freq_median, x_max, ymin=y_med_ratio, ymax=1, alpha=0.10, color="#2166AC", zorder=0)
    ax.axvspan(0, freq_median, ymin=y_med_ratio, ymax=1, alpha=0.10, color="#E69F00", zorder=0)
    ax.axvspan(freq_median, x_max, ymin=0, ymax=y_med_ratio, alpha=0.10, color="#B2182B", zorder=0)
    ax.axvspan(0, freq_median, ymin=0, ymax=y_med_ratio, alpha=0.10, color="#D3D3D3", zorder=0)

    # Plot each quadrant
    for qname, qstyle in quadrants.items():
        subset = df[df["quadrant"] == qname]
        color = COLOR_QUADRANT.get(qname, "#999999")
        if len(subset) == 0:
            continue
        display_name = qname
        ax.scatter(
            subset["frequency"],
            subset["avg_degree_norm"],
            c=color,
            marker=qstyle["marker"],
            s=qstyle["size"],
            label=f"{display_name} (n={len(subset)})",
            edgecolors="white",
            linewidth=0.3,
            zorder=qstyle["zorder"],
            alpha=0.85,
        )

    # Label Tier1 genes with separated label anchors to avoid overlap near x=1.
    tier1 = df[df["frequency"] == 1.0].sort_values("avg_degree_norm", ascending=False)
    label_y = tier1["avg_degree_norm"].to_numpy().copy()
    if len(label_y) > 1:
        min_gap = 0.035
        order = np.argsort(-label_y)
        adjusted = label_y.copy()
        last = None
        for idx in order:
            if last is None:
                last = adjusted[idx]
                continue
            if last - adjusted[idx] < min_gap:
                adjusted[idx] = last - min_gap
            last = adjusted[idx]
        label_y = np.clip(adjusted, y_min + 0.03, y_max - 0.03)

    for (_, row), y_lab in zip(tier1.iterrows(), label_y):
        ax.annotate(
            row["gene"],
            (row["frequency"], row["avg_degree_norm"]),
            xytext=(1.018, y_lab),
            textcoords="data",
            fontproperties=font_label,
            color="#333333",
            zorder=5,
            ha="left",
            va="center",
            arrowprops=dict(arrowstyle="-", color="#9aa0a6", lw=0.35),
        )

    # Axis labels
    ax.set_xlabel("跨参数组合出现频率", fontproperties=font_cn, labelpad=6)
    ax.set_ylabel("平均相对连接度 (degree / 组合内最大度值)", fontproperties=font_cn, labelpad=6)
    ax.set_xlim(-0.03, x_max)
    y_max = df["avg_degree_norm"].max() * 1.15
    ax.set_ylim(y_min, y_max)

    # Tick labels
    for label in ax.get_xticklabels():
        label.set_fontproperties(font_en)
    for label in ax.get_yticklabels():
        label.set_fontproperties(font_en)

    # Legend
    leg = ax.legend(frameon=True, fancybox=False, loc="upper left",
                    edgecolor="#cccccc", prop=font_cn_small)
    ax.grid(alpha=0.2, linestyle="--", linewidth=0.3)

    # Median divider lines
    ax.axvline(x=freq_median, color="#555555", linestyle="--", linewidth=0.8, zorder=0, alpha=0.5)
    ax.axhline(y=degree_median, color="#555555", linestyle="--", linewidth=0.8, zorder=0, alpha=0.5)

    # Median annotations
    ax.text(freq_median + 0.01, y_max * 0.93,
            f"频率中位数 = {freq_median:.2f}", fontproperties=font_cn_small,
            color="#555555", va="top")
    ax.text(0.98, degree_median + 0.02,
            f"相对连接度中位数 = {degree_median:.3f}", fontproperties=font_cn_small,
            color="#555555", ha="right", va="bottom")

    # Title
    fig.suptitle("靶点稳健性：二维象限分类（平均相对连接度）",
                 fontproperties=FontProperties(family="SimHei", size=12), y=1.03)
    save_fig(fig, "fig2_target_robustness_quadrant")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# FIGURE 3 — Pathway Rank Stability (lollipop / dumbbell)
# ══════════════════════════════════════════════════════════════════════


def fig3_pathway_rank_stability(df: pd.DataFrame):
    """Horizontal lollipop: top pathways with min/max rank range, colored by category."""
    # Filter to pathways appearing in >= 50% of combos for stability focus
    df_filt = df[df["frequency"] >= 0.5].copy()
    # Sort by median rank
    df_filt = df_filt.sort_values("avg_rank", ascending=True).head(25)

    fig, ax = plt.subplots(figsize=(9.6, 7.6), constrained_layout=True)

    y_pos = np.arange(len(df_filt))
    categories = df_filt["pathway_category"].values
    avg_ranks = df_filt["avg_rank"].values
    min_ranks = df_filt["min_rank"].values
    max_ranks = df_filt["max_rank"].values
    pathway_names = df_filt["pathway"].values
    rank_stds = df_filt["rank_std"].values

    short_names = [clean_pathway_name(p, max_len=42) for p in pathway_names]

    # Draw horizontal bars for rank range
    colors = [COLOR_CATEGORY.get(c, "#D3D3D3") for c in categories]
    edge_colors = [c for c in colors]  # same for dot edges

    for i in range(len(df_filt)):
        # Range line
        ax.plot(
            [min_ranks[i], max_ranks[i]],
            [y_pos[i], y_pos[i]],
            color=colors[i],
            linewidth=2.5,
            alpha=0.5,
            zorder=1,
        )
        # Min dot
        ax.scatter(min_ranks[i], y_pos[i], color=colors[i], s=30, zorder=2, alpha=0.6,
                   edgecolors="white", linewidth=0.3)
        # Max dot
        ax.scatter(max_ranks[i], y_pos[i], color=colors[i], s=30, zorder=2, alpha=0.6,
                   edgecolors="white", linewidth=0.3)
        # Average rank (main dot)
        ax.scatter(avg_ranks[i], y_pos[i], color=colors[i], s=60, zorder=3,
                   edgecolors="white", linewidth=0.5)

    # Category legend
    legend_patches = [
        mpatches.Patch(color=COLOR_CATEGORY["signaling_backbone"], label="信号主干"),
        mpatches.Patch(color=COLOR_CATEGORY["cross_talk"], label="交叉调控"),
        mpatches.Patch(color=COLOR_CATEGORY["other"], label="其他"),
    ]

    ax.set_yticks(y_pos)
    ax.set_yticklabels(short_names, fontsize=7.0)
    ax.tick_params(axis="y", pad=2)
    ax.set_xlabel("Rank range (min - mean - max)", fontsize=10)
    ax.set_xlim(0, max_ranks.max() + 2)
    ax.invert_yaxis()
    ax.legend(handles=legend_patches, frameon=True, fancybox=False, fontsize=8,
              loc="lower right", edgecolor="#cccccc")
    ax.grid(axis="x", alpha=0.3, linestyle="--", linewidth=0.4)

    # Annotate rank_std for the most variable pathways
    for i, std in enumerate(rank_stds):
        if std > 5:
            ax.annotate(
                f"σ={std:.1f}",
                (max_ranks[i] + 0.5, y_pos[i]),
                fontsize=6,
                color="#666666",
                va="center",
            )

    fig.suptitle("通路排名稳定性", fontsize=12, y=0.99)
    save_fig(fig, "fig3_pathway_rank_stability")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# FIGURE 4 — Pvalue Threshold Comparison Bar Chart
# ══════════════════════════════════════════════════════════════════════


def fig4_pvalue_comparison(df_full_pathway=None):
    """Bar chart comparing P=0.01 vs P=0.05 on enrichment metrics.

    If pathway_details CSV is available, reads actual data.
    Otherwise falls back to the grid-level comparison from reports.
    """
    # Attempt to load actual pvalue comparison data
    # Build metrics from grid-level data — backbone/cross_talk/direct_ox/other
    # are computed in run_sensitivity_analysis.py using analyze_enrichment()
    # which only runs for P=0.01. For the paper, we present the
    # P=0.01 vs P=0.05 comparison that was done on Q1 targets.
    #
    # These values are derived from the sensitivity analysis report (section 7)
    # which runs cached_enrichment(q1_genes, pval=0.01) and cached_enrichment(q1_genes, pval=0.05)
    # separately.
    metrics_count = ["总通路数", "信号主干", "交叉调控", "直接氧化应激", "其他"]
    metrics_ef = ["Avg EF"]

    # Values from the sensitivity analysis MD report — the actual Q1 enrichment results
    p001_counts = [64, 8, 3, 0, 53]
    p005_counts = [77, 8, 3, 0, 66]
    p001_ef = 39.84
    p005_ef = 41.96

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), gridspec_kw={"width_ratios": [3, 1]})

    # Left: side-by-side bars for count metrics only
    ax = axes[0]
    x = np.arange(len(metrics_count))
    width = 0.35

    bars1 = ax.bar(x - width / 2, p001_counts, width, label="P = 0.01", color="#2166AC",
                   edgecolor="white", linewidth=0.3)
    bars2 = ax.bar(x + width / 2, p005_counts, width, label="P = 0.05", color="#D6604D",
                   edgecolor="white", linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(metrics_count, fontsize=8)
    ax.set_ylabel("数量", fontsize=10)
    ax.legend(frameon=True, fancybox=False, fontsize=8, edgecolor="#cccccc")
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.4)

    # Value labels on bars
    for bars_list, color in [(bars1, "#2166AC"), (bars2, "#D6604D")]:
        for bar in bars_list:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5, str(int(h)),
                    ha="center", va="bottom", fontsize=7, color=color)

    # Annotate EF separately (it's on a different scale)
    ax.annotate(f"平均EF: P=0.01 -> {p001_ef:.1f}  |  P=0.05 -> {p005_ef:.1f}",
                (0.5, 0.95), xycoords="axes fraction", fontsize=8,
                color="#555555", ha="center", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.8))

    # Right: stacked bar for shared vs unique pathways
    ax2 = axes[1]
    shared = p001_counts[0]  # Shared = all P=0.01 pathways (P=0.05 is a strict superset)
    unique_005 = p005_counts[0] - p001_counts[0]

    ax2.barh(0, shared, left=0, color="#4393C3", edgecolor="white", linewidth=0.3,
             label=f"共享(n={shared})")
    ax2.barh(0, unique_005, left=shared, color="#F4A582", edgecolor="white", linewidth=0.3,
             label=f"P=0.05独有(n={unique_005})")
    ax2.set_yticks([])
    ax2.set_xlabel("通路数", fontsize=10)
    ax2.legend(frameon=True, fancybox=False, fontsize=7, loc="upper right",
               edgecolor="#cccccc")
    ax2.grid(axis="x", alpha=0.3, linestyle="--", linewidth=0.4)

    fig.suptitle("P值阈值稳健性（P=0.01 vs P=0.05）", fontsize=12, y=0.98)
    save_fig(fig, "fig4_pvalue_comparison")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# FIGURE 5 — Robust Targets vs Parametric Combinations
# ══════════════════════════════════════════════════════════════════════


def fig5_robust_vs_param(df_grid: pd.DataFrame, comp_path: Path):
    """Bar chart comparing Q1 / Tier1 / All-72-combo metrics."""
    with open(comp_path, "r", encoding="utf-8") as f:
        comp = json.load(f)

    combo = comp["combo_stats"]
    q1 = comp["q1_metrics"]
    tier1 = comp["tier1_metrics"]

    # Metrics to compare
    labels = ["显著通路数", "信号主干", "平均EF"]
    # Convert EF to comparable scale (divide by 10 for chart)
    ef_div = 10

    q1_vals = [q1["total_significant"], q1["backbone"], round(q1["avg_ef"] / ef_div, 1)]
    tier1_vals = [tier1["total_significant"], tier1["backbone"], round(tier1["avg_ef"] / ef_div, 1)]
    combo_vals = [combo["significant_pathways"]["mean"], round(combo["backbone_count"]["mean"], 1),
                  round(combo["avg_ef"]["mean"] / ef_div, 1)]

    fig, ax = plt.subplots(figsize=(7, 5))

    x = np.arange(len(labels))
    width = 0.25

    bars1 = ax.bar(x - width, q1_vals, width, label=f"Q1候选稳健靶点(n={comp['n_q1_targets']})",
                   color="#2166AC", edgecolor="white", linewidth=0.3)
    bars2 = ax.bar(x, tier1_vals, width, label=f"Tier1核心靶点(n={comp['n_tier1_targets']})",
                   color="#B2182B", edgecolor="white", linewidth=0.3)
    bars3 = ax.bar(x + width, combo_vals, width,
                   label=f"72组合均值", color="#D3D3D3",
                   edgecolor="#999999", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("数值", fontsize=10)
    ax.legend(frameon=True, fancybox=False, fontsize=7.5, edgecolor="#cccccc",
              loc="upper left")
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.4)

    # Value labels
    for bars, color in [(bars1, "#2166AC"), (bars2, "#B2182B"), (bars3, "#666666")]:
        for bar in bars:
            h = bar.get_height()
            val = int(h) if h == int(h) else h
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.3, str(val),
                    ha="center", va="bottom", fontsize=7, color=color, fontweight="bold")

    # Annotation: EF is scaled
    ax.annotate("EF值除以10用于显示", (2, max(q1_vals[2], tier1_vals[2], combo_vals[2])),
                fontsize=6.5, color="#666666", ha="center", va="bottom",
                xytext=(2, max(q1_vals[2], tier1_vals[2], combo_vals[2]) + 2),
                arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.3))

    fig.suptitle("稳健靶点子集 vs 全部参数组合", fontsize=12, y=0.98)
    save_fig(fig, "fig5_robust_vs_param_comparison")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# FIGURE 6 — Prob effect on PPI network size (supplementary)
# ══════════════════════════════════════════════════════════════════════


def fig6_prob_ppi_effect(df: pd.DataFrame):
    """Line plot: Prob threshold vs n_ppi_nodes/n_ppi_edges, grouped by STRING."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.subplots_adjust(top=0.88)

    for ax, metric, ylabel in zip(
        axes, ["n_ppi_nodes", "n_ppi_edges"], ["PPI网络节点数", "PPI网络边数"]
    ):
        for probe_val in sorted(df["prob_threshold"].unique()):
            sub = df[df["prob_threshold"] == probe_val].groupby("string_threshold")[metric].mean()
            ax.plot(
                sub.index.astype(str), sub.values, marker="o", label=f"Prob={probe_val}",
                color=COLOR_PROB.get(probe_val, "#999999"), linewidth=1.5, markersize=5,
            )

        ax.set_xlabel("STRING阈值", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.legend(frameon=True, fancybox=False, fontsize=6.5, edgecolor="#cccccc")
        ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.4)

    fig.suptitle("概率阈值对PPI网络规模的影响", fontsize=12, y=0.98)
    save_fig(fig, "fig6_prob_ppi_network_effect")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# FIGURE 7 — Target frequency distribution histogram
# ══════════════════════════════════════════════════════════════════════


def fig7_target_frequency_histogram(df: pd.DataFrame):
    """Stacked histogram of target appearance frequencies by robustness quadrant.

    Shows how targets distribute across frequency bins, color-coded by
    their Q1/Q2/Q3/Q4 classification.  Overlaid KDE shows the continuous
    distribution shape.  Key thresholds (median, Tier1) are annotated.
    """
    from scipy.stats import gaussian_kde

    fig, ax = plt.subplots(figsize=(9, 5.5))

    freqs = df["frequency"].values
    freq_median = np.median(freqs)

    # Extract quadrant prefix (Q1..Q4) – safe for CJK encoding issues
    df["_q"] = df["quadrant"].str.extract(r"(Q\d)")

    bins = np.linspace(0, 1, 11)          # 10 bins → wider, clearer bars
    x = (bins[:-1] + bins[1:]) / 2
    width = 0.088

    # Stacked bars by quadrant (bottom → top)
    prefixes = ["Q4", "Q3", "Q2", "Q1"]
    labels = ["Q4: 低频低连接", "Q3: 高频低连接",
              "Q2: 高连接低稳健", "Q1: 候选稳健"]
    colors = ["#D6604D", "#F4A582", "#92C5DE", "#2166AC"]

    bottom = np.zeros(len(x))
    for prefix, color, label in zip(prefixes, colors, labels):
        subset = df[df["_q"] == prefix]
        counts = (np.histogram(subset["frequency"], bins=bins)[0]
                  if len(subset) else np.zeros(len(x)))
        ax.bar(x, counts, width, bottom=bottom, label=label, color=color,
               edgecolor="white", linewidth=0.6, alpha=0.9, zorder=2)
        bottom += counts

    # KDE on twin axis
    kde = gaussian_kde(freqs)
    x_kde = np.linspace(0, 1, 300)
    bw = bins[1] - bins[0]
    y_kde = kde(x_kde) * len(freqs) * bw

    ax2 = ax.twinx()
    ax2.fill_between(x_kde, y_kde, alpha=0.10, color="#2166AC", zorder=1)
    ax2.plot(x_kde, y_kde, color="#2166AC", linewidth=1.5, alpha=0.5, zorder=1)
    ax2.set_ylabel("核密度估计", fontsize=9, color="#2166AC")
    ax2.tick_params(axis="y", labelcolor="#2166AC", labelsize=8)
    ax2.set_ylim(0, max(y_kde) * 1.6)

    # Median threshold line
    ax.axvline(x=freq_median, color="#B2182B", linestyle="--",
               linewidth=1.5, alpha=0.8, zorder=3)
    y_top = ax.get_ylim()[1]
    ax.text(freq_median + 0.015, y_top * 0.92,
            "中位数 = {:.3f}".format(freq_median), fontsize=9,
            color="#B2182B", fontweight="bold", zorder=4)

    # Tier1 (freq = 1.0)
    tier1 = df[df["frequency"] == 1.0]["gene"].tolist()
    if tier1:
        ax.axvline(x=1.0, color="#DAA520", linestyle="-",
                   linewidth=2, alpha=0.8, zorder=3)
        ax.text(0.97, y_top * 0.75,
                "Tier1: " + ", ".join(tier1),
                fontsize=9, color="#B8860B", fontweight="bold",
                ha="right", zorder=4)

    # Summary statistics box
    n_q1 = int((df["_q"] == "Q1").sum())
    n_q2 = int((df["_q"] == "Q2").sum())
    n_q3 = int((df["_q"] == "Q3").sum())
    n_q4 = int((df["_q"] == "Q4").sum())
    stats = ("总计 {} 个靶点   |   Q1: {}   Q2: {}   Q3: {}   Q4: {}"
             .format(len(df), n_q1, n_q2, n_q3, n_q4))
    ax.text(0.5, 1.12, stats, transform=ax.transAxes, fontsize=9,
            ha="center", va="bottom", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8f8f8",
                      edgecolor="#bbbbbb", alpha=0.95))

    ax.set_xlabel("跨参数组合出现频率", fontsize=10)
    ax.set_ylabel("靶点数量", fontsize=10)
    ax.set_xlim(-0.05, 1.05)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9,
              edgecolor="#cccccc", title="稳健性象限", title_fontsize=8)
    ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.4, zorder=0)

    fig.suptitle("靶点出现频率分布（按稳健性象限分层）", fontsize=12, y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 1.05])
    save_fig(fig, "fig7_target_frequency_distribution")
    plt.close(fig)

    # clean temporary column
    df.drop(columns=["_q"], inplace=True, errors="ignore")


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════


def main():
    print("Loading data...")
    df_grid = pd.read_csv(DATA["grid"])
    df_robust = pd.read_csv(DATA["robustness"])
    df_pathway = pd.read_csv(DATA["pathway_stability"])

    print(f"  Grid results: {len(df_grid)} combos")
    print(f"  Target robustness: {len(df_robust)} genes")
    print(f"  Pathway stability: {len(df_pathway)} pathways")

    print("\nGenerating Figure 1: Heatmap — STRING × TopK faceted by Prob...")
    fig1_heatmap(df_grid)

    print("\nGenerating Figure 2: Target robustness quadrant scatter...")
    fig2_target_robustness(df_robust)

    print("\nGenerating Figure 3: Pathway rank stability (top 25)...")
    fig3_pathway_rank_stability(df_pathway)

    print("\nGenerating Figure 4: Pvalue threshold comparison...")
    fig4_pvalue_comparison()

    print("\nGenerating Figure 5: Robust targets vs parametric combinations...")
    fig5_robust_vs_param(df_grid, DATA["comparison"])

    print("\nGenerating Figure 6: Prob threshold effect on PPI network...")
    fig6_prob_ppi_effect(df_grid)

    print("\nGenerating Figure 7: Target frequency distribution...")
    fig7_target_frequency_histogram(df_robust)

    print(f"\n{'='*60}")
    print(f"All figures saved to: {FIGURES_DIR}")
    print(f"{'='*60}")
    print("\nGenerated files:")
    for f in sorted(FIGURES_DIR.glob("*")):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
