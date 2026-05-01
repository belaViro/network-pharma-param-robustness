#!/usr/bin/env python3
"""
Fig.9: Pathway robustness quadrant plot.

X-axis: occurrence frequency across parameter combinations.
Y-axis: average enrichment strength, -log10(adjusted p-value).
Color/marker: pathway robustness quadrant.

Usage:
    cd code && python scripts/analysis/plot_pathway_quadrant.py
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle

warnings.filterwarnings("ignore", category=UserWarning)


BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
FIGURES_DIR = OUTPUTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

INPUT_CSV = OUTPUTS_DIR / "pathway_robustness.csv"

dpi = 300
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["SimHei", "Microsoft YaHei", "Arial", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.dpi": 150,
        "savefig.dpi": dpi,
        "axes.linewidth": 1.0,
        "axes.unicode_minus": False,
    }
)

QUADRANT_STYLE = {
    "Q1: 候选稳健通路": {
        "label": "Q1: 候选稳健通路",
        "color": "#3778B8",
        "marker": "o",
        "size": 58,
        "zorder": 4,
        "bg": "#DCEAF6",
        "alpha": 0.88,
    },
    "Q2: 低频高显著": {
        "label": "Q2: 低频高显著",
        "color": "#E69F00",
        "marker": "s",
        "size": 52,
        "zorder": 5,
        "bg": "#FFF3D8",
        "alpha": 0.90,
    },
    "Q3: 稳健但效应弱": {
        "label": "Q3: 稳健但效应弱",
        "color": "#B83B4A",
        "marker": "^",
        "size": 55,
        "zorder": 5,
        "bg": "#F5DDE2",
        "alpha": 0.90,
    },
    "Q4: 可忽略": {
        "label": "Q4: 可忽略",
        "color": "#D0D0D0",
        "marker": "o",
        "size": 36,
        "zorder": 3,
        "bg": "#F4F4F4",
        "alpha": 0.75,
    },
}

CATEGORY_LABELS = {
    "direct_ox": "直接OX",
    "cross_talk": "交叉调控",
    "signaling_backbone": "信号主干",
    "other": "其他",
}


def clean_pathway_name(name: str, max_len: int = 34) -> str:
    name = re.sub(r"\s*\(GO:\d+\)", "", str(name))
    replacements = {
        "transmembrane receptor protein tyrosine kinase": "TM receptor tyrosine kinase",
        "cellular response to cytokine stimulus": "cellular response to cytokine",
        "cytokine-mediated signaling pathway": "cytokine-mediated signaling",
        "AGE-RAGE signaling pathway in diabetic complications": "AGE-RAGE signaling",
        "Fluid shear stress and atherosclerosis": "Fluid shear stress",
        "Proteoglycans in cancer": "Proteoglycans in cancer",
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    if len(name) > max_len:
        name = name[: max_len - 3].rstrip() + "..."
    return name


def save_fig(fig: plt.Figure, name: str) -> None:
    for fmt in ("png", "svg"):
        path = FIGURES_DIR / f"{name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.16, format=fmt)
        print(f"Saved: {path}")


def main() -> None:
    print("Loading pathway robustness data...")
    df = pd.read_csv(INPUT_CSV)
    print(f"  {len(df)} pathways")

    freq_med = float(df["frequency"].median())
    logp_med = float(df["avg_logp"].median())
    print(f"  Frequency median: {freq_med:.4f}")
    print(f"  Avg -log10(adjusted p) median: {logp_med:.4f}")

    x_min = -0.03
    x_max = min(max(1.05, float(df["frequency"].max()) * 1.05), 1.08)
    y_min = -0.5
    y_max = float(df["avg_logp"].max()) * 1.12

    fig, ax = plt.subplots(figsize=(10.6, 7.0), constrained_layout=True)

    font_cn = FontProperties(family="SimHei", size=10)
    font_cn_bold = FontProperties(family="SimHei", size=11, weight="bold")
    font_cn_small = FontProperties(family="SimHei", size=8)

    # Draw full rectangular quadrant backgrounds so the shaded areas fill
    # the whole axis range, including the small negative margins.
    backgrounds = [
        (freq_med, logp_med, x_max - freq_med, y_max - logp_med, "Q1: 候选稳健通路"),
        (x_min, logp_med, freq_med - x_min, y_max - logp_med, "Q2: 低频高显著"),
        (freq_med, y_min, x_max - freq_med, logp_med - y_min, "Q3: 稳健但效应弱"),
        (x_min, y_min, freq_med - x_min, logp_med - y_min, "Q4: 可忽略"),
    ]
    for x, y, width, height, quadrant in backgrounds:
        ax.add_patch(
            Rectangle(
                (x, y),
                width,
                height,
                facecolor=QUADRANT_STYLE[quadrant]["bg"],
                edgecolor="none",
                alpha=0.58 if quadrant == "Q1: 候选稳健通路" else 0.44,
                zorder=0,
            )
        )

    # Q4 first, then Q1-Q3 above it.
    plot_order = [
        "Q4: 可忽略",
        "Q3: 稳健但效应弱",
        "Q2: 低频高显著",
        "Q1: 候选稳健通路",
    ]
    for quadrant in plot_order:
        sub = df[df["robustness_quadrant"] == quadrant]
        if sub.empty:
            continue
        style = QUADRANT_STYLE[quadrant]
        ax.scatter(
            sub["frequency"],
            sub["avg_logp"],
            s=style["size"],
            c=style["color"],
            marker=style["marker"],
            edgecolors="white",
            linewidth=0.45,
            alpha=style["alpha"],
            zorder=style["zorder"],
            label=f"{style['label']}  (n={len(sub)})",
        )

    ax.axvline(freq_med, color="#8A8A8A", linestyle="--", linewidth=1.0, zorder=2)
    ax.axhline(logp_med, color="#8A8A8A", linestyle="--", linewidth=1.0, zorder=2)

    ax.text(
        freq_med + 0.012,
        y_max * 0.94,
        f"频率中位数 = {freq_med:.3f}",
        fontproperties=font_cn_small,
        color="#555555",
        va="top",
    )
    ax.text(
        x_max * 0.98,
        logp_med + 0.35,
        f"显著性中位数 = {logp_med:.2f}",
        fontproperties=font_cn_small,
        color="#555555",
        ha="right",
        va="bottom",
    )

    ax.text(
        (freq_med + x_max) / 2,
        (logp_med + y_max) / 2,
        "Q1\n候选稳健",
        fontproperties=font_cn_bold,
        fontsize=16,
        ha="center",
        va="center",
        color=QUADRANT_STYLE["Q1: 候选稳健通路"]["color"],
        alpha=0.14,
        zorder=1,
    )
    ax.text(
        (x_min + freq_med) / 2,
        (logp_med + y_max) / 2,
        "Q2\n低频高显著",
        fontproperties=font_cn_bold,
        fontsize=14,
        ha="center",
        va="center",
        color=QUADRANT_STYLE["Q2: 低频高显著"]["color"],
        alpha=0.16,
        zorder=1,
    )
    ax.text(
        (freq_med + x_max) / 2,
        (y_min + logp_med) / 2,
        "Q3\n稳健但效应弱",
        fontproperties=font_cn_bold,
        fontsize=14,
        ha="center",
        va="center",
        color=QUADRANT_STYLE["Q3: 稳健但效应弱"]["color"],
        alpha=0.16,
        zorder=1,
    )
    ax.text(
        (x_min + freq_med) / 2,
        (y_min + logp_med) / 2,
        "Q4\n可忽略",
        fontproperties=font_cn_bold,
        fontsize=14,
        ha="center",
        va="center",
        color="#999999",
        alpha=0.22,
        zorder=1,
    )

    # Compact annotation: top Q1 pathways by candidate rank.
    q1_top = (
        df[df["robustness_quadrant"] == "Q1: 候选稳健通路"]
        .sort_values("candidate_rank")
        .head(8)
    )
    top_lines = []
    for _, row in q1_top.iterrows():
        cat = CATEGORY_LABELS.get(row["pathway_category"], row["pathway_category"])
        top_lines.append(
            f"{int(row['candidate_rank'])}. {clean_pathway_name(row['pathway'])} [{cat}]"
        )
    if top_lines:
        ax.text(
            0.985,
            0.965,
            "Q1代表通路:\n" + "\n".join(top_lines),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontproperties=FontProperties(family="SimHei", size=7.2),
            color="#333333",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="white",
                edgecolor="#CCCCCC",
                alpha=0.92,
            ),
            zorder=10,
        )

    ax.set_title("通路稳健性：二维象限分类（平均 -log10 adjusted p-value）", fontproperties=font_cn_bold, pad=12)
    ax.set_xlabel("跨参数组合出现频率", fontproperties=font_cn, labelpad=8)
    ax.set_ylabel("平均显著性强度 [-log10(adjusted p-value)]", fontproperties=font_cn, labelpad=8)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.grid(alpha=0.22, linestyle="--", linewidth=0.4)

    legend = ax.legend(
        loc="upper left",
        frameon=True,
        fancybox=False,
        edgecolor="#CCCCCC",
        prop=FontProperties(family="SimHei", size=8.5),
    )
    legend.get_frame().set_alpha(0.92)

    save_fig(fig, "fig9_pathway_robustness_quadrant")
    plt.close(fig)

    print(f"Figure saved to: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
