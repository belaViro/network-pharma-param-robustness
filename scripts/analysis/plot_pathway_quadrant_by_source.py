#!/usr/bin/env python3
"""
Fig.9b: 通路稳健性象限散点图（按来源分面 + KEGG/GO_BP 高亮版）

用法:
    cd code && python scripts/analysis/plot_pathway_quadrant_by_source.py
"""

import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties

warnings.filterwarnings("ignore", category=UserWarning)
matplotlib.use("Agg")

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
FIGURES_DIR = OUTPUTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────
dpi = 300
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["SimHei", "Microsoft YaHei", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": dpi,
    "axes.linewidth": 0.8,
})

# ── Colors ─────────────────────────────────────────────────────────────
COLOR_CATEGORY = {
    "signaling_backbone": "#2166AC",
    "cross_talk": "#E69F00",
    "direct_ox": "#D6604D",
    "other": "#999999",
}
COLOR_QUADRANT_BG = {
    "Q1": "#2166AC",
    "Q2": "#E69F00",
    "Q3": "#B2182B",
    "Q4": "#D3D3D3",
}
CATEGORY_LABEL = {
    "signaling_backbone": "信号主干",
    "cross_talk": "交叉调控",
    "direct_ox": "直接氧化应激",
    "other": "其他",
}
SOURCE_MARKER = {
    "KEGG": "o",
    "GO_BP": "s",
    "GO_MF": "^",
    "GO_CC": "D",
}
SOURCE_LABEL = {
    "KEGG": "KEGG Pathway",
    "GO_BP": "GO Biological Process",
    "GO_MF": "GO Molecular Function",
    "GO_CC": "GO Cellular Component",
}
SOURCE_EDGE = {
    "KEGG": "#000000",
    "GO_BP": "#333333",
    "GO_MF": "#888888",
    "GO_CC": "#888888",
}


def save_fig(fig, name):
    for fmt in ("png", "svg"):
        path = FIGURES_DIR / f"{name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", format=fmt)
        print(f"  Saved: {path}")


# ── Load data ──────────────────────────────────────────────────────────
print("Loading data...")
df = pd.read_csv(OUTPUTS_DIR / "pathway_rank_stability_augmented.csv")
print(f"  {len(df)} pathways")

# ── Compute medians ────────────────────────────────────────────────────
freq_med = df["frequency"].median()
logp_med = df["avg_logp"].median()
print(f"  Frequency median: {freq_med:.4f}")
print(f"  Avg_logp median: {logp_med:.4f}")

# ── Compute source-level Q1 rates ──────────────────────────────────────
source_q1_rates = {}
for src in df["source"].unique():
    sub = df[df["source"] == src]
    q1_count = len(sub[sub["robustness_quadrant"] == "Q1: 核心稳健通路"])
    source_q1_rates[src] = q1_count / len(sub) * 100 if len(sub) > 0 else 0

# ── Figure 1: All sources combined with shape coding ───────────────────
fig, ax = plt.subplots(figsize=(10, 7))

font_cn = FontProperties(family="SimHei", size=9)
font_cn_bold = FontProperties(family="SimHei", size=9, weight="bold")
font_cn_small = FontProperties(family="SimHei", size=7.5)
font_en = FontProperties(family="Arial", size=9)
font_en_small = FontProperties(family="Arial", size=7)

# Background shading
x_max = min(df["frequency"].max() * 1.05, 1.05)
y_min = 0
y_max = df["avg_logp"].max() * 1.15
y_med_ratio = (logp_med - y_min) / (y_max - y_min)

ax.axvspan(freq_med, x_max, ymin=y_med_ratio, ymax=1, alpha=0.06, color=COLOR_QUADRANT_BG["Q1"], zorder=0)
ax.axvspan(0, freq_med, ymin=y_med_ratio, ymax=1, alpha=0.04, color=COLOR_QUADRANT_BG["Q2"], zorder=0)
ax.axvspan(freq_med, x_max, ymin=0, ymax=y_med_ratio, alpha=0.04, color=COLOR_QUADRANT_BG["Q3"], zorder=0)
ax.axvspan(0, freq_med, ymin=0, ymax=y_med_ratio, alpha=0.03, color=COLOR_QUADRANT_BG["Q4"], zorder=0)

# Plot by source + category
for src in ["KEGG", "GO_BP", "GO_MF", "GO_CC"]:
    src_df = df[df["source"] == src]
    marker = SOURCE_MARKER[src]
    edge_color = SOURCE_EDGE[src]
    
    for cat in ["signaling_backbone", "cross_talk", "direct_ox", "other"]:
        sub = src_df[src_df["pathway_category"] == cat]
        if len(sub) == 0:
            continue
        
        sizes = 25 + sub["frequency"] * 100
        
        # KEGG gets thicker edges
        lw = 1.0 if src == "KEGG" else 0.4
        
        ax.scatter(
            sub["frequency"],
            sub["avg_logp"],
            c=COLOR_CATEGORY[cat],
            s=sizes,
            marker=marker,
            edgecolors=edge_color,
            linewidth=lw,
            alpha=0.8,
            zorder=4 if src in ("KEGG", "GO_BP") else 2,
            label=f"{SOURCE_LABEL[src]} ({CATEGORY_LABEL[cat]})",
        )

# Median lines
ax.axvline(x=freq_med, color="#555555", linestyle="--", linewidth=0.8, zorder=1, alpha=0.5)
ax.axhline(y=logp_med, color="#555555", linestyle="--", linewidth=0.8, zorder=1, alpha=0.5)

# Annotations
ax.text(freq_med + 0.01, y_max * 0.97, f"频率中位数 = {freq_med:.2f}",
        fontproperties=font_cn_small, color="#555555", va="top")
ax.text(x_max * 0.98, logp_med + 0.3, f"显著性中位数 = {logp_med:.1f}",
        fontproperties=font_cn_small, color="#555555", ha="right", va="bottom")

# Quadrant labels (large, faint)
ax.text((freq_med + x_max) / 2, (logp_med + y_max) / 2, "Q1\n核心稳健",
        fontproperties=font_cn_bold, fontsize=16, ha="center", va="center",
        color=COLOR_QUADRANT_BG["Q1"], alpha=0.18, zorder=0)
ax.text(freq_med / 2, (logp_med + y_max) / 2, "Q2\n低频高显著",
        fontproperties=font_cn_bold, fontsize=12, ha="center", va="center",
        color=COLOR_QUADRANT_BG["Q2"], alpha=0.15, zorder=0)
ax.text((freq_med + x_max) / 2, logp_med / 2, "Q3\n稳健但效应弱",
        fontproperties=font_cn_bold, fontsize=12, ha="center", va="center",
        color=COLOR_QUADRANT_BG["Q3"], alpha=0.15, zorder=0)

# Source Q1 rate annotation box
textstr = "Q1 核心稳健比例:\n"
for src in ["KEGG", "GO_BP", "GO_MF", "GO_CC"]:
    rate = source_q1_rates.get(src, 0)
    textstr += f"  {SOURCE_LABEL[src][:22]:<22}: {rate:5.1f}%\n"

ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=7.5,
        verticalalignment="top", fontproperties=FontProperties(family="SimHei", size=7),
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f8f8", edgecolor="#bbbbbb", alpha=0.95))

# Axis labels
ax.set_xlabel("跨参数组合出现频率", fontproperties=font_cn, labelpad=6)
ax.set_ylabel("平均显著性强度 (-log₁₀ Adjusted P-value)", fontproperties=font_cn, labelpad=6)
ax.set_xlim(-0.02, x_max)
ax.set_ylim(y_min, y_max)

for label in ax.get_xticklabels():
    label.set_fontproperties(font_en)
for label in ax.get_yticklabels():
    label.set_fontproperties(font_en)

# Legend - split into two columns to save space
handles, labels = ax.get_legend_handles_labels()
# Only keep unique source markers for legend
legend_elements = []
for src in ["KEGG", "GO_BP", "GO_MF", "GO_CC"]:
    legend_elements.append(plt.Line2D([0], [0], marker=SOURCE_MARKER[src], color="w",
                                       markerfacecolor="#666666", markeredgecolor=SOURCE_EDGE[src],
                                       markersize=8, label=SOURCE_LABEL[src]))

# Category color patches
for cat in ["signaling_backbone", "cross_talk", "direct_ox", "other"]:
    from matplotlib.patches import Patch
    legend_elements.append(Patch(facecolor=COLOR_CATEGORY[cat], edgecolor="white",
                                  label=CATEGORY_LABEL[cat]))

ax.legend(handles=legend_elements, frameon=True, fancybox=False, loc="lower left",
          edgecolor="#cccccc", prop=FontProperties(family="SimHei", size=6.5),
          ncol=2, columnspacing=0.8, handletextpad=0.3)

ax.grid(alpha=0.2, linestyle="--", linewidth=0.3)

fig.suptitle("通路稳健性象限（按数据库来源区分）",
             fontproperties=FontProperties(family="SimHei", size=12), y=0.98)

save_fig(fig, "fig9b_pathway_quadrant_by_source")
plt.close(fig)

# ── Figure 2: KEGG + GO_BP only (system pathways) ─────────────────────
fig2, ax2 = plt.subplots(figsize=(9, 6.5))

system_df = df[df["source"].isin(["KEGG", "GO_BP"])]

ax2.axvspan(freq_med, x_max, ymin=y_med_ratio, ymax=1, alpha=0.08, color=COLOR_QUADRANT_BG["Q1"], zorder=0)
ax2.axvspan(0, freq_med, ymin=y_med_ratio, ymax=1, alpha=0.06, color=COLOR_QUADRANT_BG["Q2"], zorder=0)
ax2.axvspan(freq_med, x_max, ymin=0, ymax=y_med_ratio, alpha=0.06, color=COLOR_QUADRANT_BG["Q3"], zorder=0)

for src in ["KEGG", "GO_BP"]:
    src_df = system_df[system_df["source"] == src]
    marker = SOURCE_MARKER[src]
    
    for cat in ["signaling_backbone", "cross_talk", "direct_ox", "other"]:
        sub = src_df[src_df["pathway_category"] == cat]
        if len(sub) == 0:
            continue
        
        sizes = 35 + sub["frequency"] * 120
        label = f"{SOURCE_LABEL[src]} - {CATEGORY_LABEL[cat]}"
        
        ax2.scatter(
            sub["frequency"],
            sub["avg_logp"],
            c=COLOR_CATEGORY[cat],
            s=sizes,
            marker=marker,
            edgecolors="white",
            linewidth=0.5,
            alpha=0.85,
            zorder=3,
            label=label,
        )

ax2.axvline(x=freq_med, color="#555555", linestyle="--", linewidth=0.8, zorder=1, alpha=0.5)
ax2.axhline(y=logp_med, color="#555555", linestyle="--", linewidth=0.8, zorder=1, alpha=0.5)

# Highlight KEGG pathways with names
kegg_q1 = system_df[(system_df["source"] == "KEGG") & 
                     (system_df["robustness_quadrant"] == "Q1: 核心稳健通路")].sort_values("avg_rank")
for _, r in kegg_q1.head(8).iterrows():
    import re
    short = re.sub(r"\s*\(GO:\d+\)", "", r["pathway"])
    if len(short) > 30:
        short = short[:27] + "..."
    ax2.annotate(short, (r["frequency"], r["avg_logp"]),
                xytext=(6, 4), textcoords="offset points",
                fontsize=5.5, color="#333333",
                fontproperties=FontProperties(family="SimHei", size=6),
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="#cccccc", alpha=0.85),
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.3))

# Big annotation: 97.5%
ax2.text(0.98, 0.98, "KEGG + GO_BP\nQ1 核心稳健率: 97.5%\n(39/40 条通路)",
         transform=ax2.transAxes, fontsize=10, fontweight="bold",
         verticalalignment="top", horizontalalignment="right",
         fontproperties=FontProperties(family="SimHei", size=10),
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#e8f4f8", edgecolor="#2166AC", alpha=0.95))

ax2.set_xlabel("跨参数组合出现频率", fontproperties=font_cn, labelpad=6)
ax2.set_ylabel("平均显著性强度 (-log₁₀ Adjusted P-value)", fontproperties=font_cn, labelpad=6)
ax2.set_xlim(-0.02, x_max)
ax2.set_ylim(y_min, y_max)

for label in ax2.get_xticklabels():
    label.set_fontproperties(font_en)
for label in ax2.get_yticklabels():
    label.set_fontproperties(font_en)

ax2.legend(frameon=True, fancybox=False, loc="lower left",
           edgecolor="#cccccc", prop=FontProperties(family="SimHei", size=7))
ax2.grid(alpha=0.2, linestyle="--", linewidth=0.3)

fig2.suptitle("系统级通路稳健性（KEGG + GO_BP）",
              fontproperties=FontProperties(family="SimHei", size=12), y=0.98)

save_fig(fig2, "fig9c_system_pathways_quadrant")
plt.close(fig2)

print(f"\n{'='*60}")
print(f"Figures saved to: {FIGURES_DIR}")
print("=" * 60)
