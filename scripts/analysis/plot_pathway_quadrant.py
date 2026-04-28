#!/usr/bin/env python3
"""
Fig.9: 通路稳健性象限散点图

X轴: 跨参数组合出现频率
Y轴: 平均显著性强度 (-log10 Adjusted P-value)
颜色: 五层分类类别 (signaling_backbone / cross_talk / direct_ox / other)
形状/背景: 二维象限 (Q1-Q4)

用法:
    cd code && python scripts/analysis/plot_pathway_quadrant.py
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
    "other": "#D3D3D3",
}
COLOR_QUADRANT_BG = {
    "Q1": "#2166AC",   # 高频高显著
    "Q2": "#E69F00",   # 低频高显著
    "Q3": "#B2182B",   # 高频低显著
    "Q4": "#D3D3D3",   # 低频低显著
}

CATEGORY_LABEL = {
    "signaling_backbone": "信号主干",
    "cross_talk": "交叉调控",
    "direct_ox": "直接氧化应激",
    "other": "其他",
}


def save_fig(fig, name):
    for fmt in ("png", "svg"):
        path = FIGURES_DIR / f"{name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.12, format=fmt)
        print(f"  Saved: {path}")


def clean_pathway_name(name: str, max_len: int = 38) -> str:
    import re

    name = re.sub(r"\s*\(GO:\d+\)", "", str(name))
    replacements = {
        "transmembrane receptor protein tyrosine kinase": "TM receptor tyrosine kinase",
        "positive regulation of": "+ regulation of",
        "regulation of": "reg. of",
        "phosphatidylinositol 3-kinase": "PI3K",
        "protein kinase B": "AKT",
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    if len(name) > max_len:
        name = name[: max_len - 3].rstrip() + "..."
    return name


# ── Load data ──────────────────────────────────────────────────────────
print("Loading data...")
df = pd.read_csv(OUTPUTS_DIR / "pathway_rank_stability_augmented.csv")
print(f"  {len(df)} pathways")

# ── Compute medians ────────────────────────────────────────────────────
freq_med = df["frequency"].median()
logp_med = df["avg_logp"].median()
print(f"  Frequency median: {freq_med:.4f}")
print(f"  Avg_logp median: {logp_med:.4f}")

# ── Figure ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8.4, 5.9), constrained_layout=True)

# Font properties
font_cn = FontProperties(family="SimHei", size=9)
font_cn_bold = FontProperties(family="SimHei", size=9, weight="bold")
font_cn_small = FontProperties(family="SimHei", size=7.5)
font_en = FontProperties(family="Arial", size=9)
font_en_small = FontProperties(family="Arial", size=7)

# Background shading for quadrants
x_max = min(df["frequency"].max() * 1.05, 1.05)
y_min = 0
y_max = df["avg_logp"].max() * 1.15
y_med_ratio = (logp_med - y_min) / (y_max - y_min)

# Q1 (高频高显著) - 右上
ax.axvspan(freq_med, x_max, ymin=y_med_ratio, ymax=1, alpha=0.08, color=COLOR_QUADRANT_BG["Q1"], zorder=0)
# Q2 (低频高显著) - 左上
ax.axvspan(0, freq_med, ymin=y_med_ratio, ymax=1, alpha=0.06, color=COLOR_QUADRANT_BG["Q2"], zorder=0)
# Q3 (高频低显著) - 右下
ax.axvspan(freq_med, x_max, ymin=0, ymax=y_med_ratio, alpha=0.06, color=COLOR_QUADRANT_BG["Q3"], zorder=0)
# Q4 (低频低显著) - 左下
ax.axvspan(0, freq_med, ymin=0, ymax=y_med_ratio, alpha=0.04, color=COLOR_QUADRANT_BG["Q4"], zorder=0)

# Plot by category
for cat in ["signaling_backbone", "cross_talk", "direct_ox", "other"]:
    sub = df[df["pathway_category"] == cat]
    if len(sub) == 0:
        continue
    
    # Size by frequency (larger = more robust)
    sizes = 30 + sub["frequency"] * 120
    
    ax.scatter(
        sub["frequency"],
        sub["avg_logp"],
        c=COLOR_CATEGORY[cat],
        s=sizes,
        label=f"{CATEGORY_LABEL[cat]} (n={len(sub)})",
        edgecolors="white",
        linewidth=0.4,
        alpha=0.8,
        zorder=3,
    )

# Median divider lines
ax.axvline(x=freq_med, color="#555555", linestyle="--", linewidth=0.8, zorder=1, alpha=0.5)
ax.axhline(y=logp_med, color="#555555", linestyle="--", linewidth=0.8, zorder=1, alpha=0.5)

# Median annotations
ax.text(freq_med + 0.01, y_max * 0.92,
        f"频率中位数 = {freq_med:.2f}",
        fontproperties=font_cn_small, color="#555555", va="top")
ax.text(x_max * 0.98, logp_med + 0.3,
        f"显著性中位数 = {logp_med:.1f}",
        fontproperties=font_cn_small, color="#555555", ha="right", va="bottom")

# Quadrant labels
ax.text((freq_med + x_max) / 2, (logp_med + y_max) / 2,
        "Q1\n候选稳健", fontproperties=font_cn_bold, fontsize=12,
        ha="center", va="center", color=COLOR_QUADRANT_BG["Q1"], alpha=0.16, zorder=0)
ax.text(freq_med / 2, (logp_med + y_max) / 2,
        "Q2\n低频高显著", fontproperties=font_cn_bold, fontsize=12,
        ha="center", va="center", color=COLOR_QUADRANT_BG["Q2"], alpha=0.14, zorder=0)
ax.text((freq_med + x_max) / 2, logp_med / 2,
        "Q3\n稳健但效应弱", fontproperties=font_cn_bold, fontsize=12,
        ha="center", va="center", color=COLOR_QUADRANT_BG["Q3"], alpha=0.14, zorder=0)

# Keep the scatter readable: list representative Q1 pathways in a compact note
# instead of attaching many overlapping labels to the same dense point cluster.
q1_df = df[df["robustness_quadrant"] == "Q1: 候选稳健通路"].sort_values("avg_rank")
label_items = [clean_pathway_name(p, max_len=34) for p in q1_df["pathway"].head(5)]
if label_items:
    note = "代表性候选通路:\n" + "\n".join(f"- {item}" for item in label_items)
    ax.text(
        0.985,
        0.965,
        note,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontproperties=FontProperties(family="SimHei", size=6.5),
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#cccccc", alpha=0.90),
    )

# Axis labels
ax.set_xlabel("跨参数组合出现频率", fontproperties=font_cn, labelpad=6)
ax.set_ylabel("平均显著性强度 [-log10(adjusted p-value)]", fontproperties=font_cn, labelpad=6)
ax.set_xlim(-0.02, x_max)
ax.set_ylim(y_min, y_max)

# Tick labels
for label in ax.get_xticklabels():
    label.set_fontproperties(font_en)
for label in ax.get_yticklabels():
    label.set_fontproperties(font_en)

# Legend
leg = ax.legend(frameon=True, fancybox=False, loc="lower right",
                edgecolor="#cccccc", prop=font_cn_small, title="五层分类", title_fontsize=7.5)
ax.grid(alpha=0.2, linestyle="--", linewidth=0.3)

# Title
fig.suptitle("通路稳健性：二维象限分类（频率 × 显著性强度）",
             fontproperties=FontProperties(family="SimHei", size=12), y=1.03)

save_fig(fig, "fig9_pathway_robustness_quadrant")
plt.close(fig)

print(f"\n{'='*60}")
print(f"Figure saved to: {FIGURES_DIR}")
print("=" * 60)
