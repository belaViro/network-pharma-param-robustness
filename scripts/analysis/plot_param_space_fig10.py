"""
Fig.10: Parameter space bubble matrix.

STRING (x) × TopK (y) with bubble size = significant GO BP/KEGG functional terms,
color = Prob threshold. Shows stability regions in parameter space.
"""

import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=UserWarning)

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
FIGURES_DIR = OUTPUTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────
dpi = 300
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["SimHei", "Microsoft YaHei", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
    }
)


def save_fig(fig, name):
    for fmt in ("png", "svg"):
        path = FIGURES_DIR / f"{name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.12, format=fmt)
        print(f"  Saved: {path}")


# ═══════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════
df = pd.read_csv(OUTPUTS_DIR / "param_grid_results.csv")

# Map Prob to color
prob_colors = {0.0: "#2166AC", 0.09: "#4393C3", 0.1: "#F4A582"}
prob_labels = {0.0: "P=0.0", 0.09: "P=0.09", 0.1: "P=0.1"}
df["color"] = df["prob_threshold"].map(prob_colors)

# Jitter x to separate Prob bubbles at same STRING
np.random.seed(42)
string_vals = sorted(df["string_threshold"].unique())
string_pos = {s: i for i, s in enumerate(string_vals)}
df["x_pos"] = df["string_threshold"].map(string_pos)

# Add small horizontal offset per Prob within same STRING
prob_offsets = {0.0: -0.18, 0.09: 0.0, 0.1: 0.18}
df["x_jitter"] = df["x_pos"] + df["prob_threshold"].map(prob_offsets)

# Bubble size: sqrt scale for area
min_sig, max_sig = df["significant_pathways"].min(), df["significant_pathways"].max()
df["bubble_size"] = ((df["significant_pathways"] - min_sig) / (max_sig - min_sig + 1)) * 500 + 80

# ═══════════════════════════════════════════════════════════════════════
# Draw
# ═══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 5.5))

# Background grid
for i in range(len(string_vals)):
    ax.axvline(i, color="#eeeeee", linewidth=0.8, zorder=0)
for k in sorted(df["top_k"].unique()):
    ax.axhline(k, color="#eeeeee", linewidth=0.8, zorder=0)

# Scatter by Prob
for prob in sorted(df["prob_threshold"].unique()):
    sub = df[df["prob_threshold"] == prob]
    ax.scatter(sub["x_jitter"], sub["top_k"], s=sub["bubble_size"],
               c=prob_colors[prob], alpha=0.75, edgecolors="white",
               linewidths=0.8, zorder=3, label=prob_labels[prob])

# Add value labels inside bubbles for extreme values
for _, row in df.iterrows():
    if row["significant_pathways"] >= 70 or row["significant_pathways"] <= 20:
        ax.text(row["x_jitter"], row["top_k"], str(row["significant_pathways"]),
                ha="center", va="center", fontsize=6.5, color="white",
                fontweight="bold", zorder=4)

ax.set_xticks(range(len(string_vals)))
ax.set_xticklabels([f"{s}" for s in string_vals])
ax.set_xlabel("STRING 置信度阈值", fontsize=10)
ax.set_ylabel("Top-K", fontsize=10)
ax.set_yticks(sorted(df["top_k"].unique()))
ax.set_ylim(22, 43)
ax.set_xlim(-0.6, len(string_vals) - 0.4)

# Legend
legend = ax.legend(loc="upper right", fontsize=8, framealpha=0.9,
                   edgecolor="#cccccc", title="Prob 阈值", title_fontsize=8)

# Title
ax.set_title("参数空间稳定性矩阵（气泡大小 = 显著功能条目数）", fontsize=12, pad=10)

# Stats annotation
stats = []
for prob in sorted(df["prob_threshold"].unique()):
    sub = df[df["prob_threshold"] == prob]
    stats.append(f"{prob_labels[prob]}: {sub['significant_pathways'].mean():.1f}±{sub['significant_pathways'].std():.1f}")
stats_text = "  |  ".join(stats)
ax.text(0.5, -0.14, stats_text, transform=ax.transAxes, fontsize=8,
        ha="center", va="top", color="#666666")

ax.grid(False)
fig.tight_layout()

save_fig(fig, "fig10_parameter_space_bubble_matrix")
plt.close(fig)
print("Done!")
