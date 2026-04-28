"""
Fig.8: Upset plot - Target set intersections across Prob thresholds.

Shows how target sets change across the current main-grid Prob thresholds.
"""

import json
import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

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
with open(OUTPUTS_DIR / "combo_targets.json") as f:
    combo_targets = json.load(f)

PROBS = [0, 0.09, 0.1]
LABELS = ["P0", "P09", "P10"]
COLORS = ["#4393C3", "#92C5DE", "#F4A582"]

prob_unions = {}
for p in PROBS:
    sets = [set(v) for k, v in combo_targets.items() if k.startswith(f"P{p}_")]
    prob_unions[p] = set().union(*sets)

# Compute intersection patterns
patterns = {}
all_genes = set().union(*prob_unions.values())
for gene in all_genes:
    membership = tuple(p for p in PROBS if gene in prob_unions[p])
    if membership not in patterns:
        patterns[membership] = []
    patterns[membership].append(gene)

# Sort by size desc
sorted_patterns = sorted(patterns.items(), key=lambda x: -len(x[1]))
# Keep all non-empty patterns (including singletons)
# But group singletons visually if too many
sorted_patterns = [(m, g) for m, g in sorted_patterns]

n_patterns = len(sorted_patterns)
n_sets = len(PROBS)

# ═══════════════════════════════════════════════════════════════════════
# Upset plot layout
# ═══════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(10, 6))
from matplotlib.gridspec import GridSpec

gs = GridSpec(2, 1, height_ratios=[2.5, 1.2], hspace=0.08)
ax_bar = fig.add_subplot(gs[0])
ax_mat = fig.add_subplot(gs[1], sharex=ax_bar)

# ── Bar chart (intersection sizes) ────────────────────────────────────
x = np.arange(n_patterns)
bar_heights = [len(genes) for _, genes in sorted_patterns]
bar_colors = []
for membership, genes in sorted_patterns:
    if len(membership) == 1:
        idx = PROBS.index(membership[0])
        bar_colors.append(COLORS[idx])
    elif set(membership) == set(PROBS):
        bar_colors.append("#2166AC")
    elif 0 in membership and 0.09 in membership and len(membership) == 2:
        bar_colors.append("#4393C3")  # P0=P0.09 only
    else:
        bar_colors.append("#74A9CF")  # mixed

bars = ax_bar.bar(x, bar_heights, color=bar_colors, edgecolor="white", linewidth=0.8, width=0.7, zorder=2)

# Bar value labels
for i, (bar, h) in enumerate(zip(bars, bar_heights)):
    ax_bar.text(bar.get_x() + bar.get_width() / 2, h + 0.5, str(h),
                ha="center", va="bottom", fontsize=9, fontweight="bold", zorder=3)

ax_bar.set_ylabel("靶点数量", fontsize=10)
ax_bar.set_xticks([])
ax_bar.set_xlim(-0.6, n_patterns - 0.4)
ax_bar.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.4, zorder=0)
ax_bar.spines["bottom"].set_visible(False)

# ── Membership matrix ─────────────────────────────────────────────────
# Build membership matrix: rows=sets, cols=patterns
for j, (membership, genes) in enumerate(sorted_patterns):
    for i, p in enumerate(PROBS):
        if p in membership:
            ax_mat.plot(j, i, "o", color=COLORS[i], markersize=10, markeredgecolor="white", markeredgewidth=0.8, zorder=3)

# Vertical connecting lines for multi-set intersections
for j, (membership, genes) in enumerate(sorted_patterns):
    if len(membership) > 1:
        indices = [PROBS.index(p) for p in membership]
        ax_mat.plot([j, j], [min(indices), max(indices)], color="#888888",
                    linewidth=2, zorder=2, solid_capstyle="round")

ax_mat.set_yticks(range(n_sets))
ax_mat.set_yticklabels(LABELS, fontsize=10)
ax_mat.set_xticks(x)
# Pattern labels below matrix
pattern_labels = []
for membership, genes in sorted_patterns:
    if len(membership) == 1:
        pattern_labels.append(f"仅{LABELS[PROBS.index(membership[0])]}")
    elif set(membership) == set(PROBS):
        pattern_labels.append("三阈值")
    else:
        names = [LABELS[PROBS.index(p)] for p in membership]
        pattern_labels.append("/".join(names))

ax_mat.set_xticklabels(pattern_labels, fontsize=8, rotation=30, ha="right")
ax_mat.set_ylim(-0.5, n_sets - 0.5)
ax_mat.invert_yaxis()
ax_mat.spines["top"].set_visible(False)
ax_mat.spines["right"].set_visible(False)
ax_mat.tick_params(axis="x", length=0)

# ── Annotations ────────────────────────────────────────────────────────
# Highlight key findings with text boxes
annotations = []
for j, (membership, genes) in enumerate(sorted_patterns):
    if set(membership) == set(PROBS):
        annotations.append((j, len(genes), f"三阈值共有\n({', '.join(genes[:3])}...)"))

for j, h, text in annotations:
    ax_bar.text(j, h + 2.5, text, ha="center", va="bottom", fontsize=7.5,
                color="#555555", style="italic",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="#fff8dc", edgecolor="#DAA520", alpha=0.8))

# Legend for bar colors
legend_elements = [
    mpatches.Patch(facecolor="#2166AC", edgecolor="white", label="3阈值共有（主网格交集）"),
    mpatches.Patch(facecolor="#4393C3", edgecolor="white", label="P0/P09共有"),
    mpatches.Patch(facecolor="#74A9CF", edgecolor="white", label="混合交集"),
]
for color, label in zip(COLORS, LABELS):
    legend_elements.append(mpatches.Patch(facecolor=color, edgecolor="white", label=f"仅{label}"))

ax_bar.legend(handles=legend_elements, loc="upper right", fontsize=7.5,
              framealpha=0.9, edgecolor="#cccccc", ncol=2)

fig.suptitle("Prob 阈值靶点集合交集（Upset 图）", fontsize=12, y=0.98)

save_fig(fig, "fig8_upset_prob_intersections")
plt.close(fig)
print("Done!")
