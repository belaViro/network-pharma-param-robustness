"""
Fig.9: Q1 core targets PPI subnetwork.

Shows the interaction network among 31 Q1 robust targets,
with node size = degree and color = frequency.
"""

import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=UserWarning)

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
SHARED_DIR = BASE_DIR / "src" / "outputs" / "shared"
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
target_df = pd.read_csv(OUTPUTS_DIR / "target_robustness.csv")
target_df["q"] = target_df["quadrant"].str.extract(r"(Q\d)")
q1_df = target_df[target_df["q"] == "Q1"].copy()
q1_genes = set(q1_df["gene"].str.upper().tolist())

# Build frequency/degree lookup
freq_map = dict(zip(q1_df["gene"].str.upper(), q1_df["frequency"]))
degree_map = dict(zip(q1_df["gene"].str.upper(), q1_df["avg_degree"]))

# Load STRING edges
edges_df = pd.read_csv(SHARED_DIR / "ppi_network" / "string_interactions_interactions.csv")
q1_edges = edges_df[(edges_df["protein1"].isin(q1_genes)) & (edges_df["protein2"].isin(q1_genes))]

# Build graph
G = nx.Graph()
for _, row in q1_edges.iterrows():
    G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

# Compute layout (seeded for reproducibility)
pos = nx.spring_layout(G, k=2.2, iterations=150, seed=42, weight="weight")

# Node attributes
node_freqs = [freq_map.get(n, 0) for n in G.nodes()]
node_degrees = [degree_map.get(n, 0) for n in G.nodes()]

# Size: map degree to marker size (scale for visibility)
min_deg, max_deg = min(node_degrees), max(node_degrees)
node_sizes = [150 + (d - min_deg) / (max_deg - min_deg) * 1200 for d in node_degrees]

# Color: frequency (blue → red)
cmap = LinearSegmentedColormap.from_list("freq", ["#92C5DE", "#2166AC", "#B2182B"])

# ═══════════════════════════════════════════════════════════════════════
# Draw
# ═══════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(7.5, 6.5))

# Draw edges with gradient color based on weight
# Low weight (0.4) -> very light gray, high weight (1.0) -> dark charcoal
min_w = 0.4
max_w = 1.0
for (u, v, d) in G.edges(data=True):
    w = d["weight"]
    # Interpolate gray: light (#E0E0E0) -> dark (#333333)
    ratio = (w - min_w) / (max_w - min_w)
    gray_val = int(224 - ratio * (224 - 51))  # 224->51
    edge_color = f"#{gray_val:02x}{gray_val:02x}{gray_val:02x}"
    alpha = 0.15 + ratio * 0.55  # 0.15 -> 0.70
    width = 0.3 + ratio * 0.6    # 0.3 -> 0.9
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v)], alpha=alpha,
                           width=width, edge_color=edge_color)

# Draw nodes
nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                                node_color=node_freqs, cmap=cmap,
                                edgecolors="white", linewidths=2.0)

# Labels: all nodes
labels = {n: n for n in G.nodes()}
nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=6.5,
                        font_weight="bold", font_color="#222222")

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min(node_freqs), vmax=max(node_freqs)))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=20, pad=0.02)
cbar.set_label("出现频率", fontsize=9)
cbar.ax.tick_params(labelsize=8)

# Title and stats
fig.suptitle("稳健靶点的PPI网络互作图", fontsize=13, fontweight="bold", y=0.98)

# Stats box
tier1_genes = [g for g, f in freq_map.items() if f == 1.0]
stats = (
    f"节点: {G.number_of_nodes()}   边: {G.number_of_edges()}\n"
    f"平均度: {2*G.number_of_edges()/G.number_of_nodes():.1f}\n"
    f"密度: {nx.density(G):.3f}\n"
    f"Tier1: {', '.join(tier1_genes)}"
)
ax.text(0.02, 0.98, stats, transform=ax.transAxes, fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f8f8f8",
                  edgecolor="#cccccc", alpha=0.9))

ax.set_axis_off()
fig.tight_layout()

save_fig(fig, "fig9_q1_ppi_subnetwork")
plt.close(fig)
print("Done!")
