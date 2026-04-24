"""
PPI Network Visualization
Generates PPI network graphs with hub-centric circular layout using matplotlib + networkx
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
from pathlib import Path
import logging


logger = logging.getLogger(__name__)

def plot_ppi_network(
    interactions_file: str = None,
    degrees_file: str = None,
    output_path: str = None,
    top_n: int = 30,
    figsize: tuple = (20, 20),
) -> str:
    """Plot PPI network with hub-centric two-layer circular layout.

    Hub nodes are placed in the inner circle (larger, red),
    regular nodes are placed in the outer circle (smaller, blue).

    Args:
        interactions_file: Path to string_interactions_interactions.csv
        degrees_file: Path to node_degrees.csv
        output_path: Output path for the figure (SVG format)
        top_n: Number of top hub genes for inner circle (default 30)
        figsize: Figure size (default 20x20)

    Returns:
        Path to saved SVG figure
    """
    # Default paths - relative to src/ directory
    src_dir = Path(__file__).parent.parent

    if interactions_file is None:
        interactions_file = src_dir / "outputs" / "shared" / "ppi_network" / "string_interactions_interactions.csv"
    else:
        interactions_file = Path(interactions_file)

    if degrees_file is None:
        degrees_file = src_dir / "outputs" / "shared" / "ppi_network" / "node_degrees.csv"
    else:
        degrees_file = Path(degrees_file)

    if output_path is None:
        output_path = src_dir / "outputs" / "shared" / "figures" / "ppi_network.svg"
    else:
        output_path = Path(output_path)

    # Read data
    interactions = pd.read_csv(interactions_file)
    degrees = pd.read_csv(degrees_file)

    # Create graph
    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

    # Add degree attribute
    degree_dict = dict(zip(degrees["Gene"], degrees["Degree"]))
    nx.set_node_attributes(G, degree_dict, "degree")

    # Get top hub genes
    top_hubs = degrees.nlargest(top_n, "Degree")["Gene"].tolist()
    hub_nodes = set(top_hubs) & set(G.nodes())
    regular_nodes = set(G.nodes()) - hub_nodes

    # Save hub nodes table
    hub_df = degrees[degrees["Gene"].isin(hub_nodes)].copy()
    hub_df = hub_df.sort_values("Degree", ascending=False).reset_index(drop=True)
    hub_df.index = hub_df.index + 1  # 1-based ranking
    hub_df.index.name = "Rank"
    hub_table_path = src_dir / "outputs" / "shared" / "ppi_network" / "hub_nodes.csv"
    hub_table_path.parent.mkdir(parents=True, exist_ok=True)
    hub_df.to_csv(hub_table_path)

    # Create hub-centric two-layer circular layout
    pos = {}

    # Inner circle for hub nodes (extended halfway toward outer)
    n_hubs = len(hub_nodes)
    if n_hubs > 0:
        hub_radius = 2.8
        hub_list = sorted(hub_nodes)
        for i, node in enumerate(hub_list):
            angle = 2 * np.pi * i / n_hubs - np.pi / 2
            pos[node] = (hub_radius * np.cos(angle), hub_radius * np.sin(angle))

    # Outer circle for regular nodes
    n_regular = len(regular_nodes)
    if n_regular > 0:
        outer_radius = 5.2
        regular_list = sorted(regular_nodes)
        for i, node in enumerate(regular_list):
            angle = 2 * np.pi * i / n_regular - np.pi / 2
            pos[node] = (outer_radius * np.cos(angle), outer_radius * np.sin(angle))

    # Create figure with white background
    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    ax.set_facecolor("white")

    # Draw edges with curved style
    for u, v, d in G.edges(data=True):
        x = [pos[u][0], pos[v][0]]
        y = [pos[u][1], pos[v][1]]
        score = d["weight"]
        if score >= 0.7:
            color = "#27ae60"
            alpha = 0.6
        elif score >= 0.5:
            color = "#f39c12"
            alpha = 0.5
        else:
            color = "#95a5a6"
            alpha = 0.3
        ax.plot(x, y, color=color, alpha=alpha, linewidth=0.8, zorder=1)

    # Draw regular nodes (outer circle) - blue
    for node in regular_nodes:
        x, y = pos[node]
        degree = G.degree(node)
        size = 80 + degree * 5
        ax.scatter(x, y, s=size, c="#3498db", edgecolors="#2980b9",
                    linewidths=0.5, alpha=0.8, zorder=2)
        ax.annotate(node, (x, y), fontsize=5, ha="center", va="center",
                    fontweight="bold", color="#2c3e50", zorder=3)

    # Draw hub nodes (inner circle) - red
    for node in hub_nodes:
        x, y = pos[node]
        degree = G.degree(node)
        size = 300 + degree * 15
        ax.scatter(x, y, s=size, c="#e74c3c", edgecolors="#c0392b",
                    linewidths=1.0, alpha=0.9, zorder=2)
        ax.annotate(node, (x, y), fontsize=6, ha="center", va="center",
                    fontweight="bold", color="white", zorder=3)

    # Draw concentric circles as guide
    circle_outer = mpatches.Circle((0, 0), 5.4, fill=False, edgecolor="#bdc3c7",
                                     linewidth=1.0, linestyle="--", alpha=0.5)
    circle_inner = mpatches.Circle((0, 0), 3.0, fill=False, edgecolor="#e74c3c",
                                     linewidth=1.0, linestyle="--", alpha=0.4)
    ax.add_patch(circle_outer)
    ax.add_patch(circle_inner)

    # Title
    ax.set_title(
        f"PPI Network (Stage 6)\nHub Nodes: {n_hubs} | Regular Nodes: {n_regular} | Edges: {G.number_of_edges()}",
        fontsize=14,
        fontweight="bold",
        pad=20,
        color="#2c3e50",
    )

    # Legend
    legend_elements = [
        mlines.Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="#e74c3c",
            markersize=14, markeredgecolor="#c0392b", markeredgewidth=1,
            label=f"Hub Genes (Top {top_n})",
        ),
        mlines.Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="#3498db",
            markersize=10, markeredgecolor="#2980b9", markeredgewidth=0.5,
            label=f"Regular Nodes ({n_regular})",
        ),
    ]
    # Node type legend (upper left)
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9,
              framealpha=0.9, edgecolor="#bdc3c7")

    # Edge color legend (lower left) - add as second legend
    edge_legend = [
        mlines.Line2D([0], [0], color="#27ae60", linewidth=2, alpha=0.6, label="Score >= 0.7"),
        mlines.Line2D([0], [0], color="#f39c12", linewidth=2, alpha=0.5, label="Score >= 0.5"),
        mlines.Line2D([0], [0], color="#95a5a6", linewidth=2, alpha=0.3, label="Score < 0.5"),
    ]
    leg2 = ax.legend(handles=edge_legend, loc="lower left", fontsize=8,
                     framealpha=0.9, edgecolor="#bdc3c2", title="Combined Score")
    ax.add_artist(leg2)  # Prevent second legend from replacing the first

    ax.axis("off")
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_aspect("equal")

    plt.tight_layout()

    # Save as SVG and PNG
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, format="svg", bbox_inches="tight", facecolor="white")

    # Also save as PNG for preview (before close)
    png_path = output_path.with_suffix(".png")
    plt.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()

    logger.info(f"Hub nodes table saved: {hub_table_path}")

    return str(output_path)


if __name__ == "__main__":
    logger.info("Generating PPI network with two-layer layout...")
    ppi_path = plot_ppi_network()
    logger.info(f"PPI network saved: {ppi_path}")