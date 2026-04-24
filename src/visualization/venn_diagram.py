"""
Venn Diagram Visualization
Generates Venn diagrams for intersection analysis using matplotlib
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import logging


logger = logging.getLogger(__name__)

def plot_venn_diagram(
    targets_file: str = None,
    disease_file: str = None,
    intersect_file: str = None,
    output_path: str = None,
) -> str:
    """Generate Venn diagram for Stage 5 intersection analysis.

    Args:
        targets_file: Path to merged_targets.csv
        disease_file: Path to merged_disease_genes.csv
        intersect_file: Path to intersected_targets.csv
        output_path: Output path for the figure

    Returns:
        Path to saved figure
    """
    # Default paths - relative to src/ directory
    src_dir = Path(__file__).parent.parent

    if targets_file is None:
        targets_file = src_dir / "outputs" / "shared" / "merged_targets.csv"
    else:
        targets_file = Path(targets_file)

    if disease_file is None:
        disease_file = src_dir / "outputs" / "shared" / "disease_genes" / "merged_disease_genes.csv"
    else:
        disease_file = Path(disease_file)

    if intersect_file is None:
        intersect_file = src_dir / "outputs" / "shared" / "intersected_targets.csv"
    else:
        intersect_file = Path(intersect_file)

    if output_path is None:
        output_path = src_dir / "outputs" / "shared" / "figures" / "venn_diagram.png"
    else:
        output_path = Path(output_path)

    # Read data
    targets_df = pd.read_csv(targets_file)
    disease_df = pd.read_csv(disease_file)
    intersect_df = pd.read_csv(intersect_file)

    swiss_targets = set(targets_df["Target_Gene"].tolist())
    disease_genes = set(disease_df["Gene_Symbol"].tolist())
    intersect_genes = set(intersect_df["Target_Gene"].tolist())

    only_swiss = len(swiss_targets - disease_genes)
    only_disease = len(disease_genes - swiss_targets)
    common = len(intersect_genes)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw two equal-sized circles
    radius = 3
    center1 = (-1.5, 0)
    center2 = (1.5, 0)

    circle1 = mpatches.Circle(center1, radius, color="#3498db", alpha=0.4, edgecolor="#2980b9", linewidth=3)
    circle2 = mpatches.Circle(center2, radius, color="#e74c3c", alpha=0.4, edgecolor="#c0392b", linewidth=3)
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Labels
    ax.text(
        center1[0] - radius * 0.4,
        radius * 0.5,
        f"Only Predicted\nTargets\n\n{only_swiss}",
        fontsize=13,
        ha="center",
        va="center",
        fontweight="bold",
        color="#2c3e50",
    )

    ax.text(
        center2[0] + radius * 0.4,
        radius * 0.5,
        f"Only Oxidative\nStress Genes\n\n{only_disease}",
        fontsize=13,
        ha="center",
        va="center",
        fontweight="bold",
        color="#2c3e50",
    )

    # Intersection
    ax.text(0, 0, f"{common}", fontsize=28, ha="center", va="center", fontweight="bold", color="#27ae60")
    ax.text(0, -0.7, "Intersection", fontsize=12, ha="center", va="top", color="#2c3e50", style="italic")

    # Set labels
    ax.text(
        center1[0],
        center1[1] + radius + 0.4,
        f"SwissTargetPrediction\n(Predicted Targets)\n{len(swiss_targets)}",
        fontsize=11,
        ha="center",
        va="bottom",
        fontweight="bold",
        color="#2980b9",
    )
    ax.text(
        center2[0],
        center2[1] + radius + 0.4,
        f"Oxidative Stress\nRelated Genes\n{len(disease_genes)}",
        fontsize=11,
        ha="center",
        va="bottom",
        fontweight="bold",
        color="#c0392b",
    )

    ax.set_xlim(-6, 6)
    ax.set_ylim(-4.5, 5)

    plt.tight_layout()

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()

    return str(output_path)


if __name__ == "__main__":
    logger.info("Generating Venn diagram...")
    venn_path = plot_venn_diagram()
    logger.info(f"Venn diagram saved: {venn_path}")