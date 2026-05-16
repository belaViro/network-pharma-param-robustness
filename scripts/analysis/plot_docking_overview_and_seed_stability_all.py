from __future__ import annotations

import csv
import shutil
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties


ROOT = Path(__file__).resolve().parents[2]
DOCKING_DIR = ROOT / "分子对接" / "results"
VALIDATION_DIR = DOCKING_DIR / "validation_runs"
PAPER_DIR = ROOT / "src" / "outputs" / "paper_figures"

SUMMARY_CSV = VALIDATION_DIR / "docking_replicates_summary.csv"
REPLICATES_CSV = VALIDATION_DIR / "docking_replicates.csv"

RECEPTORS = ["STAT3", "EGFR", "SRC", "TNF", "BCL2", "HSP90AA1"]
LIGANDS = ["Catechin", "Miquelianin", "Procyanidin_B1"]
LIGAND_LABELS = {
    "Catechin": "Catechin",
    "Miquelianin": "Miquelianin",
    "Procyanidin_B1": "Procyanidin B1",
}

TARGET_COLORS = {
    "STAT3": "#5B677A",
    "EGFR": "#2563EB",
    "SRC": "#B91C1C",
    "TNF": "#D97706",
    "BCL2": "#047857",
    "HSP90AA1": "#7C3AED",
}

LIGAND_COLORS = {
    "Catechin": "#F59E0B",
    "Miquelianin": "#8B5CF6",
    "Procyanidin_B1": "#14B8A6",
}


def setup_style() -> tuple[FontProperties, FontProperties, FontProperties, FontProperties]:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
        }
    )
    font = FontProperties(family=["Arial", "Microsoft YaHei", "SimHei"], size=9.5)
    small = FontProperties(family=["Arial", "Microsoft YaHei", "SimHei"], size=8.0)
    tiny = FontProperties(family=["Arial", "Microsoft YaHei", "SimHei"], size=7.0)
    panel = FontProperties(family=["Arial", "Microsoft YaHei", "SimHei"], size=13, weight="bold")
    return font, small, tiny, panel


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def pair_label(row: dict[str, str]) -> str:
    return f"{row['receptor']} x {LIGAND_LABELS[row['ligand']]}"


def load_data() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if not SUMMARY_CSV.exists():
        raise FileNotFoundError(SUMMARY_CSV)
    if not REPLICATES_CSV.exists():
        raise FileNotFoundError(REPLICATES_CSV)
    summary = read_csv(SUMMARY_CSV)
    replicates = read_csv(REPLICATES_CSV)
    summary = [row for row in summary if row["receptor"] in RECEPTORS and row["ligand"] in LIGANDS]
    replicates = [row for row in replicates if row["receptor"] in RECEPTORS and row["ligand"] in LIGANDS]
    if len(summary) != len(RECEPTORS) * len(LIGANDS):
        raise ValueError(f"Expected 18 summary rows, found {len(summary)}")
    if len(replicates) != len(RECEPTORS) * len(LIGANDS) * 3:
        raise ValueError(f"Expected 54 replicate rows, found {len(replicates)}")
    return summary, replicates


def matrix_from_summary(summary: list[dict[str, str]], key: str) -> np.ndarray:
    lookup = {(row["receptor"], row["ligand"]): as_float(row, key) for row in summary}
    return np.array([[lookup[(receptor, ligand)] for ligand in LIGANDS] for receptor in RECEPTORS])


def panel_label(ax: plt.Axes, label: str, panel_font: FontProperties) -> None:
    ax.text(-0.08, 1.025, label, transform=ax.transAxes, ha="left", va="bottom", fontproperties=panel_font)


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    PAPER_DIR.mkdir(parents=True, exist_ok=True)
    DOCKING_DIR.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for ext in ("png", "svg", "pdf"):
        path = PAPER_DIR / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.08)
        paths.append(path)
        shutil.copy2(path, DOCKING_DIR / path.name)
    plt.close(fig)
    return paths


def draw_heatmap(
    ax: plt.Axes,
    summary: list[dict[str, str]],
    font: FontProperties,
    small: FontProperties,
    show_title: bool = True,
    add_colorbar: bool = True,
    cbar_ax: plt.Axes | None = None,
) -> None:
    means = matrix_from_summary(summary, "mean_affinity_kcal_mol")
    sds = matrix_from_summary(summary, "sd_affinity_kcal_mol")
    cmap = LinearSegmentedColormap.from_list(
        "affinity",
        ["#FFF7BC", "#FEC44F", "#F97316", "#B91C1C", "#581C1C"],
    )
    im = ax.imshow(means, cmap=cmap.reversed(), vmin=-10.0, vmax=-5.7, aspect="auto")
    ax.set_xticks(np.arange(len(LIGANDS)))
    ax.set_xticklabels([LIGAND_LABELS[x] for x in LIGANDS], fontproperties=small, rotation=18, ha="right")
    ax.set_yticks(np.arange(len(RECEPTORS)))
    ax.set_yticklabels(RECEPTORS, fontproperties=font)
    ax.set_xlabel("Ligand", fontproperties=font)
    ax.set_ylabel("Target", fontproperties=font)
    if show_title:
        ax.set_title("6 x 3 mean Vina score", fontproperties=font, pad=7)

    for i in range(means.shape[0]):
        for j in range(means.shape[1]):
            color = "white" if means[i, j] <= -8.3 else "#111827"
            ax.text(
                j,
                i,
                f"{means[i, j]:.2f}\n±{sds[i, j]:.2f}",
                ha="center",
                va="center",
                color=color,
                fontproperties=small,
                linespacing=1.05,
            )

    if add_colorbar:
        cbar = ax.figure.colorbar(im, ax=ax if cbar_ax is None else None, cax=cbar_ax, fraction=0.045, pad=0.028)
        cbar.ax.set_title("kcal/mol", fontproperties=small, pad=3)
        cbar.ax.tick_params(labelsize=7.2)
    ax.set_box_aspect(1)


def draw_ranked_lollipop(ax: plt.Axes, summary: list[dict[str, str]], font: FontProperties, small: FontProperties, show_title: bool = True) -> None:
    ranked = sorted(summary, key=lambda row: as_float(row, "mean_affinity_kcal_mol"))
    labels = [pair_label(row) for row in ranked]
    x = [as_float(row, "mean_affinity_kcal_mol") for row in ranked]
    err = [as_float(row, "sd_affinity_kcal_mol") for row in ranked]
    y = np.arange(len(ranked))

    ax.hlines(y, xmin=0, xmax=x, color="#CBD5E1", linewidth=1.0, zorder=1)
    for idx, row in enumerate(ranked):
        ax.errorbar(
            x[idx],
            y[idx],
            xerr=err[idx],
            fmt="o",
            color=TARGET_COLORS[row["receptor"]],
            ecolor="#111827",
            elinewidth=0.9,
            capsize=2.5,
            markersize=5.0,
            zorder=3,
        )
    ax.axvline(-8.0, color="#94A3B8", linewidth=0.8, linestyle="--")
    ax.axvline(-9.0, color="#64748B", linewidth=0.8, linestyle=":")
    ax.set_yticks(y)
    ax.set_yticklabels([fill(label, 26) for label in labels], fontproperties=small)
    ax.invert_yaxis()
    ax.set_xlim(-10.1, -5.45)
    ax.set_xlabel("Mean Vina score (kcal/mol)", fontproperties=font)
    if show_title:
        ax.set_title("All 18 receptor-ligand pairs ranked", fontproperties=font, pad=7)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.text(-8.02, -0.95, "-8", ha="center", va="center", fontproperties=small, color="#64748B")
    ax.text(-9.02, -0.95, "-9", ha="center", va="center", fontproperties=small, color="#64748B")


def draw_best_ligand_by_target(
    ax: plt.Axes,
    summary: list[dict[str, str]],
    font: FontProperties,
    small: FontProperties,
    show_title: bool = True,
) -> None:
    best_rows = []
    for receptor in RECEPTORS:
        rows = [row for row in summary if row["receptor"] == receptor]
        best_rows.append(min(rows, key=lambda row: as_float(row, "mean_affinity_kcal_mol")))
    best_rows = sorted(best_rows, key=lambda row: as_float(row, "mean_affinity_kcal_mol"))
    y = np.arange(len(best_rows))
    strengths = [-as_float(row, "mean_affinity_kcal_mol") for row in best_rows]
    sds = [as_float(row, "sd_affinity_kcal_mol") for row in best_rows]

    ax.barh(
        y,
        strengths,
        color=[LIGAND_COLORS[row["ligand"]] for row in best_rows],
        edgecolor="#111827",
        linewidth=0.6,
        alpha=0.86,
        height=0.58,
    )
    ax.errorbar(strengths, y, xerr=sds, fmt="none", ecolor="#111827", elinewidth=0.8, capsize=2.4, zorder=3)
    for yi, row, strength in zip(y, best_rows, strengths):
        score = as_float(row, "mean_affinity_kcal_mol")
        text_color = "white" if row["ligand"] != "Catechin" else "#111827"
        ax.text(
            min(strength - 0.12, 9.85),
            yi,
            f"{LIGAND_LABELS[row['ligand']]}  {score:.2f}",
            ha="right",
            va="center",
            color=text_color,
            fontproperties=small,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([row["receptor"] for row in best_rows], fontproperties=small)
    ax.invert_yaxis()
    ax.set_xlim(0, 10.4)
    ax.set_xlabel("Best docking strength (-Vina score, kcal/mol)", fontproperties=font)
    if show_title:
        ax.set_title("Best ligand for each target", fontproperties=font, pad=7)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.set_box_aspect(1)


def draw_ligand_distribution(ax: plt.Axes, summary: list[dict[str, str]], font: FontProperties, small: FontProperties, show_title: bool = True) -> None:
    data = [[as_float(row, "mean_affinity_kcal_mol") for row in summary if row["ligand"] == ligand] for ligand in LIGANDS]
    positions = np.arange(1, len(LIGANDS) + 1)
    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.52,
        patch_artist=True,
        medianprops={"color": "#111827", "linewidth": 1.2},
        boxprops={"edgecolor": "#334155", "linewidth": 0.9},
        whiskerprops={"color": "#334155", "linewidth": 0.9},
        capprops={"color": "#334155", "linewidth": 0.9},
    )
    for patch, ligand in zip(box["boxes"], LIGANDS):
        color = LIGAND_COLORS[ligand]
        patch.set_facecolor(color)
        patch.set_alpha(0.56)

    for pos, values, ligand in zip(positions, data, LIGANDS):
        jitter = np.linspace(-0.12, 0.12, len(values))
        ax.scatter(np.full(len(values), pos) + jitter, values, s=22, color="#FFFFFF", edgecolor="#1F2937", linewidth=0.6, zorder=3)
        ax.text(pos, min(values) - 0.13, f"best {min(values):.2f}", ha="center", va="top", fontproperties=small, color="#374151")

    ax.set_xticks(positions)
    ax.set_xticklabels([LIGAND_LABELS[x] for x in LIGANDS], fontproperties=small, rotation=14, ha="right")
    ax.set_ylabel("Mean Vina score (kcal/mol)", fontproperties=font)
    if show_title:
        ax.set_title("Ligand-level distribution", fontproperties=font, pad=7)
    ax.set_ylim(-10.15, -5.45)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
    ax.tick_params(axis="y", labelsize=7.5)
    ax.set_box_aspect(1)


def figure_overview(summary: list[dict[str, str]], replicates: list[dict[str, str]]) -> list[Path]:
    font, small, _, panel = setup_style()
    fig = plt.figure(figsize=(13.8, 12.2))
    left = 0.075
    left_width = 0.335
    panel_size = 0.255
    cbar_width = 0.010
    b_left = 0.585
    b_bottom = 0.075
    b_width = 0.375
    b_height = 0.875

    ax_a = fig.add_axes([left, 0.695, left_width, panel_size])
    ax_c = fig.add_axes([left, 0.385, left_width, panel_size])
    ax_d = fig.add_axes([left, 0.075, left_width, panel_size])
    ax_b = fig.add_axes([b_left, b_bottom, b_width, b_height])
    cbar_ax = fig.add_axes([left + left_width + 0.006, 0.725, cbar_width, 0.195])

    draw_heatmap(ax_a, summary, font, small, show_title=False, cbar_ax=cbar_ax)
    draw_seed_points(ax_b, summary, replicates, font, small, show_title=False, legend_loc="upper right")
    draw_best_ligand_by_target(ax_c, summary, font, small, show_title=False)
    draw_ligand_distribution(ax_d, summary, font, small, show_title=False)

    for label, ax in zip(["A", "B", "C", "D"], [ax_a, ax_b, ax_c, ax_d]):
        panel_label(ax, label, panel)
    return save_figure(fig, "fig_docking_overview_all")


def draw_seed_points(
    ax: plt.Axes,
    summary: list[dict[str, str]],
    replicates: list[dict[str, str]],
    font: FontProperties,
    small: FontProperties,
    show_title: bool = True,
    legend_loc: str = "lower right",
) -> None:
    ranked = sorted(summary, key=lambda row: as_float(row, "mean_affinity_kcal_mol"))
    labels = [pair_label(row) for row in ranked]
    y_lookup = {(row["receptor"], row["ligand"]): idx for idx, row in enumerate(ranked)}
    y = np.arange(len(ranked))

    for idx, row in enumerate(ranked):
        mean = as_float(row, "mean_affinity_kcal_mol")
        sd = as_float(row, "sd_affinity_kcal_mol")
        ax.errorbar(mean, idx, xerr=sd, fmt="D", color=TARGET_COLORS[row["receptor"]], ecolor="#111827", markersize=4.7, capsize=2.6, elinewidth=0.9, zorder=4)

    seed_offsets = {1001: -0.16, 1002: 0.0, 1003: 0.16}
    seed_colors = {1001: "#FFFFFF", 1002: "#E0F2FE", 1003: "#FEF3C7"}
    for row in replicates:
        seed = int(row["seed"])
        yy = y_lookup[(row["receptor"], row["ligand"])] + seed_offsets[seed]
        ax.scatter(
            as_float(row, "affinity_kcal_mol"),
            yy,
            s=21,
            color=seed_colors[seed],
            edgecolor="#1F2937",
            linewidth=0.55,
            zorder=5,
            label=str(seed),
        )

    ax.set_yticks(y)
    ax.set_yticklabels([fill(label, 27) for label in labels], fontproperties=small)
    ax.invert_yaxis()
    ax.set_xlim(-10.1, -5.45)
    ax.set_xlabel("Vina score (kcal/mol)", fontproperties=font)
    if show_title:
        ax.set_title("All 18 pairs: seed-level points and mean ± SD", fontproperties=font, pad=7)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.tick_params(axis="x", labelsize=7.5)

    handles = []
    labels_seen = []
    for seed in [1001, 1002, 1003]:
        handles.append(ax.scatter([], [], s=24, color=seed_colors[seed], edgecolor="#1F2937", linewidth=0.55))
        labels_seen.append(f"seed {seed}")
    handles.append(ax.scatter([], [], s=26, marker="D", color="#64748B"))
    labels_seen.append("mean")
    ax.legend(handles, labels_seen, frameon=False, loc=legend_loc, fontsize=7.3, ncol=2)


def draw_sd_bars(ax: plt.Axes, summary: list[dict[str, str]], font: FontProperties, small: FontProperties) -> None:
    ranked = sorted(summary, key=lambda row: as_float(row, "sd_affinity_kcal_mol"), reverse=True)
    labels = [pair_label(row) for row in ranked]
    values = [as_float(row, "sd_affinity_kcal_mol") for row in ranked]
    y = np.arange(len(ranked))
    colors = ["#F59E0B" if row["stability"] == "acceptable" else "#64748B" for row in ranked]
    ax.barh(y, values, color=colors, alpha=0.88, height=0.64)
    ax.axvline(0.30, color="#DC2626", linewidth=0.9, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels([fill(label, 27) for label in labels], fontproperties=small)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.34)
    ax.set_xlabel("SD across 3 seeds (kcal/mol)", fontproperties=font)
    ax.set_title("Seed variability ranked", fontproperties=font, pad=8)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.text(0.296, -0.84, "stable cutoff", ha="right", va="center", fontproperties=small, color="#DC2626")
    for yi, val in enumerate(values):
        if val > 0.26:
            x_text = val - 0.012
            ha = "right"
            color = "white"
        else:
            x_text = val + 0.008
            ha = "left"
            color = "#374151"
        ax.text(x_text, yi, f"{val:.3f}", ha=ha, va="center", fontproperties=small, color=color)


def draw_sd_heatmap(ax: plt.Axes, summary: list[dict[str, str]], font: FontProperties, small: FontProperties) -> None:
    sd = matrix_from_summary(summary, "sd_affinity_kcal_mol")
    im = ax.imshow(sd, cmap="YlOrBr", vmin=0, vmax=0.35, aspect="auto")
    ax.set_xticks(np.arange(len(LIGANDS)))
    ax.set_xticklabels([LIGAND_LABELS[x] for x in LIGANDS], fontproperties=small, rotation=18, ha="right")
    ax.set_yticks(np.arange(len(RECEPTORS)))
    ax.set_yticklabels(RECEPTORS, fontproperties=small)
    ax.set_title("SD matrix", fontproperties=font, pad=7)
    for i in range(sd.shape[0]):
        for j in range(sd.shape[1]):
            color = "white" if sd[i, j] >= 0.18 else "#111827"
            ax.text(j, i, f"{sd[i, j]:.3f}", ha="center", va="center", fontproperties=small, color=color)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.028)
    cbar.ax.set_ylabel("SD", rotation=270, va="bottom", fontproperties=small)
    cbar.ax.tick_params(labelsize=7.2)


def figure_seed_stability(summary: list[dict[str, str]], replicates: list[dict[str, str]]) -> list[Path]:
    font, small, _, panel = setup_style()
    fig = plt.figure(figsize=(15.5, 9.0), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1.52, 1.0], height_ratios=[1.0, 0.88], wspace=0.15, hspace=0.24)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    draw_seed_points(ax_a, summary, replicates, font, small)
    draw_sd_bars(ax_b, summary, font, small)
    draw_sd_heatmap(ax_c, summary, font, small)
    for label, ax in zip(["A", "B", "C"], [ax_a, ax_b, ax_c]):
        panel_label(ax, label, panel)
    return save_figure(fig, "fig_docking_seed_stability_all")


def figure_seed_sd_ranked(summary: list[dict[str, str]]) -> list[Path]:
    font, small, _, _ = setup_style()
    fig, ax = plt.subplots(figsize=(7.8, 5.6), constrained_layout=True)
    draw_sd_bars(ax, summary, font, small)
    ax.set_title("Docking seed variability across all 18 pairs", fontproperties=font, pad=8)
    return save_figure(fig, "fig_docking_seed_sd_ranked")


def main() -> None:
    summary, replicates = load_data()
    overview = figure_overview(summary, replicates)
    stability = figure_seed_stability(summary, replicates)
    sd_ranked = figure_seed_sd_ranked(summary)
    print("Generated:")
    for path in [*overview, *stability, *sd_ranked]:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
