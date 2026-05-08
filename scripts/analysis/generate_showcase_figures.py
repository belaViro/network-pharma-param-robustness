"""
Generate showcase-style scientific figures for visual review.

These figures are intentionally separate from the numbered manuscript figures.
They use existing outputs and are written to:
    src/outputs/paper_figures/showcase/
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[2]
OUT_DIR = BASE_DIR / "src" / "outputs"
SENS_DIR = OUT_DIR / "sensitivity"
COMP_DIR = SENS_DIR / "compound_sensitivity"
VALID_DIR = SENS_DIR / "robustness_validation"
DOCK_SUMMARY = BASE_DIR / "分子对接" / "results" / "validation_runs" / "docking_replicates_summary.csv"
FIG_DIR = OUT_DIR / "paper_figures" / "showcase"
FIG_DIR.mkdir(parents=True, exist_ok=True)


COLORS = {
    "ink": "#111827",
    "slate": "#475569",
    "muted": "#64748B",
    "grid": "#E2E8F0",
    "blue": "#2563EB",
    "cyan": "#0891B2",
    "teal": "#0F766E",
    "green": "#15803D",
    "lime": "#65A30D",
    "amber": "#D97706",
    "orange": "#EA580C",
    "red": "#DC2626",
    "purple": "#7C3AED",
    "violet": "#6D28D9",
    "pale_blue": "#EFF6FF",
    "pale_green": "#F0FDF4",
    "pale_orange": "#FFF7ED",
    "pale_red": "#FEF2F2",
    "pale_purple": "#F5F3FF",
}

FONT = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=9)
FONT_SMALL = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=7.3)
FONT_TITLE = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=16, weight="bold")
FONT_BLOCK = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=9.2, weight="bold")

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def save(fig: plt.Figure, stem: str) -> list[str]:
    paths: list[str] = []
    for ext in ("png", "svg", "pdf"):
        path = FIG_DIR / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.08)
        paths.append(str(path.relative_to(BASE_DIR)))
    plt.close(fig)
    return paths


def write_design_philosophy() -> None:
    text = """# Systems Lattice

Systems Lattice treats scientific evidence as architecture: every count, edge, and threshold becomes a structural member in a carefully engineered visual object. Space is not decoration; it is the method by which uncertainty is separated, compressed, and made legible. The final work should feel meticulously crafted, as though every line, node, and label has been placed after repeated inspection by someone at the top of their field.

Color is used as evidence coding rather than ornament. Blue carries inputs, orange carries parameter stress, violet carries compound perturbation, green carries retained robustness, and red is reserved for audit rather than proof. The palette must look disciplined, not festive, with painstaking attention to contrast and hierarchy.

Scale communicates scientific weight. Large bands and nodes show data reduction and convergence; small labels are restrained anchors, never explanatory paragraphs. The viewer should understand the story through flow, enclosure, density, and recurrence before reading any text. The design should look like the product of deep expertise rather than a chart template.

Rhythm comes from repeated modules, gridded alignments, and controlled asymmetry. The images should reward close reading: first as striking visual objects, then as precise scientific maps. Master-level execution means no accidental overlaps, no crowded labels, no arbitrary embellishment, and no visual claim unsupported by the workflow.
"""
    (FIG_DIR / "showcase_design_philosophy.md").write_text(text, encoding="utf-8")


def ribbon(ax: plt.Axes, x0: float, y0: float, x1: float, y1: float, width: float, color: str, alpha: float = 0.45) -> None:
    dx = (x1 - x0) * 0.42
    verts = [
        (x0, y0 + width / 2),
        (x0 + dx, y0 + width / 2),
        (x1 - dx, y1 + width / 2),
        (x1, y1 + width / 2),
        (x1, y1 - width / 2),
        (x1 - dx, y1 - width / 2),
        (x0 + dx, y0 - width / 2),
        (x0, y0 - width / 2),
        (x0, y0 + width / 2),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))


def node(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, value: str, color: str, face: str) -> None:
    shadow = patches.FancyBboxPatch(
        (x + 0.04, y - 0.04),
        w,
        h,
        boxstyle="round,pad=0.03,rounding_size=0.10",
        linewidth=0,
        facecolor="#CBD5E1",
        alpha=0.32,
        zorder=2,
    )
    box = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.03,rounding_size=0.10",
        linewidth=1.4,
        edgecolor=color,
        facecolor=face,
        zorder=3,
    )
    ax.add_patch(shadow)
    ax.add_patch(box)
    ax.text(x + 0.14, y + h - 0.18, title, fontproperties=FONT_SMALL, color=COLORS["slate"], va="top", zorder=4)
    value_size = 14 if "\n" in value or len(value) > 8 else 18
    ax.text(
        x + w / 2,
        y + h * 0.40,
        value,
        fontproperties=FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=value_size, weight="bold"),
        color=color,
        ha="center",
        va="center",
        linespacing=1.15,
        zorder=4,
    )


def fig_convergence_alluvial() -> list[str]:
    fig, ax = plt.subplots(figsize=(13.8, 6.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")

    ax.text(6, 6.72, "双扰动筛选收敛图", fontproperties=FONT_TITLE, ha="center", color=COLORS["ink"])
    ax.text(
        6,
        6.35,
        "Targets and functional terms are compressed by parameter robustness and compound leave-one-out retention.",
        fontproperties=FONT,
        ha="center",
        color=COLORS["muted"],
    )

    target_nodes = [
        ("预测并集靶点", "684", 0.8, 4.70, COLORS["blue"], COLORS["pale_blue"]),
        ("疾病交集靶点", "286", 3.0, 4.70, COLORS["cyan"], "#ECFEFF"),
        ("PR-QI 靶点", "50", 5.3, 4.70, COLORS["orange"], COLORS["pale_orange"]),
        ("双重稳健靶点", "29", 7.8, 4.70, COLORS["green"], COLORS["pale_green"]),
        ("Tier 1", "STAT3", 10.0, 4.70, COLORS["red"], COLORS["pale_red"]),
    ]
    term_nodes = [
        ("GO BP/KEGG 功能条目", "70", 5.3, 1.65, COLORS["teal"], "#F0FDFA"),
        ("双重稳健功能条目", "58", 7.8, 1.65, COLORS["lime"], "#F7FEE7"),
        ("机制链条", "EGFR/SRC\nSTAT3", 10.0, 1.65, COLORS["purple"], COLORS["pale_purple"]),
    ]

    centers: dict[str, tuple[float, float]] = {}
    for title, value, x, y, c, face in target_nodes + term_nodes:
        node(ax, x, y, 1.65, 0.98, title, value, c, face)
        centers[title] = (x + 1.65, y + 0.49)
        centers[title + "_left"] = (x, y + 0.49)

    ribbon(ax, centers["预测并集靶点"][0], centers["预测并集靶点"][1], centers["疾病交集靶点_left"][0], centers["疾病交集靶点_left"][1], 0.40, COLORS["blue"], 0.28)
    ribbon(ax, centers["疾病交集靶点"][0], centers["疾病交集靶点"][1], centers["PR-QI 靶点_left"][0], centers["PR-QI 靶点_left"][1], 0.30, COLORS["cyan"], 0.30)
    ribbon(ax, centers["PR-QI 靶点"][0], centers["PR-QI 靶点"][1], centers["双重稳健靶点_left"][0], centers["双重稳健靶点_left"][1], 0.26, COLORS["orange"], 0.35)
    ribbon(ax, centers["双重稳健靶点"][0], centers["双重稳健靶点"][1], centers["Tier 1_left"][0], centers["Tier 1_left"][1], 0.18, COLORS["green"], 0.42)

    ribbon(ax, 4.65, 4.75, centers["GO BP/KEGG 功能条目_left"][0], centers["GO BP/KEGG 功能条目_left"][1], 0.18, COLORS["teal"], 0.24)
    ribbon(ax, centers["GO BP/KEGG 功能条目"][0], centers["GO BP/KEGG 功能条目"][1], centers["双重稳健功能条目_left"][0], centers["双重稳健功能条目_left"][1], 0.34, COLORS["teal"], 0.36)
    ribbon(ax, centers["双重稳健功能条目"][0], centers["双重稳健功能条目"][1], centers["机制链条_left"][0], centers["机制链条_left"][1], 0.26, COLORS["lime"], 0.36)
    ribbon(ax, centers["Tier 1"][0] - 1.65, centers["Tier 1"][1] - 0.30, centers["机制链条_left"][0], centers["机制链条_left"][1] + 0.20, 0.16, COLORS["purple"], 0.28)

    ax.text(2.0, 3.95, "disease intersection", fontproperties=FONT_SMALL, color=COLORS["muted"], ha="center")
    ax.text(5.8, 3.95, "parameter robustness", fontproperties=FONT_SMALL, color=COLORS["muted"], ha="center")
    ax.text(8.35, 3.95, "compound retention", fontproperties=FONT_SMALL, color=COLORS["muted"], ha="center")
    ax.text(0.55, 0.55, "注：GO BP 为生物过程条目，KEGG 为通路；图中统一表述为功能富集条目。", fontproperties=FONT_SMALL, color=COLORS["muted"])

    return save(fig, "showcase_01_convergence_alluvial")


def normalize(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce").fillna(0.0)
    vmin, vmax = float(values.min()), float(values.max())
    if math.isclose(vmin, vmax):
        return pd.Series([0.5] * len(values), index=values.index)
    return (values - vmin) / (vmax - vmin)


def fig_target_evidence_radar() -> list[str]:
    dual = pd.read_csv(COMP_DIR / "dual_robust_targets.csv")
    pubmed = pd.read_csv(VALID_DIR / "target_pubmed_counts.csv")
    dock = pd.read_csv(DOCK_SUMMARY)

    focus = ["STAT3", "EGFR", "SRC", "TNF", "HSP90AA1", "BCL2", "PTGS2", "MMP9"]
    df = dual[dual["gene"].isin(focus)].copy()
    df["gene"] = pd.Categorical(df["gene"], categories=focus, ordered=True)
    df = df.sort_values("gene").reset_index(drop=True)
    df = df.merge(pubmed[["gene", "pubmed_count"]], on="gene", how="left")
    docking_score = dock.groupby("receptor")["mean_affinity_kcal_mol"].min().abs().rename("docking_score")
    df = df.merge(docking_score, left_on="gene", right_index=True, how="left")

    metrics = [
        ("参数频率", "full_frequency"),
        ("网络 degree", "full_avg_degree_norm"),
        ("LOO 保留", "loo_q1_retention_rate"),
        ("最低 LOO 频率", "min_loo_frequency"),
        ("PubMed 背景", "pubmed_count"),
        ("Docking 支持", "docking_score"),
    ]
    scaled = pd.DataFrame({"gene": df["gene"].astype(str)})
    for label, col in metrics:
        values = np.log1p(df[col]) if col == "pubmed_count" else df[col]
        scaled[label] = normalize(values)
        if col == "docking_score":
            scaled.loc[df[col].isna(), label] = np.nan

    matrix = scaled[[m[0] for m in metrics]].to_numpy(dtype=float)
    composite = np.nanmean(matrix, axis=1)
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#F1F5F9")
    masked = np.ma.masked_invalid(matrix)

    fig = plt.figure(figsize=(12.5, 6.6))
    gs = fig.add_gridspec(1, 2, left=0.11, right=0.94, top=0.80, bottom=0.18, width_ratios=[4.4, 1.45], wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])
    bx = fig.add_subplot(gs[0, 1])

    fig.text(0.5, 0.93, "核心靶点多维证据矩阵", fontproperties=FONT_TITLE, ha="center", color=COLORS["ink"])
    fig.text(
        0.5,
        0.885,
        "Parameter frequency, network centrality, compound retention, literature background, and docking support.",
        fontproperties=FONT,
        ha="center",
        color=COLORS["muted"],
    )

    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(metrics)), [m[0] for m in metrics], fontproperties=FONT)
    ax.set_yticks(np.arange(len(focus)), focus, fontproperties=FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=10, weight="bold"))
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=8)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xticks(np.arange(-0.5, len(metrics), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(focus), 1), minor=True)
    ax.grid(which="minor", color="white", lw=2.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    raw_lookup = df.set_index(df["gene"].astype(str))
    for i, gene in enumerate(focus):
        for j, (label, col) in enumerate(metrics):
            raw_value = raw_lookup.loc[gene, col]
            if pd.isna(raw_value) or pd.isna(matrix[i, j]):
                text = "NA"
                color = COLORS["muted"]
            elif col == "pubmed_count":
                text = str(int(raw_lookup.loc[gene, "pubmed_count"]))
                color = "white" if matrix[i, j] > 0.55 else COLORS["ink"]
            elif col == "docking_score":
                text = f"{raw_value:.1f}"
                color = "white" if matrix[i, j] > 0.55 else COLORS["ink"]
            else:
                text = f"{raw_value:.2f}"
                color = "white" if matrix[i, j] > 0.55 else COLORS["ink"]
            ax.text(j, i, text, ha="center", va="center", fontproperties=FONT_SMALL, color=color)

    for row_index, edge in [(0, COLORS["red"]), (1, COLORS["blue"]), (2, COLORS["teal"])]:
        ax.add_patch(patches.Rectangle((-0.5, row_index - 0.5), len(metrics), 1, fill=False, edgecolor=edge, linewidth=1.6))

    bx.barh(np.arange(len(focus)), composite, color=["#DC2626", "#2563EB", "#0F766E", "#7C3AED", "#EA580C", "#15803D", "#475569", "#0891B2"], alpha=0.88)
    bx.set_xlim(0, 1)
    bx.set_ylim(ax.get_ylim())
    bx.set_xticks([0, 0.5, 1.0], ["0", "0.5", "1.0"], fontproperties=FONT_SMALL)
    bx.set_yticks(np.arange(len(focus)))
    bx.set_yticklabels([])
    bx.tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)
    bx.set_title("综合证据\n均值", fontproperties=FONT, color=COLORS["slate"], pad=10)
    bx.grid(axis="x", color=COLORS["grid"], lw=0.8)
    for spine in bx.spines.values():
        spine.set_visible(False)
    for i, val in enumerate(composite):
        bx.text(min(val + 0.025, 0.98), i, f"{val:.2f}", va="center", ha="left" if val < 0.92 else "right", fontproperties=FONT_SMALL, color=COLORS["slate"])

    cax = fig.add_axes([0.955, 0.23, 0.012, 0.44])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=7, length=0)
    cb.outline.set_visible(False)
    cb.ax.set_title("0-1", fontproperties=FONT_SMALL, color=COLORS["muted"], pad=6)

    fig.text(
        0.08,
        0.09,
        "注：PubMed 背景按 log1p 归一化；Docking 支持仅对已完成分子对接的靶点显示，NA 不代表无作用证据。",
        fontproperties=FONT_SMALL,
        color=COLORS["muted"],
    )

    return save(fig, "showcase_02_target_evidence_matrix")


def rounded_label(ax: plt.Axes, xy: tuple[float, float], w: float, h: float, text: str, color: str, face: str, size: float = 9.5) -> None:
    box = patches.FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        linewidth=1.35,
        edgecolor=color,
        facecolor=face,
        zorder=3,
    )
    ax.add_patch(box)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, fontproperties=FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=size, weight="bold"), ha="center", va="center", color=COLORS["ink"], zorder=4)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str, rad: float = 0.0, lw: float = 1.6) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=14, connectionstyle=f"arc3,rad={rad}"),
        zorder=2,
    )


def fig_mechanism_model() -> list[str]:
    fig, ax = plt.subplots(figsize=(13.8, 7.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.axis("off")

    ax.text(6, 6.68, "EGFR/SRC-STAT3 候选机制模型", fontproperties=FONT_TITLE, ha="center", color=COLORS["ink"])
    ax.text(
        6,
        6.33,
        "A convergent hypothesis supported by dual robustness, functional enrichment, external audit, and seed-stable docking.",
        fontproperties=FONT,
        ha="center",
        color=COLORS["muted"],
    )

    rounded_label(ax, (0.55, 4.55), 2.2, 1.28, "草莓叶\n活性成分集合", COLORS["blue"], COLORS["pale_blue"], 10)
    rounded_label(ax, (0.55, 2.75), 2.2, 1.18, "高影响成分\nCatechin / 强心苷类", COLORS["orange"], COLORS["pale_orange"], 8.8)
    rounded_label(ax, (3.55, 4.80), 1.50, 0.86, "EGFR", COLORS["purple"], COLORS["pale_purple"], 12)
    rounded_label(ax, (3.55, 3.70), 1.50, 0.86, "SRC", COLORS["purple"], COLORS["pale_purple"], 12)
    rounded_label(ax, (3.55, 2.60), 1.50, 0.86, "TNF / PTGS2\nMMP9", COLORS["red"], COLORS["pale_red"], 8.8)

    rounded_label(ax, (6.20, 4.78), 1.70, 0.82, "MAPK", COLORS["teal"], "#F0FDFA", 10)
    rounded_label(ax, (6.20, 3.78), 1.70, 0.82, "PI3K-Akt", COLORS["teal"], "#F0FDFA", 10)
    rounded_label(ax, (6.20, 2.78), 1.70, 0.82, "cytokine / NF-kB", COLORS["teal"], "#F0FDFA", 8.2)
    rounded_label(ax, (8.78, 3.84), 1.58, 1.08, "STAT3\nTier 1", COLORS["green"], COLORS["pale_green"], 12)
    rounded_label(ax, (10.62, 4.95), 1.08, 0.74, "炎症", COLORS["slate"], "#F8FAFC", 9.2)
    rounded_label(ax, (10.62, 3.88), 1.08, 0.74, "凋亡", COLORS["slate"], "#F8FAFC", 9.2)
    rounded_label(ax, (10.62, 2.82), 1.08, 0.74, "ROS 响应", COLORS["slate"], "#F8FAFC", 9.2)

    arrow(ax, (2.75, 5.20), (3.55, 5.23), COLORS["blue"])
    arrow(ax, (2.75, 3.34), (3.55, 4.05), COLORS["orange"], rad=0.10)
    arrow(ax, (5.05, 5.23), (6.20, 5.18), COLORS["purple"])
    arrow(ax, (5.05, 4.13), (6.20, 4.18), COLORS["purple"])
    arrow(ax, (5.05, 3.03), (6.20, 3.18), COLORS["red"])
    arrow(ax, (7.90, 5.18), (8.78, 4.52), COLORS["teal"], rad=-0.08)
    arrow(ax, (7.90, 4.18), (8.78, 4.38), COLORS["teal"])
    arrow(ax, (7.90, 3.18), (8.78, 4.15), COLORS["teal"], rad=0.08)
    arrow(ax, (10.36, 4.38), (10.62, 5.32), COLORS["green"], rad=0.04)
    arrow(ax, (10.36, 4.38), (10.62, 4.25), COLORS["green"])
    arrow(ax, (10.36, 4.38), (10.62, 3.18), COLORS["green"], rad=-0.04)

    # Evidence tags
    tags = [
        (3.36, 5.83, "docking +", COLORS["purple"]),
        (3.35, 3.42, "dual robust", COLORS["green"]),
        (6.02, 2.43, "GO BP/KEGG", COLORS["teal"]),
        (8.62, 5.13, "frequency=1.00", COLORS["green"]),
        (8.62, 3.48, "LOO=1.00", COLORS["green"]),
        (10.32, 2.35, "audit only", COLORS["red"]),
    ]
    for x, y, text, color in tags:
        ax.text(
            x,
            y,
            text,
            fontproperties=FONT_SMALL,
            color=color,
            bbox=dict(boxstyle="round,pad=0.22,rounding_size=0.08", facecolor="white", edgecolor=color, linewidth=0.8),
            zorder=5,
        )

    ax.text(
        0.55,
        0.62,
        "注：该图为候选机制模型，不表示已完成因果验证；GO BP 为生物过程条目，不等同于通路。",
        fontproperties=FONT_SMALL,
        color=COLORS["muted"],
    )

    return save(fig, "showcase_03_mechanism_model")


def main() -> None:
    write_design_philosophy()
    created = {
        "alluvial": fig_convergence_alluvial(),
        "radar": fig_target_evidence_radar(),
        "mechanism": fig_mechanism_model(),
    }
    print("Created showcase figures:")
    for key, paths in created.items():
        print(key)
        for path in paths:
            print(f"  {path}")


if __name__ == "__main__":
    main()
