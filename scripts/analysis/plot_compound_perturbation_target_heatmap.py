from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[2]
OUT_ROOT = ROOT / "src" / "outputs"
COMP_DIR = OUT_ROOT / "sensitivity" / "compound_sensitivity"
PAPER_DIR = OUT_ROOT / "paper_figures"
FIG_DIR = COMP_DIR / "figures"

DUAL_TARGETS = COMP_DIR / "dual_robust_targets.csv"
LOO_TARGETS = COMP_DIR / "leave_one_out_target_robustness.csv"
COMPOUND_INFLUENCE = COMP_DIR / "compound_influence_scores.csv"


SHORT_COMPOUND_NAMES = {
    "3beta9xi-3-beta-D-Glucopyranosyloxy-14-hydroxycard-2022-enolide": "Cardiac\nglycoside",
    "Indole-3-acrylic_acid": "Indole-3-\nacrylic acid",
    "2-Methylglutaric_acid": "2-Methyl-\nglutaric acid",
    "A-Ketoglutaric_Acid": "alpha-Keto-\nglutaric acid",
    "DL-Tryptophan": "DL-\nTryptophan",
    "Procyanidin_B1": "Procyanidin\nB1",
    "DL-Malic_acid": "DL-Malic\nacid",
    "Citric_acid": "Citric\nacid",
    "L-Allo-Isoleucine": "L-Allo-\nIsoleucine",
}


def setup_fonts() -> tuple[FontProperties, FontProperties, FontProperties, FontProperties]:
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial Unicode MS", "DejaVu Sans"]
    font_cn = FontProperties(family="SimHei", size=9)
    font_cn_small = FontProperties(family="SimHei", size=7.5)
    font_cn_tiny = FontProperties(family="SimHei", size=6.4)
    font_cn_title = FontProperties(family="SimHei", size=11.5, weight="bold")
    return font_cn, font_cn_small, font_cn_tiny, font_cn_title


def save_all(fig: plt.Figure, stem: str, paper_stem: str | None = None) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    PAPER_DIR.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for out_dir, name in [(FIG_DIR, stem), (PAPER_DIR, paper_stem or stem)]:
        for ext, kwargs in [
            ("png", {"dpi": 300, "bbox_inches": "tight"}),
            ("svg", {"bbox_inches": "tight"}),
            ("pdf", {"bbox_inches": "tight"}),
        ]:
            path = out_dir / f"{name}.{ext}"
            fig.savefig(path, **kwargs)
            outputs.append(path)
    return outputs


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in [DUAL_TARGETS, LOO_TARGETS, COMPOUND_INFLUENCE]:
        if not path.exists():
            raise FileNotFoundError(path)
    return pd.read_csv(DUAL_TARGETS), pd.read_csv(LOO_TARGETS), pd.read_csv(COMPOUND_INFLUENCE)


def build_status_matrix(
    dual: pd.DataFrame,
    loo: pd.DataFrame,
    influence: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    compound_order = influence.sort_values("compound_influence_score", ascending=False)["drop_compound"].tolist()

    row_order = dual.sort_values(
        ["is_dual_robust", "min_loo_frequency", "loo_q1_retention_rate", "full_frequency", "full_avg_degree_norm"],
        ascending=[False, False, False, False, False],
    )["gene"].tolist()

    lookup = {
        (str(row.drop_compound), str(row.gene)): row
        for row in loo.itertuples(index=False)
    }

    matrix = pd.DataFrame(index=row_order, columns=compound_order, dtype=int)
    freq_matrix = pd.DataFrame(index=row_order, columns=compound_order, dtype=float)

    for gene in row_order:
        for compound in compound_order:
            row = lookup.get((compound, gene))
            if row is None or float(row.frequency) <= 0:
                status = 0  # completely lost from TopK across that LOO parameter grid
                freq = 0.0
            elif "Q1" in str(row.quadrant):
                status = 2  # retained as PR-QI
                freq = float(row.frequency)
            else:
                status = 1  # present, but no longer PR-QI
                freq = float(row.frequency)
            matrix.loc[gene, compound] = status
            freq_matrix.loc[gene, compound] = freq

    return matrix, freq_matrix, compound_order


def compound_label(compound: str) -> str:
    return SHORT_COMPOUND_NAMES.get(compound, compound.replace("_", "\n"))


def figure_heatmap() -> list[Path]:
    dual, loo, influence = load_inputs()
    matrix, _, compound_order = build_status_matrix(dual, loo, influence)
    font_cn, font_cn_small, font_cn_tiny, font_cn_title = setup_fonts()

    values = matrix.astype(int).to_numpy()
    cmap = ListedColormap(["#C43C4B", "#E5A823", "#3778BF"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    fig, ax = plt.subplots(figsize=(10.6, 10.9), constrained_layout=True)
    ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(np.arange(len(compound_order)))
    ax.set_xticklabels([compound_label(c) for c in compound_order], rotation=45, ha="right", fontproperties=font_cn_tiny)
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index.tolist(), fontproperties=font_cn_tiny)

    dual_lookup = dual.set_index("gene")["is_dual_robust"].to_dict()
    for tick in ax.get_yticklabels():
        gene = tick.get_text()
        if bool(dual_lookup.get(gene, False)):
            tick.set_color("#1F5F99")
            tick.set_fontweight("bold")
        else:
            tick.set_color("#8A1F2D")

    for x in np.arange(-0.5, len(compound_order), 1):
        ax.axvline(x, color="white", linewidth=0.45, zorder=3)
    for y in np.arange(-0.5, len(matrix.index), 1):
        ax.axhline(y, color="white", linewidth=0.35, zorder=3)

    n_dual = int(dual["is_dual_robust"].sum())
    if 0 < n_dual < len(matrix.index):
        ax.axhline(n_dual - 0.5, color="#111111", linewidth=1.2)

    ax.set_xlabel("剔除的化合物（按 compound influence score 降序）", fontproperties=font_cn, labelpad=10)
    ax.set_title("成分 leave-one-out 下 PR-QI 靶点保留状态热图", fontproperties=font_cn_title, pad=12)

    legend_items = [
        Patch(facecolor="#3778BF", edgecolor="none", label="仍为 PR-QI"),
        Patch(facecolor="#E5A823", edgecolor="none", label="仍出现但不再为 PR-QI"),
        Patch(facecolor="#C43C4B", edgecolor="none", label="参数频率降为 0"),
        Patch(facecolor="none", edgecolor="none", label="基因名蓝色：dual-robust"),
        Patch(facecolor="none", edgecolor="none", label="基因名红色：成分敏感"),
    ]
    ax.legend(handles=legend_items, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, prop=font_cn_small)

    return save_all(fig, "fig_compound_loo_target_status_heatmap", "fig4_compound_loo_target_status_heatmap")


def figure_retention_ranking() -> list[Path]:
    dual, _, _ = load_inputs()
    font_cn, font_cn_small, font_cn_tiny, font_cn_title = setup_fonts()

    plot_df = dual.sort_values(
        ["is_dual_robust", "min_loo_frequency", "loo_q1_retention_rate", "full_frequency"],
        ascending=[True, True, True, True],
    ).reset_index(drop=True)

    colors = np.where(plot_df["is_dual_robust"], "#3778BF", "#C43C4B")
    y = np.arange(len(plot_df))

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(11.2, 11.4),
        sharey=True,
        gridspec_kw={"width_ratios": [1.05, 1.15], "wspace": 0.08},
        constrained_layout=True,
    )

    ax1.barh(y, plot_df["loo_q1_retention_rate"], color=colors, alpha=0.88, edgecolor="white", linewidth=0.4)
    ax1.axvline(0.80, color="#555555", linestyle="--", linewidth=0.9, alpha=0.7)
    ax1.set_xlim(0.78, 1.01)
    ax1.set_xlabel("LOO PR-QI 保留率", fontproperties=font_cn)
    ax1.set_yticks(y)
    ax1.set_yticklabels(plot_df["gene"], fontproperties=font_cn_tiny)
    ax1.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.25)

    ax2.scatter(plot_df["min_loo_frequency"], y, s=42, c=colors, alpha=0.90, edgecolor="white", linewidth=0.35)
    ax2.axvline(0, color="#555555", linestyle="--", linewidth=0.9, alpha=0.7)
    ax2.set_xlim(-0.02, max(0.76, float(plot_df["min_loo_frequency"].max()) * 1.08))
    ax2.set_xlabel("最低 LOO 参数频率", fontproperties=font_cn)
    ax2.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.25)
    ax2.tick_params(axis="y", left=False, labelleft=False)

    for ax in (ax1, ax2):
        for label in ax.get_xticklabels():
            label.set_fontproperties(font_cn_small)

    for tick in ax1.get_yticklabels():
        gene = tick.get_text()
        if gene in {"STAT3", "EGFR", "SRC", "TNF", "HSP90AA1", "BCL2"}:
            tick.set_fontweight("bold")
            tick.set_color("#1F5F99")
        elif not bool(plot_df.loc[plot_df["gene"] == gene, "is_dual_robust"].iloc[0]):
            tick.set_color("#8A1F2D")

    legend_items = [
        Patch(facecolor="#3778BF", label="dual-robust"),
        Patch(facecolor="#C43C4B", label="成分敏感"),
    ]
    ax2.legend(handles=legend_items, loc="lower right", frameon=False, prop=font_cn_small)
    fig.suptitle("全成分 PR-QI 靶点的成分扰动保留率排序", fontproperties=font_cn_title, y=0.995)

    return save_all(fig, "fig_compound_loo_target_retention_ranking", "figS_compound_loo_target_retention_ranking")


def main() -> None:
    outputs = figure_heatmap()
    outputs.extend(figure_retention_ranking())
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
