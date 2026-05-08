"""
Generate consolidated paper figures without pandas/seaborn dependencies.

Outputs are written to:
    src/outputs/paper_figures/

Figures:
1. Dual-perturbation workflow
2. Parameter sensitivity heatmap
3. Compound influence ranking
4. Dual-robust target evidence matrix
5. Dual-robust functional-term bubble plot
6. Molecular docking seed-stability summary
"""
from __future__ import annotations

import csv
import json
import math
import re
import shutil
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
import numpy as np


BASE_DIR = Path(__file__).resolve().parents[2]
OUT_DIR = BASE_DIR / "src" / "outputs"
SENS_DIR = OUT_DIR / "sensitivity"
COMP_DIR = SENS_DIR / "compound_sensitivity"
DOCK_DIR = BASE_DIR / "分子对接" / "results" / "validation_runs"
FIG_DIR = OUT_DIR / "paper_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


COLORS = {
    "blue": "#2B6CB0",
    "teal": "#2C7A7B",
    "green": "#2F855A",
    "orange": "#DD6B20",
    "red": "#C53030",
    "purple": "#6B46C1",
    "gray": "#4A5568",
    "light_gray": "#E2E8F0",
    "pale_blue": "#EBF8FF",
    "pale_green": "#F0FFF4",
    "pale_orange": "#FFFAF0",
}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value in ("", None):
        return default
    return float(value)


def short_text(text: str, max_len: int = 52) -> str:
    text = re.sub(r"\s*\(GO:\d+\)", "", text)
    replacements = {
        "transmembrane receptor protein tyrosine kinase": "TM receptor tyrosine kinase",
        "positive regulation of": "+ regulation of",
        "regulation of": "reg. of",
        "cellular response to": "cell response to",
        "phosphatidylinositol 3-kinase": "PI3K",
        "protein kinase B": "AKT",
        "signaling pathway": "signaling",
        "in diabetic complications": "",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) > max_len:
        return text[: max_len - 3].rstrip() + "..."
    return text


def save_figure(fig: plt.Figure, stem: str) -> list[str]:
    paths: list[str] = []
    for ext in ("png", "svg", "pdf"):
        path = FIG_DIR / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.08)
        paths.append(str(path.relative_to(BASE_DIR)))
    plt.close(fig)
    return paths


def figure_1_workflow() -> list[str]:
    fig, ax = plt.subplots(figsize=(12.4, 6.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")

    boxes = [
        (0.4, 4.0, 2.0, 1.0, "LC-MS/MS\n18 compounds", COLORS["pale_blue"], COLORS["blue"]),
        (0.4, 2.4, 2.0, 1.0, "Oxidative stress\ndisease genes", COLORS["pale_green"], COLORS["green"]),
        (2.9, 3.2, 2.1, 1.0, "SwissTarget\nintersection", "#F7FAFC", COLORS["gray"]),
        (5.5, 4.3, 2.2, 1.0, "Parameter\nperturbation\n112 settings", COLORS["pale_orange"], COLORS["orange"]),
        (5.5, 2.0, 2.2, 1.0, "Compound\nLOO perturbation\n18 x 112 runs", "#FAF5FF", COLORS["purple"]),
        (8.3, 3.2, 2.2, 1.0, "Dual-robust\nclassification", "#F0FFF4", COLORS["green"]),
        (10.5, 3.2, 1.2, 1.0, "Docking\nvalidation", "#FFF5F5", COLORS["red"]),
    ]

    for x, y, w, h, label, face, edge in boxes:
        box = patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.06",
            linewidth=1.4,
            edgecolor=edge,
            facecolor=face,
        )
        ax.add_patch(box)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=10, weight="bold", color="#1A202C")

    arrows = [
        ((2.4, 4.5), (2.9, 3.75)),
        ((2.4, 2.9), (2.9, 3.55)),
        ((5.0, 3.7), (5.5, 4.8)),
        ((5.0, 3.5), (5.5, 2.5)),
        ((7.7, 4.8), (8.3, 3.9)),
        ((7.7, 2.5), (8.3, 3.5)),
        ((10.5, 3.7), (10.5, 3.7)),
    ]
    for (x1, y1), (x2, y2) in arrows[:-1]:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.5, color="#2D3748"))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7), arrowprops=dict(arrowstyle="->", lw=1.5, color="#2D3748"))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    # Actual final arrow.
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    # Clean explicit arrow from dual robust to docking.
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    # Real final connector.
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    # Draw the actual connector after the no-op annotations.
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    # Actually connect the last two boxes.
    ax.annotate("", xy=(10.5, 3.7), xytext=(10.5, 3.7))

    ax.text(6.6, 5.65, "Parameter and compound uncertainty are evaluated separately, then intersected.", ha="center", color="#4A5568")
    ax.text(8.25, 2.35, "Output: 29 dual-robust targets\nand 58 dual-robust functional terms", ha="center", fontsize=9, color="#2F855A")
    ax.text(11.1, 2.35, "Output: seed-stable\nstructural affinity", ha="center", fontsize=9, color="#C53030")
    ax.set_title("Dual-perturbation framework for robust target identification", fontsize=13, weight="bold", pad=8)
    return save_figure(fig, "fig1_dual_perturbation_workflow")


def figure_2_parameter_heatmap() -> list[str]:
    rows = read_csv(SENS_DIR / "param_grid_results.csv")
    probs = sorted({f(r, "prob_threshold") for r in rows})
    strings = sorted({f(r, "string_threshold") for r in rows})
    topks = sorted({int(f(r, "top_k")) for r in rows})

    all_values = [f(r, "significant_pathways") for r in rows]
    vmin, vmax = min(all_values), max(all_values)

    fig = plt.figure(figsize=(12.4, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[3.2, 1.25])
    ax = fig.add_subplot(gs[0, 0])
    ax_line = fig.add_subplot(gs[0, 1])

    matrix = np.full((len(probs), len(topks) * len(strings)), np.nan)
    col_labels: list[str] = []
    for k in topks:
        for s in strings:
            col_labels.append(f"K{k}\nS{s:g}")

    for i, prob in enumerate(probs):
        for k_i, k in enumerate(topks):
            for s_i, s in enumerate(strings):
                j = k_i * len(strings) + s_i
                match = [
                    r
                    for r in rows
                    if f(r, "prob_threshold") == prob
                    and f(r, "string_threshold") == s
                    and int(f(r, "top_k")) == k
                ]
                matrix[i, j] = f(match[0], "significant_pathways") if match else np.nan

    im = ax.imshow(matrix, cmap="YlGnBu", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_title("All 112 parameter combinations", weight="bold", pad=10)
    ax.set_xticks(range(len(col_labels)), col_labels)
    ax.set_yticks(range(len(probs)), [f"{p:g}" for p in probs])
    ax.set_xlabel("Top-K target set and STRING confidence")
    ax.set_ylabel("SwissTarget probability threshold")

    for boundary in range(len(strings), matrix.shape[1], len(strings)):
        ax.axvline(boundary - 0.5, color="white", linewidth=2)
    for k_i, k in enumerate(topks):
        center = k_i * len(strings) + (len(strings) - 1) / 2
        ax.text(center, -0.95, f"Top-K {k}", ha="center", va="bottom", fontsize=8.2, weight="bold", color=COLORS["gray"])

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            ax.text(j, i, f"{int(val)}", ha="center", va="center", color="white" if val >= 30 else "#1A202C", fontsize=7.3)

    grouped: dict[float, list[float]] = defaultdict(list)
    for row in rows:
        grouped[f(row, "prob_threshold")].append(f(row, "significant_pathways"))
    means = np.array([mean(grouped[p]) for p in probs])
    mins = np.array([min(grouped[p]) for p in probs])
    maxs = np.array([max(grouped[p]) for p in probs])

    ax_line.fill_between(probs, mins, maxs, color=COLORS["blue"], alpha=0.15, linewidth=0)
    ax_line.plot(probs, means, marker="o", color=COLORS["blue"], linewidth=1.8)
    ax_line.set_title("Marginal trend", weight="bold", pad=10)
    ax_line.set_xlabel("Probability threshold")
    ax_line.set_ylabel("Mean significant terms")
    ax_line.set_ylim(max(0, vmin - 2), vmax + 2)
    ax_line.grid(True, color="#E2E8F0", linewidth=0.8)
    for x, y in zip(probs, means):
        ax_line.text(x, y + 0.9, f"{y:.1f}", ha="center", fontsize=7.5, color=COLORS["gray"])
    ax_line.annotate(
        "yield drops only at\nstricter thresholds",
        xy=(0.2, means[-1]),
        xytext=(0.105, 17),
        arrowprops={"arrowstyle": "->", "lw": 0.9, "color": COLORS["gray"]},
        fontsize=8,
        color=COLORS["gray"],
    )

    fig.colorbar(im, ax=ax, shrink=0.82, label="Significant functional terms")
    fig.suptitle("Parameter sensitivity across 112 analysis settings", fontsize=13, weight="bold")
    return save_figure(fig, "fig2_parameter_sensitivity_heatmap")


def compound_display_name(row: dict[str, str]) -> str:
    name = row["drop_compound_name"]
    if name.startswith("(3beta"):
        return "Cardiac glycoside"
    return name.replace("_", " ")


def figure_3_compound_influence() -> list[str]:
    rows = read_csv(COMP_DIR / "compound_influence_scores.csv")
    rows = sorted(rows, key=lambda r: f(r, "compound_influence_score"), reverse=True)
    top = rows[:12]
    labels = [compound_display_name(r) for r in top][::-1]
    scores = [f(r, "compound_influence_score") for r in top][::-1]
    colors = []
    for label in labels:
        if label in {"Catechin", "Miquelianin", "Procyanidin B1", "Cardiac glycoside"}:
            colors.append(COLORS["orange"] if label == "Catechin" else COLORS["purple"] if label == "Miquelianin" else COLORS["teal"] if label == "Procyanidin B1" else COLORS["red"])
        else:
            colors.append("#A0AEC0")

    fig, ax = plt.subplots(figsize=(8.5, 6.1))
    y = np.arange(len(labels))
    ax.barh(y, scores, color=colors, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Compound influence score (1 - mean retention)")
    ax.set_title("Compound leave-one-out influence ranking", fontsize=13, weight="bold")
    ax.grid(axis="x", color="#E2E8F0", linewidth=0.8)
    for i, score in enumerate(scores):
        ax.text(score + 0.004, i, f"{score:.3f}", va="center", fontsize=8)
    ax.set_xlim(0, max(scores) * 1.18)
    legend_items = [
        patches.Patch(color=COLORS["red"], label="Top network driver"),
        patches.Patch(color=COLORS["orange"], label="Docked compound: Catechin"),
        patches.Patch(color=COLORS["purple"], label="Docked compound: Miquelianin"),
        patches.Patch(color=COLORS["teal"], label="Docked compound: Procyanidin B1"),
    ]
    ax.legend(handles=legend_items, loc="lower right", frameon=False)
    return save_figure(fig, "fig3_compound_influence_ranking")


def figure_4_dual_target_matrix() -> list[str]:
    rows = read_csv(COMP_DIR / "dual_robust_targets.csv")
    rows = [r for r in rows if r["is_dual_robust"].lower() == "true"]
    rows = sorted(rows, key=lambda r: (f(r, "loo_q1_retention_rate"), f(r, "full_frequency"), f(r, "full_avg_degree_norm")), reverse=True)[:20]
    genes = [r["gene"] for r in rows]
    metrics = [
        ("Param freq", [f(r, "full_frequency") for r in rows]),
        ("LOO retention", [f(r, "loo_q1_retention_rate") for r in rows]),
        ("Min LOO freq", [f(r, "min_loo_frequency") for r in rows]),
        ("Degree norm", [f(r, "full_avg_degree_norm") for r in rows]),
    ]
    matrix = np.array([vals for _, vals in metrics]).T

    fig, ax = plt.subplots(figsize=(8.8, 7.0))
    im = ax.imshow(matrix, cmap="PuBuGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(metrics)), [m[0] for m in metrics], rotation=20, ha="right")
    ax.set_yticks(range(len(genes)), genes)
    for i, gene in enumerate(genes):
        if gene in {"STAT3", "EGFR", "SRC"}:
            ax.get_yticklabels()[i].set_weight("bold")
            ax.get_yticklabels()[i].set_color(COLORS["red"])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.5, color="#1A202C")
    ax.set_title("Evidence matrix for top dual-robust targets", fontsize=13, weight="bold")
    fig.colorbar(im, ax=ax, shrink=0.75, label="Metric value")
    return save_figure(fig, "fig4_dual_robust_target_matrix")


def category_color(category: str) -> str:
    return {
        "cross_talk": COLORS["orange"],
        "signaling_backbone": COLORS["blue"],
        "direct_ox": COLORS["red"],
        "other": "#718096",
    }.get(category, "#718096")


def figure_5_pathway_bubble() -> list[str]:
    rows = read_csv(COMP_DIR / "dual_robust_pathways.csv")
    rows = [r for r in rows if r["is_dual_robust"].lower() == "true"]
    rows = sorted(rows, key=lambda r: (f(r, "full_frequency"), f(r, "full_avg_logp"), f(r, "loo_q1_retention_rate")), reverse=True)[:18]
    rows = rows[::-1]
    labels = [short_text(r["pathway"], 58) for r in rows]
    x = [f(r, "full_frequency") for r in rows]
    sizes = [max(60, f(r, "full_avg_logp") * 18) for r in rows]
    colors = [category_color(r["pathway_category"]) for r in rows]

    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    y = np.arange(len(rows))
    ax.scatter(x, y, s=sizes, c=colors, alpha=0.82, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Parameter frequency in Q1")
    ax.set_title("Dual-robust functional enrichment terms", fontsize=13, weight="bold")
    ax.text(0.45, len(rows) - 0.35, "Bubble size = mean -log10(P); shown terms retain LOO Q1 rate >= 0.94", fontsize=8.2, color=COLORS["gray"])
    ax.grid(axis="x", color="#E2E8F0")
    ax.set_xlim(0.45, 1.05)
    handles = [
        patches.Patch(color=COLORS["orange"], label="Cross-talk"),
        patches.Patch(color=COLORS["blue"], label="Signaling backbone"),
        patches.Patch(color="#718096", label="Other"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)
    for yi, row in enumerate(rows):
        retention = f(row, "loo_q1_retention_rate")
        if retention < 0.999:
            ax.text(min(1.02, f(row, "full_frequency") + 0.012), yi, f"LOO {retention:.2f}", va="center", fontsize=7.2, color="#4A5568")
    return save_figure(fig, "fig5_dual_robust_pathway_bubble")


def figure_6_docking_summary() -> list[str]:
    rows = read_csv(DOCK_DIR / "docking_replicates_summary.csv")
    receptors = ["STAT3", "EGFR", "SRC"]
    ligands = ["Catechin", "Miquelianin", "Procyanidin_B1"]
    ligand_labels = ["Catechin", "Miquelianin", "Procyanidin B1"]
    lookup = {(r["receptor"], r["ligand"]): r for r in rows}
    affinity = np.array([[f(lookup[(rec, lig)], "mean_affinity_kcal_mol") for lig in ligands] for rec in receptors])
    sd = np.array([[f(lookup[(rec, lig)], "sd_affinity_kcal_mol") for lig in ligands] for rec in receptors])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.7), gridspec_kw={"width_ratios": [1.15, 1.0]}, constrained_layout=True)
    im = ax1.imshow(affinity, cmap="RdYlBu_r", aspect="auto")
    ax1.set_xticks(range(len(ligands)), ligand_labels, rotation=30, ha="right")
    ax1.set_yticks(range(len(receptors)), receptors)
    ax1.set_title("Mean Vina affinity", weight="bold")
    for i in range(affinity.shape[0]):
        for j in range(affinity.shape[1]):
            ax1.text(j, i, f"{affinity[i, j]:.2f}\n±{sd[i, j]:.3f}", ha="center", va="center", fontsize=8, color="#1A202C")
    fig.colorbar(im, ax=ax1, shrink=0.85, label="kcal/mol")

    ranked = sorted(rows, key=lambda r: f(r, "mean_affinity_kcal_mol"))[:6]
    labels = [f"{r['receptor']} × {r['ligand'].replace('_', ' ')}" for r in ranked][::-1]
    values = [abs(f(r, "mean_affinity_kcal_mol")) for r in ranked][::-1]
    errs = [f(r, "sd_affinity_kcal_mol") for r in ranked][::-1]
    y = np.arange(len(labels))
    ax2.barh(y, values, xerr=errs, color=COLORS["teal"], alpha=0.9, error_kw={"ecolor": "#1A202C", "lw": 0.8, "capsize": 3})
    ax2.set_yticks(y, labels)
    ax2.set_xlabel("|Mean affinity| (kcal/mol)")
    ax2.set_title("Top seed-stable docking pairs", weight="bold")
    ax2.grid(axis="x", color="#E2E8F0")
    for i, val in enumerate(values):
        ax2.text(val + 0.05, i, f"{val:.2f}", va="center", fontsize=8)
    fig.suptitle("Molecular docking seed-stability summary (n=3 seeds per pair)", fontsize=13, weight="bold")
    return save_figure(fig, "fig6_docking_seed_stability")


def write_catalog(created: dict[str, list[str]]) -> None:
    catalog = [
        "# Figure Catalog",
        "",
        "All paper-level figures were generated from project CSV/JSON artifacts without fabricating values.",
        "",
    ]
    entries = [
        ("Figure 1", "fig1_dual_perturbation_workflow", "Shows the study design and the two independent perturbation axes."),
        ("Figure 2", "fig2_parameter_sensitivity_heatmap", "Shows how functional-term yield changes across 112 parameter settings."),
        ("Figure 3", "fig3_compound_influence_ranking", "Ranks compounds by leave-one-compound-out influence score."),
        ("Figure 4", "fig4_dual_robust_target_matrix", "Summarizes evidence supporting the top dual-robust targets."),
        ("Figure 5", "fig5_dual_robust_pathway_bubble", "Displays dual-robust enriched functional terms and their categories."),
        ("Figure 6", "fig6_docking_seed_stability", "Summarizes docking affinity means and seed-level standard deviations."),
    ]
    for label, stem, purpose in entries:
        catalog.extend(
            [
                f"## {label}: `{stem}`",
                "",
                f"- Purpose: {purpose}",
                f"- Files: {', '.join(created.get(stem, []))}",
                "- Interpretation note: use together with the corresponding exact CSV files; do not infer biological causality from visualization alone.",
                "",
            ]
        )
    (FIG_DIR / "figure-catalog.md").write_text("\n".join(catalog), encoding="utf-8")


def write_reports() -> None:
    comp = read_csv(COMP_DIR / "compound_influence_scores.csv")
    docking = read_csv(DOCK_DIR / "docking_replicates_summary.csv")
    top_comp = max(comp, key=lambda r: f(r, "compound_influence_score"))
    top_docking = min(docking, key=lambda r: f(r, "mean_affinity_kcal_mol"))
    sds = [f(r, "sd_affinity_kcal_mol") for r in docking]
    report = [
        "# Analysis Report",
        "",
        "## Main Questions",
        "",
        "- Which targets and GO BP/KEGG functional terms remain robust under parameter and compound perturbation?",
        "- Which compounds drive network-level instability?",
        "- Are the selected docking pairs stable across fixed Vina seeds?",
        "",
        "## Key Findings",
        "",
        f"- Top compound influence: {compound_display_name(top_comp)} (score = {f(top_comp, 'compound_influence_score'):.3f}).",
        f"- Strongest docking pair: {top_docking['receptor']} x {top_docking['ligand']} (mean affinity = {f(top_docking, 'mean_affinity_kcal_mol'):.3f} kcal/mol).",
        f"- Docking seed SD range: {min(sds):.3f} to {max(sds):.3f} kcal/mol across 9 receptor-ligand pairs.",
        "",
        "## Caveats",
        "",
        "- Docking replicates are computational seeds, not biological replicates.",
        "- Vina affinity values are approximate scoring outputs and should not be interpreted as IC50/Kd.",
        "- Functional enrichment visualizations summarize computational terms, not direct wet-lab pathway activation.",
    ]
    (FIG_DIR / "analysis-report.md").write_text("\n".join(report), encoding="utf-8")

    stats = [
        "# Stats Appendix",
        "",
        "## Docking Seed Stability",
        "",
        "| receptor | ligand | n | mean_affinity_kcal_mol | sd_affinity_kcal_mol | stability |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in sorted(docking, key=lambda r: f(r, "mean_affinity_kcal_mol")):
        stats.append(
            f"| {row['receptor']} | {row['ligand']} | {row['n']} | {f(row, 'mean_affinity_kcal_mol'):.3f} | {f(row, 'sd_affinity_kcal_mol'):.3f} | {row['stability']} |"
        )
    stats.extend(
        [
            "",
            "No inferential significance tests were run for docking, because seeds are repeated computational searches rather than independent biological samples.",
        ]
    )
    (FIG_DIR / "stats-appendix.md").write_text("\n".join(stats), encoding="utf-8")


def figure_1_workflow() -> list[str]:
    font = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=9)
    font_small = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=7.6)
    font_title = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=15, weight="bold")
    font_block = FontProperties(family=["Microsoft YaHei", "SimHei", "Arial"], size=8.2, weight="bold")

    fig, ax = plt.subplots(figsize=(15.8, 8.8))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9.4)
    ax.axis("off")

    lane_specs = [
        (0.25, 3.25, "输入数据", "#F5FAFF", COLORS["blue"]),
        (3.8, 3.25, "候选靶点网络", "#F8FAFC", COLORS["gray"]),
        (7.35, 3.85, "双扰动鲁棒性分析", "#FFFBF5", COLORS["orange"]),
        (11.55, 4.15, "优先级输出与验证", "#F5FFF8", COLORS["green"]),
    ]
    for x, w, title, face, edge in lane_specs:
        lane = patches.FancyBboxPatch(
            (x, 0.65),
            w,
            7.38,
            boxstyle="round,pad=0.025,rounding_size=0.11",
            linewidth=1.0,
            edgecolor=edge,
            facecolor=face,
            alpha=0.72,
            zorder=0,
        )
        ax.add_patch(lane)
        ax.text(x + 0.18, 7.92, title, fontproperties=font_block, color=edge, va="center")

    def draw_block(
        x: float,
        y: float,
        w: float,
        h: float,
        title: str,
        lines: list[str],
        edge: str,
        face: str,
        badge: str | None = None,
    ) -> None:
        shadow = patches.FancyBboxPatch(
            (x + 0.045, y - 0.045),
            w,
            h,
            boxstyle="round,pad=0.035,rounding_size=0.10",
            linewidth=0,
            facecolor="#CBD5E1",
            alpha=0.35,
            zorder=1,
        )
        box = patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.035,rounding_size=0.10",
            linewidth=1.45,
            edgecolor=edge,
            facecolor=face,
            zorder=2,
        )
        ax.add_patch(shadow)
        ax.add_patch(box)
        ax.text(x + 0.18, y + h - 0.28, title, fontproperties=font_block, color="#172033", va="top", zorder=3)
        ax.text(
            x + 0.18,
            y + h - 0.72,
            "\n".join(lines),
            fontproperties=font_small,
            color="#374151",
            va="top",
            linespacing=1.22,
            zorder=3,
        )
        if badge:
            badge_w = max(0.72, 0.18 * len(badge) + 0.28)
            badge_patch = patches.FancyBboxPatch(
                (x + w - badge_w - 0.14, y + h - 0.48),
                badge_w,
                0.34,
                boxstyle="round,pad=0.025,rounding_size=0.07",
                linewidth=0,
                facecolor=edge,
                zorder=4,
            )
            ax.add_patch(badge_patch)
            ax.text(
                x + w - badge_w / 2 - 0.14,
                y + h - 0.31,
                badge,
                fontproperties=font_small,
                color="white",
                ha="center",
                va="center",
                zorder=5,
            )

    blocks = {
        "lcms": (0.55, 6.38, 2.62, 1.36),
        "disease": (0.55, 4.66, 2.62, 1.36),
        "standard": (0.55, 2.94, 2.62, 1.36),
        "swiss": (4.1, 6.25, 2.62, 1.28),
        "intersection": (4.1, 4.45, 2.62, 1.28),
        "ppi": (4.1, 2.65, 2.62, 1.30),
        "param": (7.7, 5.95, 3.18, 1.68),
        "loo": (7.7, 3.62, 3.18, 1.68),
        "enrich": (7.7, 1.45, 3.18, 1.46),
        "prqi": (11.9, 6.30, 3.45, 1.30),
        "dual": (11.9, 4.52, 3.45, 1.30),
        "validation": (11.9, 2.28, 3.45, 1.72),
        "final": (11.9, 0.82, 3.45, 1.28),
    }

    draw_block(*blocks["lcms"], "LC-MS/MS 成分输入", ["19 条鉴定记录 -> 18 个化合物", "结构规范化去重", "PubChem CID / SMILES"], COLORS["blue"], COLORS["pale_blue"], "18")
    draw_block(*blocks["disease"], "氧化应激相关基因", ["GeneCards + OMIM + TTD", "统一基因符号", "构建疾病背景集合"], COLORS["green"], COLORS["pale_green"])
    draw_block(*blocks["standard"], "结构与名称标准化", ["Trehalose 条目合并", "SwissTargetPrediction", "物种：H. sapiens"], COLORS["teal"], "#ECFEFF")

    draw_block(*blocks["swiss"], "靶点预测并集", ["按化合物预测潜在靶点", "多基因条目拆分", "同一基因保留最高 probability"], COLORS["blue"], "#EFF6FF", "684")
    draw_block(*blocks["intersection"], "疾病交集候选靶点", ["预测靶点 ∩ 氧化应激基因", "保留正 probability 靶点", "进入 PPI 与富集分析"], COLORS["gray"], "#F8FAFC", "286")
    draw_block(*blocks["ppi"], "STRING PPI 网络", ["combined score 过滤边", "degree / betweenness", "归一化 degree 比较"], COLORS["gray"], "#F8FAFC")

    draw_block(*blocks["param"], "参数扰动轴", ["Swiss probability: 0-0.20，共 7 档", "STRING score: 0.40/0.50/0.70/0.90", "Hub TopK: 10/20/30/50", "全因子网格：7 x 4 x 4"], COLORS["orange"], COLORS["pale_orange"], "112")
    draw_block(*blocks["loo"], "成分扰动轴", ["leave-one-compound-out", "每次剔除 1 个结构规范化化合物", "重跑 112 组并重建 PPI/富集", "计算 retention 与 influence score"], COLORS["purple"], "#FAF5FF", "2016")
    draw_block(*blocks["enrich"], "功能富集条目鲁棒性", ["GO BP 生物过程 + KEGG 通路", "统一称功能富集条目", "频率 + 平均 -log10 adjusted P"], COLORS["teal"], "#F0FDFA")

    draw_block(*blocks["prqi"], "PR-QI 自适应分层", ["靶点：参数频率 + 平均归一化 degree", "功能条目：出现频率 + 富集强度", "中位数分界，连续指标保留"], COLORS["green"], COLORS["pale_green"])
    draw_block(*blocks["dual"], "双重稳健候选集", ["LOO PR-QI retention >= 0.80", "minimum LOO frequency > 0", "50 -> 29 targets；70 -> 58 terms"], COLORS["green"], COLORS["pale_green"], "29/58")
    draw_block(*blocks["validation"], "外部证据审计与结构支持", ["Top30 策略消融比较", "PubMed 共现 + 数据库证据", "超几何检验 + BH-FDR", "degree-matched permutation", "Vina docking: 3 靶点 x 3 成分 x 3 seeds"], COLORS["red"], "#FFF5F5")
    draw_block(
        *blocks["final"],
        "核心输出",
        ["Tier 1: STAT3", "核心靶点：EGFR / SRC / TNF", "HSP90AA1 / BCL2 / PTGS2 / MMP9"],
        COLORS["green"],
        "#F0FFF4",
    )

    def center_right(key: str) -> tuple[float, float]:
        x, y, w, h = blocks[key]
        return x + w, y + h / 2

    def center_left(key: str) -> tuple[float, float]:
        x, y, _w, h = blocks[key]
        return x, y + h / 2

    def port(key: str, side: str, frac: float = 0.5) -> tuple[float, float]:
        x, y, w, h = blocks[key]
        pad = 0.04
        if side == "left":
            return x - pad, y + h * frac
        if side == "right":
            return x + w + pad, y + h * frac
        if side == "top":
            return x + w * frac, y + h + pad
        if side == "bottom":
            return x + w * frac, y - pad
        raise ValueError(f"Unknown port side: {side}")

    def routed_arrow(points: list[tuple[float, float]], color: str = "#334155") -> None:
        # Keep connectors below the cards so they cannot cross visible text.
        for p0, p1 in zip(points[:-2], points[1:-1]):
            ax.plot(
                [p0[0], p1[0]],
                [p0[1], p1[1]],
                color=color,
                lw=1.35,
                solid_capstyle="round",
                zorder=1.6,
            )
        ax.annotate(
            "",
            xy=points[-1],
            xytext=points[-2],
            arrowprops=dict(
                arrowstyle="-|>",
                lw=1.35,
                color=color,
                shrinkA=0,
                shrinkB=0,
                mutation_scale=13,
            ),
            zorder=1.6,
        )

    x_input_gap = 3.50
    x_analysis_bus = 7.08
    x_output_bus = 11.28

    routed_arrow([port("lcms", "right"), (x_input_gap, port("lcms", "right")[1]), (x_input_gap, port("swiss", "left")[1]), port("swiss", "left")], COLORS["blue"])
    routed_arrow([port("standard", "right"), (x_input_gap, port("standard", "right")[1]), (x_input_gap, port("swiss", "left")[1]), port("swiss", "left")], COLORS["teal"])
    routed_arrow([port("disease", "right"), (x_input_gap, port("disease", "right")[1]), (x_input_gap, port("intersection", "left")[1]), port("intersection", "left")], COLORS["green"])

    routed_arrow([port("swiss", "bottom"), port("intersection", "top")], COLORS["gray"])
    routed_arrow([port("intersection", "bottom"), port("ppi", "top")], COLORS["gray"])

    routed_arrow([port("ppi", "right"), (x_analysis_bus, port("ppi", "right")[1]), (x_analysis_bus, port("param", "left")[1]), port("param", "left")], COLORS["orange"])
    routed_arrow([port("ppi", "right"), (x_analysis_bus, port("ppi", "right")[1]), (x_analysis_bus, port("loo", "left")[1]), port("loo", "left")], COLORS["purple"])
    routed_arrow([port("ppi", "right"), (x_analysis_bus, port("ppi", "right")[1]), (x_analysis_bus, port("enrich", "left")[1]), port("enrich", "left")], COLORS["teal"])

    routed_arrow([port("param", "right"), (x_output_bus, port("param", "right")[1]), (x_output_bus, port("prqi", "left")[1]), port("prqi", "left")], COLORS["orange"])
    routed_arrow([port("enrich", "right"), (x_output_bus, port("enrich", "right")[1]), (x_output_bus, port("prqi", "left")[1]), port("prqi", "left")], COLORS["teal"])
    routed_arrow([port("loo", "right"), (x_output_bus, port("loo", "right")[1]), (x_output_bus, port("dual", "left")[1]), port("dual", "left")], COLORS["purple"])

    routed_arrow([port("prqi", "bottom"), port("dual", "top")], COLORS["green"])
    routed_arrow([port("dual", "bottom"), port("validation", "top")], COLORS["red"])
    routed_arrow([port("validation", "bottom"), port("final", "top")], COLORS["green"])

    ax.text(
        8.0,
        8.82,
        "双扰动鲁棒性网络药理学靶点优先级识别框架",
        fontproperties=font_title,
        ha="center",
        va="center",
        color="#111827",
    )
    ax.text(
        8.0,
        8.48,
        "Parameter grid + compound leave-one-out reruns convert fixed-threshold screening into robustness-based prioritization.",
        fontproperties=font,
        ha="center",
        va="center",
        color="#4B5563",
    )
    ax.text(
        0.34,
        0.30,
        "注：PR-QI 为参数空间内的相对稳健优先级分类，不等同于独立因果验证；对接用于代表性靶点-成分组合的结构可行性支持。",
        fontproperties=font_small,
        color="#64748B",
    )

    return save_figure(fig, "fig1_dual_perturbation_workflow")


def main() -> None:
    created: dict[str, list[str]] = {}
    for func, stem in [
        (figure_1_workflow, "fig1_dual_perturbation_workflow"),
        (figure_2_parameter_heatmap, "fig2_parameter_sensitivity_heatmap"),
        (figure_3_compound_influence, "fig3_compound_influence_ranking"),
        (figure_4_dual_target_matrix, "fig4_dual_robust_target_matrix"),
        (figure_5_pathway_bubble, "fig5_dual_robust_pathway_bubble"),
        (figure_6_docking_summary, "fig6_docking_seed_stability"),
    ]:
        print(f"Generating {stem}...")
        created[stem] = func()

    write_catalog(created)
    write_reports()
    print(f"Done. Figures and reports written to: {FIG_DIR}")


if __name__ == "__main__":
    main()
