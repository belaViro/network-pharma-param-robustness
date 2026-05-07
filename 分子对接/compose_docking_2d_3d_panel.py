from __future__ import annotations

import csv
import shutil
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
RESULTS_DIR = ROOT / "results"
PAPER_DIR = ROOT.parent / "src" / "outputs" / "paper_figures"

RECEPTORS = ["STAT3", "EGFR", "SRC"]
LIGANDS = ["Catechin", "Miquelianin", "Procyanidin_B1"]
LIGAND_LABELS = {
    "Catechin": "Catechin",
    "Miquelianin": "Miquelianin",
    "Procyanidin_B1": "Procyanidin B1",
}

OUT_BASENAME = "fig8_docking_2d_3d_binding_modes"


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        Path("C:/Windows/Fonts/msyhbd.ttc" if bold else "C:/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/simhei.ttf"),
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
    ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=font)
    return box[2] - box[0], box[3] - box[1]


def nonwhite_bbox(img: Image.Image, focus: tuple[int, int, int, int], threshold: int = 244) -> tuple[int, int, int, int]:
    cropped = img.crop(focus).convert("RGB")
    arr = np.asarray(cropped)
    mask = np.any(arr < threshold, axis=2)
    ys, xs = np.where(mask)
    if len(xs) == 0 or len(ys) == 0:
        return focus

    left = int(xs.min()) + focus[0]
    upper = int(ys.min()) + focus[1]
    right = int(xs.max()) + focus[0] + 1
    lower = int(ys.max()) + focus[1] + 1
    return left, upper, right, lower


def expand_box(box: tuple[int, int, int, int], img_size: tuple[int, int], pad: int) -> tuple[int, int, int, int]:
    left, upper, right, lower = box
    width, height = img_size
    return (
        max(0, left - pad),
        max(0, upper - pad),
        min(width, right + pad),
        min(height, lower + pad),
    )


def crop_panel_image(path: Path, kind: str) -> Image.Image:
    img = Image.open(path).convert("RGB")
    width, height = img.size

    if kind == "3d":
        focus = (90, 95, width - 120, height - 120)
        pad = 65
    else:
        focus = (35, 55, width - 35, height - 55)
        pad = 35

    box = expand_box(nonwhite_bbox(img, focus), img.size, pad)
    return img.crop(box)


def fit_image(img: Image.Image, box_size: tuple[int, int]) -> Image.Image:
    max_w, max_h = box_size
    scale = min(max_w / img.width, max_h / img.height)
    new_size = (max(1, int(img.width * scale)), max(1, int(img.height * scale)))
    return img.resize(new_size, Image.Resampling.LANCZOS)


def read_energy_table() -> dict[tuple[str, str], tuple[float, float]]:
    summary_path = RESULTS_DIR / "validation_runs" / "docking_replicates_summary.csv"
    energies: dict[tuple[str, str], tuple[float, float]] = {}
    if not summary_path.exists():
        return energies

    with summary_path.open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            receptor = row["receptor"]
            ligand = row["ligand"]
            energies[(receptor, ligand)] = (
                float(row["mean_affinity_kcal_mol"]),
                float(row["sd_affinity_kcal_mol"]),
            )
    return energies


def draw_centered_text(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    font: ImageFont.ImageFont,
    fill: str,
) -> None:
    x1, y1, x2, y2 = box
    w, h = text_size(draw, text, font)
    draw.text((x1 + (x2 - x1 - w) / 2, y1 + (y2 - y1 - h) / 2), text, font=font, fill=fill)


def draw_legend(draw: ImageDraw.ImageDraw, x: int, y: int, font: ImageFont.ImageFont) -> None:
    items = [
        ("H-bond", "#2CA25F"),
        ("Hydrophobic", "#B45ACB"),
        ("Pi-Pi", "#42B7B9"),
    ]
    cursor = x
    for label, color in items:
        draw.rounded_rectangle((cursor, y, cursor + 44, y + 18), radius=7, fill=color, outline="#777777")
        draw.text((cursor + 52, y - 2), label, font=font, fill="#333333")
        label_w, _ = text_size(draw, label, font)
        cursor += 52 + label_w + 34


def compose() -> list[Path]:
    energies = read_energy_table()

    title_font = load_font(50, bold=True)
    header_font = load_font(31, bold=True)
    small_font = load_font(23)
    label_font = load_font(25, bold=True)
    tiny_font = load_font(20)

    margin_x = 80
    margin_top = 118
    margin_bottom = 86
    gap = 22
    cell_w = 1660
    cell_h = 1050
    header_h = 96
    panel_label_h = 34
    pane_gap = 16

    canvas_w = margin_x * 2 + cell_w * len(LIGANDS) + gap * (len(LIGANDS) - 1)
    canvas_h = margin_top + cell_h * len(RECEPTORS) + gap * (len(RECEPTORS) - 1) + margin_bottom
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    draw = ImageDraw.Draw(canvas)

    title = "分子对接结合构象：3D 姿态与 2D 相互作用"
    title_w, _ = text_size(draw, title, title_font)
    draw.text(((canvas_w - title_w) / 2, 34), title, font=title_font, fill="#111111")

    panel_idx = 0
    for row_idx, receptor in enumerate(RECEPTORS):
        for col_idx, ligand in enumerate(LIGANDS):
            panel_idx += 1
            x = margin_x + col_idx * (cell_w + gap)
            y = margin_top + row_idx * (cell_h + gap)
            x2 = x + cell_w
            y2 = y + cell_h

            draw.rounded_rectangle((x, y, x2, y2), radius=10, fill="#FFFFFF", outline="#D0D4DA", width=2)
            draw.rectangle((x, y, x2, y + header_h), fill="#F3F6F9")

            ligand_label = LIGAND_LABELS[ligand]
            mean_sd = energies.get((receptor, ligand))
            if mean_sd:
                mean, sd = mean_sd
                subtitle = f"{receptor} x {ligand_label}    Vina score = {mean:.2f} +/- {sd:.2f} kcal/mol"
            else:
                subtitle = f"{receptor} x {ligand_label}"

            letter = chr(ord("A") + panel_idx - 1)
            draw.text((x + 24, y + 27), letter, font=header_font, fill="#111111")
            draw.text((x + 76, y + 33), subtitle, font=small_font, fill="#263238")

            content_y = y + header_h + 18
            pane_w = (cell_w - 60 - pane_gap) // 2
            pane_h = cell_h - header_h - 52
            pane1 = (x + 24, content_y, x + 24 + pane_w, content_y + pane_h)
            pane2 = (x + 24 + pane_w + pane_gap, content_y, x + 24 + pane_w + pane_gap + pane_w, content_y + pane_h)

            for pane, pane_label in [(pane1, "3D pose"), (pane2, "2D interaction")]:
                px1, py1, px2, py2 = pane
                draw.rectangle((px1, py1, px2, py2), fill="#FFFFFF", outline="#E1E5EA", width=1)
                draw_centered_text(draw, (px1, py1, px2, py1 + panel_label_h), pane_label, label_font, "#334155")

            paths = [
                (RESULTS_DIR / f"{receptor}_{ligand}_3d.png", "3d", pane1),
                (RESULTS_DIR / f"{receptor}_{ligand}_2d_interaction.png", "2d", pane2),
            ]

            for path, kind, pane in paths:
                if not path.exists():
                    continue
                img = crop_panel_image(path, kind)
                px1, py1, px2, py2 = pane
                fitted = fit_image(img, (px2 - px1 - 34, py2 - py1 - panel_label_h - 28))
                paste_x = px1 + (px2 - px1 - fitted.width) // 2
                paste_y = py1 + panel_label_h + 10 + (py2 - py1 - panel_label_h - 18 - fitted.height) // 2
                canvas.paste(fitted, (paste_x, paste_y))

    footer = "Panels show the best docking pose for each receptor-ligand pair; energy labels summarize three seed replicates."
    draw.text((margin_x, canvas_h - 54), footer, font=tiny_font, fill="#4B5563")
    draw_legend(draw, canvas_w - 760, canvas_h - 54, tiny_font)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PAPER_DIR.mkdir(parents=True, exist_ok=True)

    result_png = RESULTS_DIR / f"{OUT_BASENAME}.png"
    result_pdf = RESULTS_DIR / f"{OUT_BASENAME}.pdf"
    paper_png = PAPER_DIR / f"{OUT_BASENAME}.png"
    paper_pdf = PAPER_DIR / f"{OUT_BASENAME}.pdf"

    canvas.save(result_png, dpi=(300, 300))
    canvas.save(result_pdf, "PDF", resolution=300.0)
    for src, dst in [(result_png, paper_png), (result_pdf, paper_pdf)]:
        try:
            shutil.copyfile(src, dst)
        except PermissionError:
            print(f"Warning: could not copy to paper figure directory: {dst}")

    outputs = [result_png, result_pdf, paper_png, paper_pdf]
    return outputs


def main() -> None:
    for path in compose():
        print(path)


if __name__ == "__main__":
    main()
