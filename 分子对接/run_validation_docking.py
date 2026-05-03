"""
Run seed-stability validation for the molecular docking module.

This script performs 9 core receptor-ligand dockings with three fixed seeds
and writes CSV, JSON, and Markdown summaries under results/validation_runs/.
"""
from __future__ import annotations

import csv
import json
import re
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
RECEPTORS_DIR = BASE_DIR / "receptors"
LIGANDS_DIR = BASE_DIR / "ligands"
RESULTS_DIR = BASE_DIR / "results"
VALIDATION_DIR = RESULTS_DIR / "validation_runs"
POSES_DIR = VALIDATION_DIR / "poses"
LOGS_DIR = VALIDATION_DIR / "logs"

VINA_EXE = BASE_DIR / "vina_1.2.7_win.exe"

SEEDS = [1001, 1002, 1003]
EXHAUSTIVENESS = 16
CPU = 4
NUM_MODES = 9
BOX_SIZE = (20, 20, 20)


@dataclass(frozen=True)
class Receptor:
    name: str
    pdb_id: str
    domain: str
    center: tuple[float, float, float]
    pdbqt: Path


@dataclass(frozen=True)
class Ligand:
    name: str
    pdbqt: Path
    kind: str = "core"


RECEPTORS = {
    "STAT3": Receptor(
        name="STAT3",
        pdb_id="6NJS",
        domain="SH2",
        center=(13.028, 55.637, 0.250),
        pdbqt=RECEPTORS_DIR / "STAT3.pdbqt",
    ),
    "EGFR": Receptor(
        name="EGFR",
        pdb_id="1M17",
        domain="kinase",
        center=(22.014, 0.253, 52.794),
        pdbqt=RECEPTORS_DIR / "EGFR.pdbqt",
    ),
    "SRC": Receptor(
        name="SRC",
        pdb_id="2H8H",
        domain="kinase",
        center=(20.255, 20.924, 58.253),
        pdbqt=RECEPTORS_DIR / "SRC.pdbqt",
    ),
}

CORE_LIGANDS = {
    "Catechin": LIGANDS_DIR / "Catechin.pdbqt",
    "Miquelianin": LIGANDS_DIR / "Miquelianin.pdbqt",
    "Procyanidin_B1": LIGANDS_DIR / "Procyanidin_B1.pdbqt",
}


def ensure_inputs() -> None:
    required = [VINA_EXE]
    required.extend(receptor.pdbqt for receptor in RECEPTORS.values())
    required.extend(CORE_LIGANDS.values())
    missing = [path for path in required if not path.exists()]
    if missing:
        formatted = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(f"Required docking input is missing:\n{formatted}")


def parse_best_affinity(stdout: str) -> float:
    for line in stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        match = re.match(r"^1\s+(-?\d+(?:\.\d+)?)", stripped)
        if match:
            return float(match.group(1))
    raise ValueError("Could not parse Vina best affinity from stdout")


def run_vina(receptor: Receptor, ligand: Ligand, seed: int) -> dict[str, object]:
    output_pdbqt = POSES_DIR / f"{receptor.name}_{ligand.name}_seed{seed}.pdbqt"
    log_file = LOGS_DIR / f"{receptor.name}_{ligand.name}_seed{seed}.log"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(VINA_EXE),
        "--receptor",
        str(receptor.pdbqt),
        "--ligand",
        str(ligand.pdbqt),
        "--center_x",
        str(receptor.center[0]),
        "--center_y",
        str(receptor.center[1]),
        "--center_z",
        str(receptor.center[2]),
        "--size_x",
        str(BOX_SIZE[0]),
        "--size_y",
        str(BOX_SIZE[1]),
        "--size_z",
        str(BOX_SIZE[2]),
        "--out",
        str(output_pdbqt),
        "--cpu",
        str(CPU),
        "--num_modes",
        str(NUM_MODES),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--seed",
        str(seed),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    log_file.write_text(
        "COMMAND:\n"
        + " ".join(cmd)
        + "\n\nSTDOUT:\n"
        + result.stdout
        + "\n\nSTDERR:\n"
        + result.stderr,
        encoding="utf-8",
    )
    if result.returncode != 0:
        raise RuntimeError(f"Vina failed for {receptor.name} x {ligand.name}; see {log_file}")

    return {
        "receptor": receptor.name,
        "ligand": ligand.name,
        "kind": ligand.kind,
        "seed": seed,
        "affinity_kcal_mol": parse_best_affinity(result.stdout),
        "pose_file": str(output_pdbqt.relative_to(BASE_DIR)),
        "log_file": str(log_file.relative_to(BASE_DIR)),
    }


def summarize_replicates(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[float]] = {}
    for row in rows:
        key = (str(row["receptor"]), str(row["ligand"]), str(row["kind"]))
        grouped.setdefault(key, []).append(float(row["affinity_kcal_mol"]))

    summaries: list[dict[str, object]] = []
    for (receptor, ligand, kind), values in sorted(grouped.items()):
        sd = statistics.stdev(values) if len(values) > 1 else 0.0
        if sd < 0.3:
            stability = "stable"
        elif sd <= 0.6:
            stability = "acceptable"
        else:
            stability = "unstable"
        summaries.append(
            {
                "receptor": receptor,
                "ligand": ligand,
                "kind": kind,
                "n": len(values),
                "mean_affinity_kcal_mol": round(statistics.mean(values), 3),
                "sd_affinity_kcal_mol": round(sd, 3),
                "best_affinity_kcal_mol": round(min(values), 3),
                "worst_affinity_kcal_mol": round(max(values), 3),
                "stability": stability,
            }
        )
    return summaries


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    def format_value(value: object) -> str:
        if isinstance(value, float):
            return f"{value:.3f}"
        return str(value)

    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = [
        "| " + " | ".join(format_value(row.get(column, "")) for column in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, divider, *body])


def write_markdown_report(summaries: list[dict[str, object]], output_path: Path) -> None:
    ranked_rows = sorted(summaries, key=lambda row: float(row["mean_affinity_kcal_mol"]))
    lines = [
        "# 分子对接 seed 稳定性验证报告",
        "",
        "## 方法概览",
        "",
        "- AutoDock Vina: v1.2.7",
        f"- exhaustiveness: {EXHAUSTIVENESS}",
        f"- seeds: {', '.join(str(seed) for seed in SEEDS)}",
        f"- box size: {BOX_SIZE[0]} x {BOX_SIZE[1]} x {BOX_SIZE[2]} A",
        "- core combinations: 3 receptors x 3 ligands = 9",
        "- replicate dockings: 9 combinations x 3 seeds = 27",
        "",
        "## 核心化合物重复对接汇总",
        "",
        markdown_table(
            ranked_rows,
            [
                "receptor",
                "ligand",
                "n",
                "mean_affinity_kcal_mol",
                "sd_affinity_kcal_mol",
                "best_affinity_kcal_mol",
                "stability",
            ],
        ),
        "",
        "## 判读标准",
        "",
        "- sd < 0.3 kcal/mol: stable",
        "- 0.3 <= sd <= 0.6 kcal/mol: acceptable",
        "- sd > 0.6 kcal/mol: unstable",
        "",
        "## 局限",
        "",
        "- Vina scoring is approximate and cannot replace enzyme or cell experiments.",
        "- Protein receptors are treated as rigid structures.",
        "- Solvent and water-mediated interactions are not explicitly modeled.",
        "- This report evaluates seed-level sampling stability only.",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_inputs()
    for directory in [VALIDATION_DIR, POSES_DIR, LOGS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    ligands = [Ligand(name=name, pdbqt=path) for name, path in CORE_LIGANDS.items()]
    rows: list[dict[str, object]] = []

    for receptor in RECEPTORS.values():
        for ligand in ligands:
            for seed in SEEDS:
                print(f"Docking {receptor.name} x {ligand.name} seed={seed}")
                rows.append(run_vina(receptor, ligand, seed))

    summaries = summarize_replicates(rows)

    write_csv(
        VALIDATION_DIR / "docking_replicates.csv",
        rows,
        ["receptor", "ligand", "kind", "seed", "affinity_kcal_mol", "pose_file", "log_file"],
    )
    write_csv(
        VALIDATION_DIR / "docking_replicates_summary.csv",
        summaries,
        [
            "receptor",
            "ligand",
            "kind",
            "n",
            "mean_affinity_kcal_mol",
            "sd_affinity_kcal_mol",
            "best_affinity_kcal_mol",
            "worst_affinity_kcal_mol",
            "stability",
        ],
    )

    metadata = {
        "tool": "AutoDock Vina v1.2.7",
        "settings": {
            "exhaustiveness": EXHAUSTIVENESS,
            "cpu": CPU,
            "num_modes": NUM_MODES,
            "box_size_angstrom": BOX_SIZE,
            "seeds": SEEDS,
        },
        "receptors": {
            name: {
                "pdb_id": receptor.pdb_id,
                "domain": receptor.domain,
                "center": receptor.center,
            }
            for name, receptor in RECEPTORS.items()
        },
        "replicates": rows,
        "summary": summaries,
    }
    (VALIDATION_DIR / "validation_summary.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_markdown_report(summaries, VALIDATION_DIR / "validation_report.md")

    print(f"Validation outputs written to: {VALIDATION_DIR}")


if __name__ == "__main__":
    main()
