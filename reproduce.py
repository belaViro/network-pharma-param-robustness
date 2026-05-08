#!/usr/bin/env python3
"""One-command reproduction orchestrator for the manuscript workflow.

Default behavior runs the full analysis workflow from packaged repository inputs:
traditional baseline, parameter perturbation, compound leave-one-out perturbation,
external evidence audit, final manuscript figures, and supplementary tables.

Use ``--core`` for the faster parameter-robustness workflow only.
Use ``--include-acquisition`` only when fresh PubChem/SwissTarget acquisition is
needed; it requires network access and the optional acquisition dependencies.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "src" / "outputs"
PYTHON = sys.executable


@dataclass(frozen=True)
class Stage:
    name: str
    command: list[str]


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the network pharmacology dual-perturbation workflow. "
            "The default mode is the full manuscript workflow from packaged inputs."
        )
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--full",
        dest="mode",
        action="store_const",
        const="full",
        default="full",
        help="run the full manuscript workflow (default)",
    )
    mode.add_argument(
        "--core",
        dest="mode",
        action="store_const",
        const="core",
        help="run only the core parameter-robustness workflow",
    )
    parser.add_argument("--skip-install", action="store_true", help="do not run pip install -e .")
    parser.add_argument(
        "--include-acquisition",
        action="store_true",
        help="also rerun PubChem/SwissTarget acquisition before downstream analyses",
    )
    parser.add_argument(
        "--skip-traditional",
        action="store_true",
        help="in full mode, skip the traditional fixed-threshold baseline branch",
    )
    parser.add_argument(
        "--compound-force",
        action="store_true",
        help="force rerun completed compound LOO checkpoints",
    )
    parser.add_argument(
        "--pubmed-mode",
        choices=["skip", "cache", "query"],
        default="cache",
        help="PubMed handling for robustness validation; cache mode avoids network calls",
    )
    parser.add_argument(
        "--pubmed-scope",
        choices=["target_sets", "background", "universe"],
        default="target_sets",
        help="genes queried when --pubmed-mode=query",
    )
    parser.add_argument(
        "--n-permutations",
        type=positive_int,
        default=1000,
        help="degree-matched support permutations in the validation bundle",
    )
    parser.add_argument("--dry-run", action="store_true", help="print planned commands without running them")
    return parser.parse_args()


def command_text(command: Iterable[str]) -> str:
    return " ".join(str(item) for item in command)


def run_stage(index: int, total: int, stage: Stage, dry_run: bool) -> dict[str, object]:
    print()
    print(f"[{index}/{total}] {stage.name}")
    print(f"    {command_text(stage.command)}")

    started = time.perf_counter()
    if dry_run:
        return {
            "stage": stage.name,
            "command": stage.command,
            "status": "dry_run",
            "duration_seconds": 0.0,
        }

    subprocess.run(stage.command, cwd=ROOT, check=True)
    duration = round(time.perf_counter() - started, 3)
    print(f"    completed in {duration}s")
    return {
        "stage": stage.name,
        "command": stage.command,
        "status": "completed",
        "duration_seconds": duration,
    }


def install_stage(include_acquisition: bool) -> Stage:
    package = ".[acquisition]" if include_acquisition else "."
    return Stage("Install project dependencies", [PYTHON, "-m", "pip", "install", "-e", package, "--quiet"])


def acquisition_stages() -> list[Stage]:
    return [
        Stage(
            "Fresh PubChem and SwissTargetPrediction acquisition",
            [PYTHON, "scripts/run_pipeline.py"],
        )
    ]


def core_stages() -> list[Stage]:
    return [
        Stage("Parameter grid sensitivity analysis", [PYTHON, "scripts/analysis/run_sensitivity_analysis.py"]),
        Stage("Robust-vs-parametric enrichment comparison", [PYTHON, "scripts/analysis/compare_robust_vs_parametric.py"]),
        Stage("GO BP/KEGG functional-term quadrant augmentation", [PYTHON, "scripts/analysis/augment_pathway_robustness.py"]),
        Stage("Sensitivity figure regeneration", [PYTHON, "scripts/analysis/plot_sensitivity_figures.py"]),
        Stage("Functional-term quadrant figure", [PYTHON, "scripts/analysis/plot_pathway_quadrant.py"]),
        Stage("Functional-term quadrant figure by source", [PYTHON, "scripts/analysis/plot_pathway_quadrant_by_source.py"]),
        Stage("Parameter-space stability figure", [PYTHON, "scripts/analysis/plot_param_space_fig10.py"]),
        Stage("Top-set overlap figure", [PYTHON, "scripts/analysis/plot_upset_fig8.py"]),
        Stage("Network robustness figure", [PYTHON, "scripts/analysis/plot_network_fig9.py"]),
    ]


def traditional_stage() -> Stage:
    return Stage(
        "Traditional fixed-threshold baseline branch",
        [PYTHON, "scripts/run_traditional_pipeline.py"],
    )


def post_core_full_stages(args: argparse.Namespace) -> list[Stage]:
    stages: list[Stage] = []

    compound_command = [PYTHON, "scripts/analysis/run_compound_sensitivity_analysis.py"]
    if args.compound_force:
        compound_command.append("--force")
    stages.append(Stage("Compound leave-one-out perturbation analysis", compound_command))

    stages.append(
        Stage(
            "External evidence, Top30 ablation, and degree-matched validation bundle",
            [
                PYTHON,
                "scripts/analysis/robustness_validation_bundle.py",
                "--pubmed-mode",
                args.pubmed_mode,
                "--pubmed-scope",
                args.pubmed_scope,
                "--n-permutations",
                str(args.n_permutations),
            ],
        )
    )
    stages.append(Stage("Final manuscript figure generation", [PYTHON, "scripts/analysis/generate_paper_figures.py"]))
    stages.append(Stage("Supplementary table export", [PYTHON, "scripts/analysis/export_supplementary_tables.py"]))
    return stages


def build_stages(args: argparse.Namespace) -> list[Stage]:
    stages: list[Stage] = []
    if not args.skip_install:
        stages.append(install_stage(args.include_acquisition))
    if args.include_acquisition:
        stages.extend(acquisition_stages())
    if args.mode == "full":
        if not args.skip_traditional:
            stages.append(traditional_stage())
    stages.extend(core_stages())
    if args.mode == "full":
        stages.extend(post_core_full_stages(args))
    return stages


def write_manifest(args: argparse.Namespace, results: list[dict[str, object]], started_at: str) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTPUT_DIR / "reproduction_manifest.json"
    payload = {
        "workflow": "full_manuscript" if args.mode == "full" else "core_parameter_robustness",
        "started_at": started_at,
        "finished_at": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version,
        "arguments": vars(args),
        "stages": results,
    }
    manifest_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def main() -> int:
    if sys.version_info < (3, 10):
        print("Error: Python >= 3.10 is required.", file=sys.stderr)
        return 1

    args = parse_args()
    started_at = datetime.now().isoformat(timespec="seconds")
    stages = build_stages(args)

    print("Network pharmacology dual-perturbation reproduction")
    print(f"Mode: {args.mode}")
    print(f"Project root: {ROOT}")
    print(f"Stages: {len(stages)}")
    if args.mode == "full" and args.pubmed_mode == "cache":
        print("PubMed mode: cache (no network calls for PubMed counts)")

    results: list[dict[str, object]] = []
    for index, stage in enumerate(stages, start=1):
        results.append(run_stage(index, len(stages), stage, args.dry_run))

    manifest_path = write_manifest(args, results, started_at)
    print()
    print("Reproduction workflow completed." if not args.dry_run else "Dry run completed.")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
