#!/usr/bin/env python3
"""
Full Pipeline - 完整流水线 (Stage 1-6 + 传统方法评估)

Stage 1: PubChem 查询 SMILES
Stage 2: SwissTarget 靶点预测
Stage 3: 合并去重靶点
Stage 4: 疾病基因交集
Stage 5: 交集分析
Stage 6: PPI 网络构建
Stage 7: 传统方法评估 (Degree Top-30 → 富集分析 → 命中率评估)

Usage:
    python scripts/run_full_pipeline.py
"""

import sys
import datetime
from pathlib import Path
from typing import List, Dict

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.pipeline import (
    stage1_fetch_all,
    stage2_predict_all,
    merge_and_deduplicate_targets,
    merge_disease_genes,
    intersect_targets_with_disease_genes,
    build_ppi_network,
)
from src.config import STRAWBERRY_LEAF_MOLECULES
from src.enrichment import run_enrichment_analysis, plot_all_enrichment
from scripts.evaluate_traditional import evaluate_traditional_method
from src.logging_config import setup_logging
import logging

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent


def run_full_pipeline(
    molecule_names: List[str],
    max_workers: int = 4,
    run_traditional: bool = True,
) -> Dict:
    """完整流水线 (Stage 1-6 + 传统方法评估)

    Args:
        molecule_names: 分子名称列表
        max_workers: Stage 2 最大进程数
        run_traditional: 是否运行传统方法评估

    Returns:
        包含各阶段结果的字典
    """
    start_time = datetime.datetime.now()
    results = {}

    # =========================================================================
    # Stage 1: PubChem 查询
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"  Stage 1: PubChem 查询 ({len(molecule_names)} 个分子)")
    logger.info(f"{'='*60}")
    stage1_start = datetime.datetime.now()
    stage1_results = stage1_fetch_all(molecule_names)
    stage1_duration = datetime.datetime.now() - stage1_start

    stage1_success = sum(1 for r in stage1_results if r["status"] == "成功")
    stage1_failed = sum(1 for r in stage1_results if r["status"] != "成功")
    logger.info(f"\n  Stage 1 完成: {stage1_success} 成功, {stage1_failed} 失败 ({stage1_duration})")

    results["stage1"] = {
        "duration": stage1_duration,
        "success": stage1_success,
        "failed": stage1_failed,
    }

    # =========================================================================
    # Stage 2: SwissTarget 靶点预测
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"  Stage 2: SwissTarget 靶点预测")
    logger.info(f"{'='*60}")
    swiss_output_dir = str(PROJECT_ROOT / "src/outputs/shared/swiss_targets")
    Path(swiss_output_dir).mkdir(parents=True, exist_ok=True)

    stage2_start = datetime.datetime.now()
    stage2_results = stage2_predict_all(stage1_results, swiss_output_dir, max_workers)
    stage2_duration = datetime.datetime.now() - stage2_start

    stage2_success = sum(1 for r in stage2_results if r["status"] == "成功")
    stage2_failed = sum(1 for r in stage2_results if r["status"] != "成功")
    logger.info(f"\n  Stage 2 完成: {stage2_success} 成功, {stage2_failed} 失败 ({stage2_duration})")

    results["stage2"] = {
        "duration": stage2_duration,
        "success": stage2_success,
        "failed": stage2_failed,
    }

    # =========================================================================
    # Stage 3: 合并去重
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"  Stage 3: 合并去重靶点")
    logger.info(f"{'='*60}")
    stage3_start = datetime.datetime.now()
    stage3_result = merge_and_deduplicate_targets()
    stage3_duration = datetime.datetime.now() - stage3_start
    logger.info(f"\n  Stage 3 完成 ({stage3_duration})")

    results["stage3"] = {"duration": stage3_duration, "result": stage3_result}

    # =========================================================================
    # Stage 4: 疾病基因交集
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"  Stage 4: 疾病基因交集")
    logger.info(f"{'='*60}")
    stage4_start = datetime.datetime.now()
    stage4_result = merge_disease_genes()
    stage4_duration = datetime.datetime.now() - stage4_start
    logger.info(f"\n  Stage 4 完成 ({stage4_duration})")

    results["stage4"] = {"duration": stage4_duration, "result": stage4_result}

    # =========================================================================
    # Stage 5: 交集分析
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"  Stage 5: 交集分析")
    logger.info(f"{'='*60}")
    stage5_start = datetime.datetime.now()
    stage5_result = intersect_targets_with_disease_genes()
    stage5_duration = datetime.datetime.now() - stage5_start
    logger.info(f"\n  Stage 5 完成 ({stage5_duration})")

    results["stage5"] = {"duration": stage5_duration, "result": stage5_result}

    # =========================================================================
    # Stage 6: PPI 网络构建
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"  Stage 6: PPI 网络构建")
    logger.info(f"{'='*60}")
    stage6_start = datetime.datetime.now()
    stage6_result = build_ppi_network()
    stage6_duration = datetime.datetime.now() - stage6_start
    logger.info(f"\n  Stage 6 完成 ({stage6_duration})")

    results["stage6"] = {"duration": stage6_duration, "result": stage6_result}

    # =========================================================================
    # 分支: 传统方法评估 (Degree Top-30)
    # =========================================================================
    if run_traditional:
        logger.info(f"\n{'='*60}")
        logger.info(f"  传统方法评估 (Degree Top-30 → 富集分析 → 命中率评估)")
        logger.info(f"{'='*60}")

        eval_start = datetime.datetime.now()
        eval_result = evaluate_traditional_method(top_k=30)
        eval_duration = datetime.datetime.now() - eval_start

        results["traditional_eval"] = {
            "evaluation": eval_result,
            "duration": eval_duration,
        }
        logger.info(f"\n  传统方法评估完成 ({eval_duration})")

    # =========================================================================
    # 生成报告
    # =========================================================================
    end_time = datetime.datetime.now()
    total_duration = end_time - start_time

    report_path = PROJECT_ROOT / f"src/outputs/traditional/evaluation/pipeline_report_{end_time.strftime('%Y%m%d_%H%M%S')}.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"# 完整流水线报告\n\n")
        f.write(f"**执行时间**: {start_time.strftime('%Y-%m-%d %H:%M:%S')} - {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**总耗时**: {total_duration}\n\n")
        f.write(f"---\n\n")
        f.write(f"## 各阶段耗时\n\n")
        f.write(f"| Stage | 描述 | 耗时 |\n")
        f.write(f"|-------|------|------|\n")
        f.write(f"| 1 | PubChem 查询 | {results['stage1']['duration']} |\n")
        f.write(f"| 2 | SwissTarget 预测 | {results['stage2']['duration']} |\n")
        f.write(f"| 3 | 合并去重 | {results['stage3']['duration']} |\n")
        f.write(f"| 4 | 疾病基因交集 | {results['stage4']['duration']} |\n")
        f.write(f"| 5 | 交集分析 | {results['stage5']['duration']} |\n")
        f.write(f"| 6 | PPI 网络构建 | {results['stage6']['duration']} |\n")

        if run_traditional:
            f.write(f"| 7 | 传统方法评估 | {results['traditional_eval']['duration']} |\n")

    logger.info(f"\n{'='*60}")
    logger.info(f"  完整流水线完成")
    logger.info(f"{'='*60}")
    logger.info(f"  总耗时: {total_duration}")
    logger.info(f"  报告已保存: {report_path}")

    results["total_duration"] = total_duration
    results["report_path"] = str(report_path)

    return results


if __name__ == "__main__":
    setup_logging()
    strawberry_leaf_molecules = STRAWBERRY_LEAF_MOLECULES
    logger.info(f"\n{'#'*60}")
    logger.info(f"# 完整流水线 - 草莓叶活性成分靶点筛选")
    logger.info(f"# 去重分子数: {len(strawberry_leaf_molecules)}")
    logger.info(f"{'#'*60}\n")

    results = run_full_pipeline(
        strawberry_leaf_molecules,
        max_workers=1,
        run_traditional=True,
    )
