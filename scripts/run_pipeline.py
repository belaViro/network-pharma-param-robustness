#!/usr/bin/env python3
"""
SwissTargetPrediction 自动化靶点预测流水线

两阶段流水线：
  Stage 1: PubChem API 查询（ThreadPoolExecutor）获取 SMILES
  Stage 2: SwissTargetPrediction 靶点预测（ProcessPoolExecutor）下载 CSV
"""

import sys
import datetime
from pathlib import Path
from typing import List, Dict

# 添加项目根目录到 path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.pipeline import stage1_fetch_all, stage2_predict_all
from src.config import STRAWBERRY_LEAF_MOLECULES
from src.logging_config import setup_logging
import logging

logger = logging.getLogger(__name__)

# 项目根目录
PROJECT_ROOT = Path(__file__).parent.parent


def run_pipeline(
    molecule_names: List[str],
    output_dir: str = None,
    max_workers: int = None,
) -> Dict:
    """完整流水线主函数

    Args:
        molecule_names: 分子名称列表
        output_dir: 输出目录（默认 src/outputs/shared/swiss_targets/）
        max_workers: Stage 2 最大进程数（默认 4）

    Returns:
        包含 stage1_results、stage2_results、report_path 的字典
    """
    if output_dir is None:
        output_dir = str(PROJECT_ROOT / "src" / "outputs" / "shared" / "swiss_targets")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    start_time = datetime.datetime.now()

    # ========== Stage 1 ==========
    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 1: PubChem 查询 ({len(molecule_names)} 个分子)")
    logger.info(f"{'='*50}")
    stage1_start = datetime.datetime.now()
    stage1_results = stage1_fetch_all(molecule_names)
    stage1_duration = datetime.datetime.now() - stage1_start

    stage1_success = sum(1 for r in stage1_results if r["status"] == "成功")
    stage1_failed = sum(1 for r in stage1_results if r["status"] != "成功")
    logger.info(f"\n  Stage 1 完成: {stage1_success} 成功, {stage1_failed} 失败 ({stage1_duration})")

    # ========== Stage 2 ==========
    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 2: Swiss 靶点预测")
    logger.info(f"{'='*50}")
    stage2_start = datetime.datetime.now()
    stage2_results = stage2_predict_all(stage1_results, output_dir, max_workers)
    stage2_duration = datetime.datetime.now() - stage2_start

    stage2_success = sum(1 for r in stage2_results if r["status"] == "成功")
    stage2_failed = sum(1 for r in stage2_results if r["status"] != "成功")
    logger.info(f"\n  Stage 2 完成: {stage2_success} 成功, {stage2_failed} 失败 ({stage2_duration})")

    # ========== 生成报告 ==========
    end_time = datetime.datetime.now()
    total_duration = end_time - start_time

    report_path = str(PROJECT_ROOT / "src" / "outputs" / "reports" / f"pipeline_report_{end_time.strftime('%Y%m%d_%H%M%S')}.md")

    final_results = []
    for r1 in stage1_results:
        r2 = next((r for r in stage2_results if r["smiles"] == r1.get("smiles")), None)
        final_results.append({
            "name": r1["name"],
            "query_name": r1.get("query_name", r1["name"]),
            "smiles": r1.get("smiles", ""),
            "cid": r1.get("cid"),
            "stage1_status": r1["status"],
            "stage2_status": r2["status"] if r2 else "未执行",
            "csv": r2["csv"] if r2 else "",
        })

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"# SwissTargetPrediction 靶点预测流水线报告\n\n")
        f.write(f"**执行时间**: {start_time.strftime('%Y-%m-%d %H:%M:%S')} - {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**总耗时**: {total_duration}\n")
        f.write(f"  - Stage 1 (PubChem查询): {stage1_duration}\n")
        f.write(f"  - Stage 2 (Swiss预测): {stage2_duration}\n\n")
        f.write(f"**总分子数**: {len(molecule_names)}\n")
        f.write(f"**Stage 1 成功**: {stage1_success} | **Stage 1 失败**: {stage1_failed}\n")
        f.write(f"**Stage 2 成功**: {stage2_success} | **Stage 2 失败**: {stage2_failed}\n\n")
        f.write("---\n\n")
        f.write("## 详细结果\n\n")
        f.write("| 序号 | 分子名称 | PubChem查询名 | CID | Stage 1 | Stage 2 | CSV |\n")
        f.write("|------|---------|--------------|-----|---------|--------|-----|\n")
        for i, r in enumerate(final_results, 1):
            cid = r.get("cid")
            if cid:
                cid_link = f"[{cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})"
            else:
                cid_link = "-"
            s1_icon = "✓" if r["stage1_status"] == "成功" else "✗"
            s2_icon = "✓" if r["stage2_status"] == "成功" else "✗"
            csv_name = Path(r["csv"]).name if r["csv"] else "-"
            f.write(f"| {i} | {r['name']} | {r['query_name']} | {cid_link} | {s1_icon} | {s2_icon} | {csv_name} |\n")

    logger.info(f"\n{'='*50}")
    logger.info(f"  流水线完成")
    logger.info(f"{'='*50}")
    logger.info(f"  总耗时: {total_duration}")
    logger.info(f"  报告已保存: {report_path}")

    return {
        "stage1_results": stage1_results,
        "stage2_results": stage2_results,
        "report_path": report_path,
    }


def run_from_names(names: List[str], output_dir: str = None, max_workers: int = 4):
    """从分子名称列表运行流水线"""
    return run_pipeline(names, output_dir, max_workers)


def run_from_smiles(smiles_list: List[str], output_dir: str = None, max_workers: int = 4):
    """从 SMILES 列表直接运行 Stage 2（跳过 Stage 1）"""
    if output_dir is None:
        output_dir = str(PROJECT_ROOT / "src" / "outputs" / "shared" / "swiss_targets")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    stage1_results = [{"name": f"mol_{i}", "smiles": s, "status": "成功"} for i, s in enumerate(smiles_list)]
    return stage2_predict_all(stage1_results, output_dir, max_workers)


if __name__ == "__main__":
    setup_logging()
    strawberry_leaf_molecules = STRAWBERRY_LEAF_MOLECULES
    logger.info(f"\n{'#'*50}")
    logger.info(f"# SwissTargetPrediction 靶点预测流水线")
    logger.info(f"# 草莓叶活性成分: {len(strawberry_leaf_molecules)} 个去重分子")
    logger.info(f"{'#'*50}\n")
    run_pipeline(strawberry_leaf_molecules, max_workers=1)
