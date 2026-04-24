"""
传统方法完整流程 (Stage 3-6 + 分支A)

Stage 3: 合并靶点
Stage 4: 合并疾病基因
Stage 5: 交集分析
Stage 6: PPI网络构建
分支A: Degree Top-30 → 富集分析 → 命中率评估 → 图表生成

Usage:
    python scripts/run_traditional_pipeline.py
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.logging_config import setup_logging

logger = logging.getLogger(__name__)

from src.pipeline.stage3_merge_targets import merge_and_deduplicate_targets
from src.pipeline.stage4_merge_disease_genes import merge_disease_genes
from src.pipeline.stage5_intersection_analysis import intersect_targets_with_disease_genes
from src.pipeline.stage6_ppi_network import build_ppi_network
from src.enrichment.enrichment_analysis import run_enrichment_analysis
from src.visualization import plot_all_enrichment
from src.visualization.ppi_network import plot_ppi_network
from scripts.evaluate_traditional import evaluate_traditional_method


PROJECT_ROOT = Path(__file__).parent.parent
SHARED_DIR = PROJECT_ROOT / "src/outputs/shared"
TRADITIONAL_DIR = PROJECT_ROOT / "src/outputs/traditional"


def run_traditional_pipeline() -> dict:
    """运行传统方法完整流程

    Returns:
        包含各阶段结果的字典
    """
    start_time = datetime.now()
    results = {}

    logger.info(f"\n{'#'*60}")
    logger.info(f"# 传统方法流程 - Stage 3-6 + 分支A")
    logger.info(f"{'#'*60}\n")

    # =========================================================================
    # Stage 3: 合并靶点
    # =========================================================================
    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 3: 合并靶点")
    logger.info(f"{'='*50}")
    stage3_start = datetime.now()
    stage3_result = merge_and_deduplicate_targets()
    results["stage3"] = stage3_result
    logger.info(f"  完成: {stage3_result['unique_targets']} 个去重靶点")

    # =========================================================================
    # Stage 4: 合并疾病基因
    # =========================================================================
    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 4: 合并疾病基因")
    logger.info(f"{'='*50}")
    stage4_start = datetime.now()
    stage4_result = merge_disease_genes()
    results["stage4"] = stage4_result
    logger.info(f"  完成: {stage4_result['total_genes']} 个疾病基因")

    # =========================================================================
    # Stage 5: 交集分析
    # =========================================================================
    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 5: 交集分析")
    logger.info(f"{'='*50}")
    stage5_start = datetime.now()
    stage5_result = intersect_targets_with_disease_genes()
    results["stage5"] = stage5_result
    logger.info(f"  完成: {stage5_result['intersected_targets']} 个交集靶点")

    # =========================================================================
    # Stage 6: PPI网络构建
    # =========================================================================
    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 6: PPI网络构建")
    logger.info(f"{'='*50}")
    stage6_start = datetime.now()
    stage6_result = build_ppi_network()
    results["stage6"] = stage6_result
    logger.info(f"  完成: {stage6_result['network_nodes']} 节点, {stage6_result['interactions']} 边")

    # 生成韦恩图 (Stage 5 交集可视化，共享)
    logger.info(f"\n  生成韦恩图...")
    from src.visualization.venn_diagram import plot_venn_diagram
    SHARED_FIGURES_DIR = SHARED_DIR / "figures"
    SHARED_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    venn_path = plot_venn_diagram()
    logger.info(f"  韦恩图已保存: {venn_path}")

    # =========================================================================
    # 分支A: 传统方法 (Degree Top-30 → 富集分析 → 评估)
    # =========================================================================
    logger.info(f"\n{'='*50}")
    logger.info(f"  分支A: 传统方法 (Degree Top-30)")
    logger.info(f"{'='*50}")

    # 读取Top-30 Hub基因
    import pandas as pd
    hub_df = pd.read_csv(SHARED_DIR / "ppi_network/hub_nodes.csv")
    gene_list = hub_df["Gene"].tolist()
    logger.info(f"\n  读取 {len(gene_list)} 个Hub基因")

    # 富集分析 (输出到 traditional/enrichment/)
    logger.info(f"\n  运行富集分析...")
    enrichment_result = run_enrichment_analysis(
        gene_list=gene_list,
        pvalue_cutoff=0.01,
        top_n=20,
        method="traditional",  # 输出到 traditional/enrichment/
    )
    results["enrichment"] = enrichment_result
    logger.info(f"  富集分析完成")

    # 绘制富集图表
    logger.info(f"\n  生成富集分析图表...")
    TRADITIONAL_FIGURES_DIR = TRADITIONAL_DIR / "figures"
    TRADITIONAL_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    if enrichment_result["status"] == "成功":
        plot_all_enrichment(enrichment_result["significant"], output_dir=str(TRADITIONAL_FIGURES_DIR))
        logger.info(f"  富集图表已保存到: {TRADITIONAL_FIGURES_DIR}")

    # 评估
    logger.info(f"\n  运行评估...")
    eval_result = evaluate_traditional_method(top_k=30)
    results["evaluation"] = eval_result
    logger.info(f"\n  评估结果:")
    logger.info(f"    - 氧化应激通路数: {eval_result['oxidative_stress_pathway_count']}")
    logger.info(f"    - 平均 EF: {eval_result['avg_enrichment_factor']}")
    logger.info(f"    - 占比: {eval_result['oxidative_stress_ratio']:.1f}%")

    # =========================================================================
    # 生成报告
    # =========================================================================
    end_time = datetime.now()
    total_duration = end_time - start_time

    report_path = TRADITIONAL_DIR / f"evaluation/pipeline_report_{end_time.strftime('%Y%m%d_%H%M%S')}.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"# 传统方法流程报告\n\n")
        f.write(f"**执行时间**: {start_time.strftime('%Y-%m-%d %H:%M:%S')} - {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**总耗时**: {total_duration}\n\n")
        f.write(f"---\n\n")
        f.write(f"## 数据统计\n\n")
        f.write(f"| Stage | 描述 | 结果 |\n")
        f.write(f"|-------|------|------|\n")
        f.write(f"| 3 | 合并靶点 | {stage3_result['unique_targets']} 个去重靶点 |\n")
        f.write(f"| 4 | 合并疾病基因 | {stage4_result['total_genes']} 个疾病基因 |\n")
        f.write(f"| 5 | 交集分析 | {stage5_result['intersected_targets']} 个交集靶点 |\n")
        f.write(f"| 6 | PPI网络 | {stage6_result['network_nodes']} 节点, {stage6_result['interactions']} 边 |\n")
        f.write(f"| A | 氧化应激通路 | {eval_result['oxidative_stress_pathway_count']} 个 (占比 {eval_result['oxidative_stress_ratio']:.1f}%) |\n")
        f.write(f"| A | 平均EF | {eval_result['avg_enrichment_factor']:.2f} |\n")

    logger.info(f"\n{'='*50}")
    logger.info(f"  传统方法流程完成")
    logger.info(f"{'='*50}")
    logger.info(f"  总耗时: {total_duration}")
    logger.info(f"  报告已保存: {report_path}")

    results["total_duration"] = total_duration
    results["report_path"] = str(report_path)

    return results


if __name__ == "__main__":
    setup_logging()
    results = run_traditional_pipeline()