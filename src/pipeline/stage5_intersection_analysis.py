"""
Stage 5 - 靶点与疾病基因交集分析
将 SwissTargetPrediction 靶点与氧化应激疾病基因取交集
"""

import pandas as pd
from pathlib import Path
from typing import Dict
import logging


logger = logging.getLogger(__name__)


def intersect_targets_with_disease_genes(
    targets_path: str = None,
    disease_genes_path: str = None,
    output_path: str = None,
) -> Dict:
    """靶点与疾病基因取交集

    Args:
        targets_path: 靶点文件路径 (merged_targets.csv)
        disease_genes_path: 疾病基因文件路径 (merged_disease_genes.csv)
        output_path: 输出文件路径

    Returns:
        包含统计信息的字典
    """
    if targets_path is None:
        from src.config import MERGED_TARGETS_CSV
        targets_path = MERGED_TARGETS_CSV
    else:
        targets_path = Path(targets_path)

    if disease_genes_path is None:
        from src.config import MERGED_DISEASE_GENES_CSV
        disease_genes_path = MERGED_DISEASE_GENES_CSV
    else:
        disease_genes_path = Path(disease_genes_path)

    if output_path is None:
        from src.config import INTERSECTED_TARGETS_CSV
        output_path = INTERSECTED_TARGETS_CSV
    else:
        output_path = Path(output_path)

    try:
        logger.info("Stage 5: 靶点与疾病基因交集分析")

        # 加载靶点数据
        targets_df = pd.read_csv(targets_path)
        logger.info(f"靶点数 (SwissTarget): {len(targets_df)}")

        # 加载疾病基因数据
        disease_df = pd.read_csv(disease_genes_path)
        logger.info(f"疾病基因数 (GeneCards+OMIM+TTD): {len(disease_df)}")

        # 提取基因符号集合
        target_genes = set(targets_df["Target_Gene"].str.upper())
        disease_genes = set(disease_df["Gene_Symbol"].str.upper())

        # 取交集
        intersection = target_genes & disease_genes
        logger.info(f"交集靶点数: {len(intersection)}")

        # 筛选出交集的靶点详情
        intersected_df = targets_df[
            targets_df["Target_Gene"].str.upper().isin(intersection)
        ].copy()

        # 按概率降序排列
        intersected_df = intersected_df.sort_values(
            "Probability", ascending=False
        ).reset_index(drop=True)

        # 添加疾病来源信息
        def get_disease_source(gene):
            row = disease_df[disease_df["Gene_Symbol"].str.upper() == gene.upper()]
            if len(row) > 0:
                sources = row["Source"].unique().tolist()
                return ", ".join(sources)
            return ""

        intersected_df["Disease_Source"] = intersected_df["Target_Gene"].apply(get_disease_source)

        # 保存结果
        output_path.parent.mkdir(parents=True, exist_ok=True)
        intersected_df.to_csv(output_path, index=False, encoding="utf-8-sig")

        # 统计
        source_stats = intersected_df["Source_Compound"].value_counts().to_dict()
        disease_source_stats = intersected_df["Disease_Source"].value_counts().to_dict()

        logger.info(f"Stage 5 完成: {len(intersected_df)} 个交集靶点")

        return {
            "status": "成功",
            "total_targets": len(targets_df),
            "disease_genes": len(disease_df),
            "intersected_targets": len(intersected_df),
            "output_path": str(output_path),
            "source_stats": source_stats,
            "disease_source_stats": disease_source_stats,
            "intersected_df": intersected_df,
        }

    except FileNotFoundError as e:
        logger.error(f"[Stage5] 文件未找到: {e}")
        return {"status": "错误", "error": f"文件未找到: {e}", "intersected_targets": 0}
    except Exception as e:
        logger.exception(f"[Stage5] 交集分析失败: {e}")
        return {"status": "错误", "error": str(e), "intersected_targets": 0}


if __name__ == "__main__":
    result = intersect_targets_with_disease_genes()

    if result["status"] == "成功":
        logger.info(f"\n交集靶点列表 (前20个):")
        logger.info(result["intersected_df"][["Target_Gene", "Source_Compound", "Probability", "Disease_Source"]].head(20).to_string(index=False))
