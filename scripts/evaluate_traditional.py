"""
评估传统方法 (Degree Top-30) 的富集效果

Usage:
    python scripts/evaluate_traditional.py
"""

import pandas as pd
from pathlib import Path
from src.logging_config import setup_logging
import logging

logger = logging.getLogger(__name__)

# 氧化应激关键词
ANTIOXIDANT_KEYWORDS = [
    'oxidative stress',
    'ROS',
    'ferroptosis',
    'glutathione',
    'catalase',
    'peroxidase',
    'superoxide',
    'NF-kappa',
    'MAPK',
    'p53',
    'apoptosis',
    'NADPH oxidase',
    'thioredoxin',
    'lipid peroxidation',
    'autophagy',
    'mitophagy',
    'AGE-RAGE',
    'IL-17',
    'TNF',
    'NOD-like',
]


def is_oxidative_stress_pathway(term: str) -> bool:
    """判断通路是否与氧化应激相关"""
    term_lower = str(term).lower()
    return any(kw.lower() in term_lower for kw in ANTIOXIDANT_KEYWORDS)


def evaluate_traditional_method(
    kegg_path: str = None,
    degrees_path: str = None,
    top_k: int = 30,
) -> dict:
    """评估传统方法 (Degree Top-K)

    Args:
        kegg_path: KEGG富集结果CSV路径
        degrees_path: 节点度数CSV路径
        top_k: 取Top-K个靶点

    Returns:
        评估指标字典
    """
    if kegg_path is None:
        kegg_path = Path(__file__).parent.parent / "src/outputs/traditional/enrichment/KEGG_enrichment.csv"

    if degrees_path is None:
        degrees_path = Path(__file__).parent.parent / "src/outputs/shared/ppi_network/node_degrees.csv"

    # 读取数据
    kegg_df = pd.read_csv(kegg_path, sep='\t')
    degrees_df = pd.read_csv(degrees_path)

    # Top-K by Degree
    top_k_genes = degrees_df.sort_values('Degree', ascending=False).head(top_k)
    top_k_gene_list = top_k_genes['Gene'].tolist()

    logger.info("=" * 60)
    logger.info(f"传统方法评估 (Degree Top-{top_k})")
    logger.info("=" * 60)
    logger.info(f"\nTop-{top_k} Hub Genes:")
    for i, (_, row) in enumerate(top_k_genes.iterrows(), 1):
        logger.info(f"  {i:2d}. {row['Gene']:<15} (Degree: {row['Degree']})")

    # 筛选氧化应激通路
    kegg_df['Is_Oxidative'] = kegg_df['Term'].apply(is_oxidative_stress_pathway)
    oxidative_df = kegg_df[kegg_df['Is_Oxidative']]

    # 统计
    total_significant = len(kegg_df)
    oxidative_count = len(oxidative_df)
    oxidative_ratio = oxidative_count / total_significant * 100 if total_significant > 0 else 0

    # 计算平均EF（氧化应激通路）
    avg_ef = oxidative_df['EF'].mean() if len(oxidative_df) > 0 else 0
    max_ef = oxidative_df['EF'].max() if len(oxidative_df) > 0 else 0
    min_pvalue = kegg_df['Adjusted P-value'].min() if len(kegg_df) > 0 else 1

    results = {
        'method': 'Degree',
        'top_k': top_k,
        'top_k_genes': top_k_gene_list,
        'total_significant_pathways': total_significant,
        'oxidative_stress_pathway_count': oxidative_count,
        'oxidative_stress_ratio': oxidative_ratio,
        'avg_enrichment_factor': round(avg_ef, 2),
        'max_enrichment_factor': round(max_ef, 2),
        'min_pvalue': min_pvalue,
    }

    # 打印结果
    logger.info("\n" + "=" * 60)
    logger.info("评估结果")
    logger.info("=" * 60)
    logger.info(f"显著通路总数: {total_significant}")
    logger.info(f"氧化应激相关通路数: {oxidative_count}")
    logger.info(f"氧化应激通路占比: {oxidative_ratio:.1f}%")
    logger.info(f"平均 Enrichment Factor: {avg_ef:.2f}")
    logger.info(f"最大 Enrichment Factor: {max_ef:.2f}")
    logger.info(f"最小 Adjusted P-value: {min_pvalue:.2e}")

    # 打印氧化应激通路详情
    if len(oxidative_df) > 0:
        logger.info("\n氧化应激相关通路详情:")
        logger.info("-" * 80)
        logger.info(f"{'Rank':<5} {'Term':<45} {'k/N':<10} {'EF':<8} {'P-value'}")
        logger.info("-" * 80)
        for i, (_, row) in enumerate(oxidative_df.iterrows(), 1):
            term = row['Term'][:43] + '..' if len(row['Term']) > 45 else row['Term']
            logger.info(f"{i:<5} {term:<45} {row['Overlap']:<10} {row['EF']:<8.2f} {row['Adjusted P-value']:.2e}")

    return results


if __name__ == "__main__":
    setup_logging()
    results = evaluate_traditional_method(top_k=30)