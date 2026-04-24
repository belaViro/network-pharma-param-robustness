"""
Stage 8 - 富集分析
对PPI网络核心靶点进行GO和KEGG富集分析
"""

import pandas as pd
import gseapy as gp
from pathlib import Path
from typing import Dict, List


def run_enrichment_analysis(
    gene_list: List[str],
    output_dir: str = None,
    pvalue_cutoff: float = 0.01,
    top_n: int = 20,
    background: int = 20500,  # 人类全基因组蛋白编码基因数
    method: str = "shared",  # "shared" or "traditional"
) -> Dict:
    """运行GO和KEGG富集分析

    Args:
        gene_list: 基因符号列表
        output_dir: 输出目录
        pvalue_cutoff: P-value筛选阈值 (默认0.01)
        top_n: 展示前N条通路/条目
        background: 背景基因总数 (默认20500)
        method: 方法标识 ("shared" or "traditional")

    Returns:
        包含富集结果的字典
    """
    if output_dir is None:
        base_dir = Path(__file__).parent.parent / "outputs"
        if method == "shared":
            output_dir = base_dir / "shared" / "enrichment"
        else:
            output_dir = base_dir / method / "enrichment"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"\n{'='*50}")
    logger.info(f"  Stage 8: 富集分析")
    logger.info(f"{'='*50}")
    logger.info(f"  基因数: {len(gene_list)}")
    logger.info(f"  P-value阈值: {pvalue_cutoff}")
    logger.info(f"  背景基因数: {background}")

    # GO富集分析 (Biological Process, Molecular Function, Cellular Component)
    logger.info(f"\n  运行GO富集分析...")

    go_bp = gp.enrichr(
        gene_list=gene_list,
        gene_sets='GO_Biological_Process_2021',
        organism='hsapiens',
        cutoff=pvalue_cutoff,
        outdir=None,
        background=background,
    )

    go_mf = gp.enrichr(
        gene_list=gene_list,
        gene_sets='GO_Molecular_Function_2021',
        organism='hsapiens',
        cutoff=pvalue_cutoff,
        outdir=None,
        background=background,
    )

    go_cc = gp.enrichr(
        gene_list=gene_list,
        gene_sets='GO_Cellular_Component_2021',
        organism='hsapiens',
        cutoff=pvalue_cutoff,
        outdir=None,
        background=background,
    )

    # KEGG富集分析
    logger.info(f"  运行KEGG富集分析...")

    kegg = gp.enrichr(
        gene_list=gene_list,
        gene_sets='KEGG_2021_Human',
        organism='hsapiens',
        cutoff=pvalue_cutoff,
        outdir=None,
        background=background,
    )

    # 保存结果
    results = {
        'GO_BP': go_bp,
        'GO_MF': go_mf,
        'GO_CC': go_cc,
        'KEGG': kegg,
    }

    # 提取显著结果并计算Enrichment Factor
    significant_results = {}
    n_genes = len(gene_list)

    for name, result in results.items():
        if result.res2d is not None and len(result.res2d) > 0:
            df = result.res2d[result.res2d['Adjusted P-value'] < pvalue_cutoff].copy()
            if len(df) > 0:
                # 计算Enrichment Factor: EF = (k/n) / (N/M)
                # Overlap格式: "k/N" 其中k=命中数, N=通路总基因数
                df['k'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
                df['N'] = df['Overlap'].apply(lambda x: int(x.split('/')[1]))
                df['EF'] = (df['k'] / n_genes) / (df['N'] / background)

                df = df.sort_values('Adjusted P-value').head(top_n)
                significant_results[name] = df

                # 保存到CSV
                df.to_csv(output_dir / f"{name}_enrichment.csv", index=False, sep='\t')
                logger.info(f"  {name}: {len(df)} 条显著结果")
            else:
                logger.info(f"  {name}: 无显著结果")
                significant_results[name] = pd.DataFrame()

    logger.info(f"\n  富集结果已保存到: {output_dir}")

    return {
        'status': '成功',
        'gene_count': len(gene_list),
        'results': results,
        'significant': significant_results,
        'output_dir': str(output_dir),
    }


import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    from src.visualization import plot_all_enrichment

    # 读取hub节点
    hub_df = pd.read_csv(Path(__file__).parent.parent / "outputs" / "shared" / "ppi_network" / "hub_nodes.csv")
    gene_list = hub_df['Gene'].tolist()

    logger.info(f"读取到 {len(gene_list)} 个Hub基因")

    # 运行富集分析 (输出到 traditional/enrichment/)
    result = run_enrichment_analysis(gene_list, pvalue_cutoff=0.01, method="traditional")

    # 绘制图表 (输出到 traditional/figures/)
    if result['status'] == '成功':
        traditional_figures_dir = Path(__file__).parent.parent / "outputs" / "traditional" / "figures"
        plot_all_enrichment(result['significant'], output_dir=str(traditional_figures_dir))


__all__ = [
    "run_enrichment_analysis",
]