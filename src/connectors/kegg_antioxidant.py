"""
KEGG抗氧化通路筛选接口
专门用于识别与抗氧化/氧化应激强相关的KEGG通路
"""

import re

import pandas as pd
import gseapy as gp
from pathlib import Path
from typing import List, Dict, Tuple, Optional


# KEGG中与抗氧化直接相关的通路关键词
ANTIOXIDANT_KEYWORDS = [
    'oxidative stress',      # 氧化应激
    'antioxidant',           # 抗氧化
    'ROS',                   # 活性氧
    'heme oxygenase',         # 血红素加氧酶
    'glutathione',            # 谷胱甘肽
    'peroxidase',             # 过氧化物酶
    'superoxide',             # 超氧化物
    'catalase',               # 过氧化氢酶
    'glutathione peroxidase', # GPX
    'glutathione reductase',  # GSR
    'NADPH oxidase',          # NOX
    'thioredoxin',            # 硫氧还蛋白
    'paraoxonase',            # 对氧磷酶
    'lipid peroxidation',     # 脂质过氧化
    'apoptosis',              # 凋亡 (氧化应激可诱导)
    'autophagy',              # 自噬
    'ferroptosis',            # 铁死亡
    'mitophagy',              # 线粒体自噬
    'AGE-RAGE',               # 糖基化终末产物
    'advanced glycosylation', # 糖基化
    'AGE/',                   # AGE-RAGE通路
    'IL-17',                  # IL-17信号通路
    'TNF',                    # 肿瘤坏死因子
    'NF-kappa',               # NF-κB
    'NF-kappaB',              # NF-κB
    'MAPK',                   # MAPK信号
    'PI3K-Akt',               # PI3K-Akt
    'apoptosis',              # 细胞凋亡
    'p53',                    # p53信号通路
    'cell cycle',             # 细胞周期
    'DNA damage',             # DNA损伤
    'DNA repair',             # DNA修复
    'base excision repair',   # 碱基切除修复
    'mismatch repair',        # 错配修复
    'homologous recombination', # 同源重组
    'nucleotide excision repair', # 核苷酸切除修复
    'cancer',                 # 癌症通路(氧化应激参与)
    'diabetes',               # 糖尿病(氧化应激相关)
    'atherosclerosis',        # 动脉粥样硬化
]


def get_antioxidant_pathways(
    gene_list: List[str],
    pvalue_cutoff: float = 0.05,
    top_n: int = 50,
    min_overlap: int = 2,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    获取与抗氧化强相关的KEGG通路

    Args:
        gene_list: 基因符号列表
        pvalue_cutoff: P-value阈值 (默认0.05)
        top_n: 返回前N条最显著的通路
        min_overlap: 最少重叠基因数

    Returns:
        (显著通路DataFrame, 抗氧化相关通路DataFrame)
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"  KEGG抗氧化通路筛选")
    logger.info(f"{'='*60}")
    logger.info(f"  输入基因数: {len(gene_list)}")
    logger.info(f"  P-value阈值: {pvalue_cutoff}")

    # 运行KEGG富集分析
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets='KEGG_2021_Human',
        organism='Human',
        cutoff=pvalue_cutoff,
        outdir=None,
        background=20500,
    )

    if enr.res2d is None or len(enr.res2d) == 0:
        logger.info("  未找到显著富集通路")
        return pd.DataFrame(), pd.DataFrame()

    # 解析结果
    df = enr.res2d.copy()
    df['k'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
    df['N'] = df['Overlap'].apply(lambda x: int(x.split('/')[1]))
    df['EF'] = (df['k'] / len(gene_list)) / (df['N'] / 20500)

    # 筛选显著通路
    sig_df = df[df['Adjusted P-value'] < pvalue_cutoff].copy()
    sig_df = sig_df[sig_df['k'] >= min_overlap]
    sig_df = sig_df.sort_values('Adjusted P-value').head(top_n)

    # 筛选抗氧化相关通路
    def is_antioxidant_related(term: str) -> bool:
        term_lower = term.lower()
        return any(
            re.search(r"\b" + re.escape(kw.lower()) + r"\b", term_lower)
            for kw in ANTIOXIDANT_KEYWORDS
        )

    anti_df = sig_df[sig_df['Term'].apply(is_antioxidant_related)].copy()

    logger.info(f"\n  总显著通路: {len(sig_df)} 条")
    logger.info(f"  抗氧化相关: {len(anti_df)} 条")

    return sig_df, anti_df


def filter_by_keywords(
    pathways_df: pd.DataFrame,
    include_keywords: List[str],
    exclude_keywords: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    根据关键词筛选通路

    Args:
        pathways_df: 通路结果DataFrame
        include_keywords: 必须包含的关键词
        exclude_keywords: 必须不包含的关键词

    Returns:
        筛选后的DataFrame
    """
    if exclude_keywords is None:
        exclude_keywords = []

    def _keyword_match(text: str, keywords: List[str]) -> bool:
        text_lower = text.lower()
        return any(
            re.search(r"\b" + re.escape(kw.lower()) + r"\b", text_lower)
            for kw in keywords
        )

    mask = pathways_df['Term'].apply(
        lambda x: _keyword_match(x, include_keywords)
        and not _keyword_match(x, exclude_keywords)
    )

    return pathways_df[mask].copy()


def save_results(
    sig_df: pd.DataFrame,
    anti_df: pd.DataFrame,
    output_dir: Optional[Path] = None,
) -> Dict[str, str]:
    """保存结果到CSV"""
    if output_dir is None:
        output_dir = Path(__file__).parent.parent.parent / "outputs" / "shared" / "enrichment"

    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}

    if len(sig_df) > 0:
        sig_path = output_dir / "kegg_significant_pathways.csv"
        sig_df.to_csv(sig_path, index=False, sep='\t')
        paths['significant'] = str(sig_path)
        logger.info(f"\n  显著通路已保存: {sig_path}")

    if len(anti_df) > 0:
        anti_path = output_dir / "kegg_antioxidant_pathways.csv"
        anti_df.to_csv(anti_path, index=False, sep='\t')
        paths['antioxidant'] = str(anti_path)
        logger.info(f"  抗氧化通路已保存: {anti_path}")

    return paths


def print_summary(df: pd.DataFrame, title: str = "通路列表") -> None:
    """打印通路摘要"""
    if len(df) == 0:
        logger.info(f"\n  无数据")
        return

    logger.info(f"\n{title} (共{len(df)}条):")
    logger.info("-" * 80)
    logger.info(f"{'Rank':<5} {'Term':<50} {'k/N':<8} {'P-value':<12} {'EF':<8}")
    logger.info("-" * 80)

    for i, (_, row) in enumerate(df.iterrows(), 1):
        term = row['Term'][:48] + '..' if len(row['Term']) > 50 else row['Term']
        logger.info(f"{i:<5} {term:<50} {row['Overlap']:<8} {row['Adjusted P-value']:<12.2e} {row['EF']:<8.2f}")


import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    from pathlib import Path

    # 示例：从hub节点读取基因
    hub_path = Path(__file__).parent.parent.parent / "outputs" / "shared" / "ppi_network" / "hub_nodes.csv"

    if hub_path.exists():
        hub_df = pd.read_csv(hub_path)
        gene_list = hub_df['Gene'].tolist()
        logger.info(f"从hub节点读取 {len(gene_list)} 个基因")
    else:
        # 备选：从merged disease genes
        disease_path = Path(__file__).parent.parent.parent / "outputs" / "shared" / "disease_genes" / "merged_disease_genes.csv"
        if disease_path.exists():
            disease_df = pd.read_csv(disease_path)
            gene_list = disease_df['Gene'].tolist()
            logger.info(f"从疾病基因读取 {len(gene_list)} 个基因")
        else:
            logger.info("未找到基因文件")
            exit(1)

    # 抗氧化通路分析
    sig_df, anti_df = get_antioxidant_pathways(gene_list, pvalue_cutoff=0.05)

    # 打印结果
    print_summary(sig_df, "显著KEGG通路")
    print_summary(anti_df, "抗氧化相关通路")

    # 保存
    save_results(sig_df, anti_df)