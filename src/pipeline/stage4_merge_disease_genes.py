"""
Stage 4 - 合并疾病相关基因数据库
合并 GeneCards, OMIM, TTD 的 oxidative stress 相关基因
"""

import pandas as pd
from pathlib import Path
from typing import Dict


def load_genecards(file_path: str) -> pd.DataFrame:
    """加载 GeneCards 数据"""
    df = pd.read_csv(file_path)
    df = df.rename(columns={
        "Gene Symbol": "Gene_Symbol",
        "Uniprot ID": "UniProt_ID",
        "Relevance score": "Relevance_Score",
    })
    df = df[["Gene_Symbol", "UniProt_ID", "Relevance_Score", "Description"]].copy()
    df["Source"] = "GeneCards"
    df = df.dropna(subset=["Gene_Symbol"])
    # 清理基因符号
    df["Gene_Symbol"] = df["Gene_Symbol"].str.strip()
    return df


def load_omim(file_path: str) -> pd.DataFrame:
    """加载 OMIM 数据

    OMIM 格式特点:
    - 前5行是元数据
    - * 开头的是基因，# 开头的是疾病表型
    - 基因符号在 Title 字段中，用分号分隔（如 "TITLE; GENE_SYMBOL"）
    """
    # 跳过前5行（元数据），不使用header，手动指定列名
    df = pd.read_csv(file_path, skiprows=5, header=None)
    df.columns = ['MIM Number', 'Title', 'Included Titles', 'Cytogenetic Location', 'Genomic Coordinates', 'Entrez Gene ID']

    # 只保留基因（* 开头），排除疾病（# 开头）
    df = df[df["MIM Number"].astype(str).str.startswith("*")]

    # 从 Title 字段提取基因符号（分号后是基因符号）
    def extract_gene_symbol(title):
        if pd.isna(title):
            return None
        parts = str(title).split(";")
        if len(parts) >= 2:
            gene = parts[-1].strip()
            # 去除末尾的可能疾病标注
            gene = gene.split(",")[0].strip()
            return gene
        return None

    df["Gene_Symbol"] = df["Title"].apply(extract_gene_symbol)
    df = df.dropna(subset=["Gene_Symbol"])

    # 提取 Entrez Gene ID
    df["Entrez_ID"] = pd.to_numeric(df["Entrez Gene ID"], errors="coerce")

    result = df[["Gene_Symbol", "Entrez_ID", "Title"]].copy()
    result.columns = ["Gene_Symbol", "Entrez_ID", "Description"]
    result["Source"] = "OMIM"
    result["Relevance_Score"] = None  # OMIM 没有评分
    result["UniProt_ID"] = None

    return result


def load_ttd(file_path: str) -> pd.DataFrame:
    """加载 TTD 数据

    TTD 格式特点:
    - Gene Names 列可能包含多个基因，用逗号分隔
    """
    df = pd.read_excel(file_path)

    # 从 Gene Names 列提取所有基因
    all_genes = []
    for _, row in df.iterrows():
        gene_names_str = row["Gene Names"]
        if pd.isna(gene_names_str) or gene_names_str == "N.A.":
            continue
        # 分割多个基因
        genes = str(gene_names_str).split(",")
        for gene in genes:
            gene = gene.strip()
            if gene and gene != "N.A.":
                all_genes.append({
                    "Gene_Symbol": gene,
                    "UniProt_ID": None,
                    "Relevance_Score": None,
                    "Description": row["Target Name"],
                    "Source": "TTD",
                })

    return pd.DataFrame(all_genes)


def merge_disease_genes(
    input_dir: str = None,
    output_path: str = None,
) -> Dict:
    """合并所有疾病相关基因数据库

    Args:
        input_dir: 输入目录，默认 src/outputs/shared/disease_genes
        output_path: 输出文件路径

    Returns:
        包含统计信息的字典
    """
    import logging
    logger = logging.getLogger(__name__)

    if input_dir is None:
        from src.config import DISEASE_GENES_DIR
        input_dir = DISEASE_GENES_DIR
    else:
        input_dir = Path(input_dir)

    if output_path is None:
        from src.config import MERGED_DISEASE_GENES_CSV
        output_path = MERGED_DISEASE_GENES_CSV
    else:
        output_path = Path(output_path)

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        logger.info("Stage 4: 合并疾病相关基因数据库")

        all_records = []

        # 1. 加载 GeneCards
        genecards_file = input_dir / "GeneCards-SearchResults.csv"
        if genecards_file.exists():
            try:
                genecards_df = load_genecards(str(genecards_file))
                all_records.append(genecards_df)
                logger.info(f"GeneCards: {len(genecards_df)} 个基因")
            except Exception as e:
                logger.error(f"GeneCards 加载失败: {e}")
        else:
            logger.warning(f"GeneCards: 文件不存在 ({genecards_file})")

        # 2. 加载 OMIM
        omim_file = input_dir / "OMIM-Entry-Search.csv"
        if omim_file.exists():
            try:
                omim_df = load_omim(str(omim_file))
                all_records.append(omim_df)
                logger.info(f"OMIM: {len(omim_df)} 个基因")
            except Exception as e:
                logger.error(f"OMIM 加载失败: {e}")
        else:
            logger.warning(f"OMIM: 文件不存在 ({omim_file})")

        # 3. 加载 TTD
        ttd_file = input_dir / "TTD_Search_Results.xlsx"
        if ttd_file.exists():
            try:
                ttd_df = load_ttd(str(ttd_file))
                all_records.append(ttd_df)
                logger.info(f"TTD: {len(ttd_df)} 个基因")
            except Exception as e:
                logger.error(f"TTD 加载失败: {e}")
        else:
            logger.warning(f"TTD: 文件不存在 ({ttd_file})")

        if not all_records:
            logger.error("[Stage4] 错误: 未找到任何数据文件")
            return {"status": "无文件", "error": "未找到任何数据文件", "total_genes": 0, "unique_genes": 0}

        # 合并
        merged_df = pd.concat(all_records, ignore_index=True)
        logger.info(f"合并后总记录: {len(merged_df)}")

        # 按基因符号去重，优先保留有 UniProt_ID 的记录
        source_priority = {"GeneCards": 0, "OMIM": 1, "TTD": 2}
        merged_df["_priority"] = merged_df["Source"].map(source_priority).fillna(99)
        merged_df = merged_df.sort_values(["_priority", "UniProt_ID"], na_position="last")
        deduplicated_df = merged_df.drop_duplicates(subset=["Gene_Symbol"], keep="first")
        deduplicated_df = deduplicated_df.drop(columns=["_priority"])

        # 按基因符号排序
        deduplicated_df = deduplicated_df.sort_values("Gene_Symbol").reset_index(drop=True)

        # 保存
        deduplicated_df.to_csv(output_path, index=False, encoding="utf-8-sig")

        # 统计各数据库贡献
        source_contribution = deduplicated_df["Source"].value_counts().to_dict()

        logger.info(f"Stage 4 完成: {len(deduplicated_df)} 个去重基因")

        return {
            "status": "成功",
            "total_genes": len(deduplicated_df),
            "output_path": str(output_path),
            "source_contribution": source_contribution,
            "deduplicated_df": deduplicated_df,
        }

    except Exception as e:
        logger.exception(f"[Stage4] 合并疾病基因失败: {e}")
        return {"status": "错误", "error": str(e), "total_genes": 0, "unique_genes": 0}


if __name__ == "__main__":
    result = merge_disease_genes()

    if result["status"] == "成功":
        logger.info(f"\n合并后的基因列表示例 (前10个):")
        logger.info(result["deduplicated_df"][["Gene_Symbol", "Source", "Description"]].head(10).to_string(index=False))