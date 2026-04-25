"""
富集结果解析与分析帮助函数

提供从富集结果中提取分类指标、汇总统计等工具。
"""

from typing import Dict, List

import pandas as pd

from src.sensitivity.classification import classify_pathway


def analyze_enrichment(enrichment, include_sources=None) -> Dict:
    """从富集结果中提取分类指标。

    Parameters
    ----------
    enrichment : dict
        富集分析结果字典。
    include_sources : list[str], optional
        指定只统计哪些来源的数据库（如 ["KEGG", "GO_BP"]）。
        默认为 None，表示统计全部来源。

    Returns
    -------
    dict with keys:
        total_significant, direct_ox, cross_talk, backbone, other,
        avg_ef, min_pvalue, pathway_details
    """
    category_counts = {
        "direct_ox": 0,
        "cross_talk": 0,
        "signaling_backbone": 0,
        "other": 0,
    }
    total_significant = 0
    avg_ef = 0.0
    min_pval = 1.0
    pathway_details = []

    if enrichment["status"] == "成功" and "significant" in enrichment:
        for source, df_sig in enrichment["significant"].items():
            # 如果指定了 include_sources，跳过不在列表中的来源
            if include_sources is not None and source not in include_sources:
                continue
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                total_significant += len(df_sig)
                for _, row_sig in df_sig.iterrows():
                    term_original = str(row_sig.get("Term", ""))
                    pathway_category = classify_pathway(term_original, source)
                    category_counts[pathway_category] += 1
                    ef = row_sig.get("EF", 0)
                    pval = row_sig.get("Adjusted P-value", 1.0)
                    pathway_details.append({
                        "pathway": term_original,
                        "source": source,
                        "pathway_category": pathway_category,
                        "pvalue": pval,
                        "ef": ef,
                    })
                    if isinstance(ef, (int, float)) and isinstance(pval, (int, float)):
                        avg_ef += ef
                        if pval < min_pval:
                            min_pval = pval
        if total_significant > 0:
            avg_ef /= total_significant

    return {
        "total_significant": total_significant,
        "direct_ox": category_counts["direct_ox"],
        "cross_talk": category_counts["cross_talk"],
        "backbone": category_counts["signaling_backbone"],
        "other": category_counts["other"],
        "avg_ef": round(avg_ef, 2),
        "min_pvalue": min_pval,
        "pathway_details": pathway_details,
    }
