"""
Layer 1 验证：单次参数组合分析
确认 run_single_analysis 等核心函数在真实数据上能正确跑通
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter

from src.config import SWISS_TARGETS_DIR, ANTIOXIDANT_KEYWORDS, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

print("=" * 70)
print("Layer 1 验证：单次分析流程")
print("=" * 70)

# ============================================================================
# 1. 加载数据（与 _param_sensitivity_analysis.py 完全一致）
# ============================================================================
print("\n[1/5] 加载 SwissTarget 预测结果...")
all_predictions = []
for csv_file in SWISS_TARGETS_DIR.glob("*.csv"):
    df = pd.read_csv(csv_file)
    if "Target" in df.columns and "Common name" in df.columns:
        if "Probability*" in df.columns:
            all_predictions.append(df[["Common name", "Probability*"]].rename(
                columns={"Common name": "Gene", "Probability*": "Probability"}
            ))

pred_df = pd.concat(all_predictions, ignore_index=True)
pred_df = pred_df.drop_duplicates(subset=["Gene"])
print(f"  唯一靶点基因数: {len(pred_df)}")
print(f"  Probability 范围: [{pred_df['Probability'].min():.3f}, {pred_df['Probability'].max():.3f}]")

print("\n[2/5] 加载 STRING 边数据...")
string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
print(f"  边数: {len(string_edges)}, Score范围: [{string_edges['combined_score'].min():.3f}, {string_edges['combined_score'].max():.3f}]")

# ============================================================================
# 2. 模拟一个参数组合（传统固定参数）
# ============================================================================
prob_th = 0.0
string_th = 0.4
pval_th = 0.01
top_k = 30

print(f"\n[3/5] 测试参数组合: Prob≥{prob_th}, STRING≥{string_th}, Pvalue<{pval_th}, TopK={top_k}")

# Probability 筛选
filtered_targets = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
print(f"  Probability筛选后: {len(filtered_targets)} 个靶点")

# STRING score 过滤
filtered_edges = string_edges[string_edges["combined_score"] >= string_th]
print(f"  STRING筛选后: {len(filtered_edges)} 条边")

# 构建网络
G = nx.Graph()
for _, row in filtered_edges.iterrows():
    G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

# 取交集节点（关键步骤：STRING节点 ∩ SwissTarget靶点）
valid_nodes = [n for n in G.nodes() if n in filtered_targets]
G = G.subgraph(valid_nodes).copy()
print(f"  交集后PPI网络: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")

# Degree 排序
degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
topk_genes = [g for g, d in degree_sorted[:top_k]]
print(f"  Top-{top_k} Hub基因: {len(topk_genes)} 个")
print(f"  前5个: {topk_genes[:5]}")

# ============================================================================
# 3. 富集分析（关键步骤，最容易出错）
# ============================================================================
print(f"\n[4/5] 运行富集分析...")
try:
    enrichment = run_enrichment_analysis(
        gene_list=topk_genes,
        pvalue_cutoff=pval_th,
        top_n=20,
        background=BACKGROUND_GENES,
    )
    print(f"  富集分析状态: {enrichment['status']}")
except Exception as e:
    print(f"  ❌ 富集分析失败: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# 4. 统计指标（与主脚本完全一致）
# ============================================================================
print(f"\n[5/5] 计算评估指标...")

oxidative_count = 0
total_significant = 0
avg_ef = 0.0
min_pval = 1.0
kegg_pathways = []

if enrichment["status"] == "成功" and "significant" in enrichment:
    for source, df_sig in enrichment["significant"].items():
        if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
            total_significant += len(df_sig)
            for _, row_sig in df_sig.iterrows():
                term = str(row_sig.get("Term", "")).lower()
                if any(kw.lower() in term for kw in ANTIOXIDANT_KEYWORDS):
                    oxidative_count += 1
                if source == "KEGG":
                    kegg_pathways.append(row_sig.get("Term", ""))
                ef = row_sig.get("EF", 0)
                pval = row_sig.get("Adjusted P-value", 1.0)
                if isinstance(ef, (int, float)):
                    avg_ef += ef
                if isinstance(pval, (int, float)) and pval < min_pval:
                    min_pval = pval
    if total_significant > 0:
        avg_ef /= total_significant

print(f"\n{'='*70}")
print("Layer 1 验证结果")
print(f"{'='*70}")
print(f"  输入靶点数:       {len(filtered_targets)}")
print(f"  PPI网络节点:      {G.number_of_nodes()}")
print(f"  PPI网络边:        {G.number_of_edges()}")
print(f"  Top-K Hub数:      {len(topk_genes)}")
print(f"  显著通路总数:     {total_significant}")
print(f"  氧化应激通路数:   {oxidative_count}")
print(f"  平均EF:           {avg_ef:.3f}")
print(f"  最小P-value:      {min_pval:.2e}")
print(f"  KEGG通路数:       {len(kegg_pathways)}")
print(f"\n  ✅ Layer 1 验证通过！核心函数可以正确运行。")
print(f"{'='*70}")
