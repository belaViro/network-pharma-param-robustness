"""
对照实验：稳健靶点方法 vs 传统单一参数方法
对比三组靶点的富集效果
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from src.config import SWISS_TARGETS_DIR, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# 加载数据
all_predictions = []
for csv_file in SWISS_TARGETS_DIR.glob("*.csv"):
    df = pd.read_csv(csv_file)
    if "Target" in df.columns and "Common name" in df.columns and "Probability*" in df.columns:
        all_predictions.append(df[["Common name", "Probability*"]].rename(
            columns={"Common name": "Gene", "Probability*": "Probability"}
        ))

pred_df = pd.concat(all_predictions, ignore_index=True).drop_duplicates(subset=["Gene"])
string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")

# 加载靶点稳健性结果
robust_df = pd.read_csv("src/outputs/sensitivity/mvp_v2_target_robustness.csv")

# ============================================================================
# 准备三组靶点
# ============================================================================

print("=" * 70)
print("对照实验：稳健靶点 vs 传统靶点")
print("=" * 70)

# 对照组1：传统方法 - 固定参数 Prob=0.05, STRING=0.4, TopK=30
def get_traditional_targets():
    filtered = pred_df[pred_df["Probability"] >= 0.05]["Gene"].tolist()
    edges = string_edges[string_edges["combined_score"] >= 0.4]
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["protein1"], row["protein2"])
    valid = [n for n in G.nodes() if n in filtered]
    G = G.subgraph(valid).copy()
    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    return [g for g, d in degree_sorted[:30]]

# 对照组2：传统方法 - 固定参数 Prob=0.0, STRING=0.4, TopK=30
def get_loose_targets():
    filtered = pred_df[pred_df["Probability"] >= 0.0]["Gene"].tolist()
    edges = string_edges[string_edges["combined_score"] >= 0.4]
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["protein1"], row["protein2"])
    valid = [n for n in G.nodes() if n in filtered]
    G = G.subgraph(valid).copy()
    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    return [g for g, d in degree_sorted[:30]]

# 实验组A：全梯度稳健靶点（频率=100%）
def get_full_robust_targets():
    full = robust_df[robust_df["frequency"] == 1.0].sort_values("avg_degree", ascending=False)
    return full["gene"].tolist()[:30]  # 取Top30

# 实验组B：Q1稳健靶点（高频+高拓扑）
def get_q1_targets():
    q1 = robust_df[robust_df["quadrant"].str.contains("Q1", na=False)].sort_values("frequency", ascending=False)
    return q1["gene"].tolist()[:30]

# 获取各组靶点
trad_targets = get_traditional_targets()
loose_targets = get_loose_targets()
robust_targets = get_full_robust_targets()
q1_targets = get_q1_targets()

print(f"\n对照组1（传统方法，Prob=0.05）: {len(trad_targets)}个靶点")
print(f"  {', '.join(trad_targets[:10])}...")

print(f"\n对照组2（传统方法，Prob=0.0）: {len(loose_targets)}个靶点")
print(f"  {', '.join(loose_targets[:10])}...")

print(f"\n实验组A（全梯度稳健，频率=100%）: {len(robust_targets)}个靶点")
print(f"  {', '.join(robust_targets[:10])}...")

print(f"\n实验组B（Q1稳健靶点）: {len(q1_targets)}个靶点")
print(f"  {', '.join(q1_targets[:10])}...")

# ============================================================================
# 富集分析对比
# ============================================================================

def analyze_enrichment(targets, label):
    """对给定靶点做富集分析，返回核心指标"""
    enrichment = run_enrichment_analysis(
        gene_list=targets,
        pvalue_cutoff=0.01,
        top_n=20,
        background=BACKGROUND_GENES
    )
    
    results = {
        "label": label,
        "n_targets": len(targets),
        "total_sig": 0,
        "backbone_count": 0,
        "cross_talk_count": 0,
        "mapk_rank": None,
        "pi3k_rank": None,
        "mapk_pvalue": None,
        "pi3k_pvalue": None,
        "avg_logp": 0,
        "top5_pathways": []
    }
    
    if enrichment["status"] == "成功" and "significant" in enrichment:
        all_pathways = []
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                results["total_sig"] += len(df_sig)
                for idx, row in df_sig.iterrows():
                    term = str(row.get("Term", ""))
                    pval = row.get("Adjusted P-value", row.get("P-value", 1.0))
                    logp = -np.log10(max(pval, 1e-300))
                    all_pathways.append({
                        "term": term,
                        "pvalue": pval,
                        "logp": logp,
                        "source": source
                    })
        
        # 排序并取Top5
        all_pathways.sort(key=lambda x: x["pvalue"])
        results["top5_pathways"] = [p["term"] for p in all_pathways[:5]]
        
        # 计算平均-logP
        if all_pathways:
            results["avg_logp"] = np.mean([p["logp"] for p in all_pathways[:20]])
        
        # 找MAPK和PI3K-Akt
        for i, p in enumerate(all_pathways):
            if "MAPK signaling pathway" in p["term"]:
                results["mapk_rank"] = i + 1
                results["mapk_pvalue"] = p["pvalue"]
            if "PI3K-Akt signaling pathway" in p["term"]:
                results["pi3k_rank"] = i + 1
                results["pi3k_pvalue"] = p["pvalue"]
        
        # 统计主干通路和交叉调控通路
        backbone_keywords = ['MAPK', 'PI3K-Akt', 'Rap1', 'Ras signaling', 'receptor tyrosine kinase']
        cross_talk = {
            "AGE-RAGE signaling pathway in diabetic complications",
            "IL-17 signaling pathway",
            "Fluid shear stress and atherosclerosis"
        }
        import re
        for p in all_pathways:
            term_lower = p["term"].lower()
            for kw in backbone_keywords:
                if re.search(r'\b' + re.escape(kw.lower()) + r'\b', term_lower):
                    results["backbone_count"] += 1
                    break
            if p["term"] in cross_talk:
                results["cross_talk_count"] += 1
    
    return results

# 跑四组富集分析
print("\n" + "=" * 70)
print("开始富集分析对比...")
print("=" * 70)

results_list = []
for targets, label in [
    (trad_targets, "传统方法(Prob=0.05)"),
    (loose_targets, "传统方法(Prob=0.0)"),
    (robust_targets, "全梯度稳健靶点"),
    (q1_targets, "Q1稳健靶点")
]:
    print(f"\n分析: {label} ({len(targets)}个靶点)")
    r = analyze_enrichment(targets, label)
    results_list.append(r)
    print(f"  显著通路数: {r['total_sig']}, 主干通路: {r['backbone_count']}, 交叉调控: {r['cross_talk_count']}")
    print(f"  MAPK排名: {r['mapk_rank']}, PI3K-Akt排名: {r['pi3k_rank']}")
    print(f"  平均-logP: {r['avg_logp']:.1f}")
    print(f"  Top5: {', '.join(r['top5_pathways'][:3])}...")

# ============================================================================
# 汇总对比表
# ============================================================================

print("\n" + "=" * 70)
print("汇总对比表")
print("=" * 70)

comparison = pd.DataFrame(results_list)
print(comparison[["label", "n_targets", "total_sig", "backbone_count", "cross_talk_count", 
                   "mapk_rank", "pi3k_rank", "avg_logp"]].to_string(index=False))

# ============================================================================
# 跨参数稳定性对比
# ============================================================================

print("\n" + "=" * 70)
print("跨参数稳定性对比")
print("=" * 70)

# 传统方法：Prob=0.05 vs Prob=0.0 的靶点重合度
trad_set = set(trad_targets)
loose_set = set(loose_targets)
robust_set = set(robust_targets)
q1_set = set(q1_targets)

print(f"\n靶点重合度:")
print(f"  传统(Prob=0.05) vs 传统(Prob=0.0): {len(trad_set & loose_set)}/{len(trad_set | loose_set)} = {len(trad_set & loose_set)/len(trad_set | loose_set)*100:.1f}%")
print(f"  传统(Prob=0.05) vs 全梯度稳健: {len(trad_set & robust_set)}/{len(trad_set | robust_set)} = {len(trad_set & robust_set)/len(trad_set | robust_set)*100:.1f}%")
print(f"  传统(Prob=0.0) vs Q1稳健: {len(loose_set & q1_set)}/{len(loose_set | q1_set)} = {len(loose_set & q1_set)/len(loose_set | q1_set)*100:.1f}%")

print(f"\n结论:")
print(f"  传统方法换参数（0.05→0.0），靶点重合度仅 {len(trad_set & loose_set)/len(trad_set | loose_set)*100:.1f}%")
print(f"  说明传统方法结果高度依赖参数选择")
print(f"  稳健靶点是跨参数共识的核心")
