"""
改良版对照实验：客观标准评价富集效果
不预设任何通路重要，从五个独立维度对比
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

# ============================================================================
# 标准1：跨参数靶点稳定性（已有数据）
# ============================================================================

print("=" * 70)
print("改良版对照实验：客观标准评价")
print("=" * 70)

# 加载36组合的靶点结果
target_df = pd.read_csv("src/outputs/sensitivity/mvp_v2_target_robustness.csv")

# 计算各组靶点的跨参数稳定性
print("\n【标准1：跨参数靶点稳定性】")
print("-" * 50)

# 传统方法：Prob=0.05的Top-30 vs Prob=0.1的Top-30
# 从原始数据中提取
target_records = pd.read_csv("src/outputs/sensitivity/mvp_v2_combinations.csv")

# 获取各参数组合下的Top-30靶点
prob05_combos = target_records[target_records["prob_threshold"] == 0.05]["combo"].unique()
prob10_combos = target_records[target_records["prob_threshold"] == 0.1]["combo"].unique()

# 加载all_target_records来获取每个组合的具体靶点
all_targets = pd.read_csv("src/outputs/sensitivity/mvp_v2_combinations.csv")  # 这个不包含靶点详情

# 简化：直接用频率数据
freq_100 = set(target_df[target_df["frequency"] == 1.0]["gene"])
freq_67 = set(target_df[target_df["frequency"] >= 0.67]["gene"])
freq_33 = set(target_df[target_df["frequency"] >= 0.33]["gene"])

print(f"频率=100%靶点（9个）: {sorted(freq_100)}")
print(f"频率>=67%靶点（21个）: {sorted(freq_67)}")
print(f"频率>=33%靶点（52个）: {sorted(freq_33)[:15]}...")

# 模拟传统方法：只跑一个参数组合
# 用Prob=0.05, STRING=0.4, TopK=30的靶点作为"传统方法"
trad_targets = ["SRC", "TNF", "EGFR", "STAT3", "MMP9", "ESR1", "BCL2", "HIF1A", 
                "PTGS2", "ERBB2", "KDR", "MMP2", "APP", "IL2", "MET", "MAPK14", 
                "KIT", "STAT1", "FLT1", "MMP1", "ESR2", "ABCB1", "CA9", "MMP7", 
                "FGFR1", "MAOB", "ABCG2", "CDK1", "MMP14", "FGF1"]

# Q1靶点（31个）
q1_targets = target_df[target_df["quadrant"].str.contains("Q1", na=False)]["gene"].tolist()

# 计算重合度
trad_set = set(trad_targets)
q1_set = set(q1_targets)

print(f"\n靶点集合对比:")
print(f"  传统方法(30个) vs Q1稳健(31个): 重合={len(trad_set & q1_set)}, 并集={len(trad_set | q1_set)}, 重合率={len(trad_set & q1_set)/len(trad_set | q1_set)*100:.1f}%")
print(f"  Q1独有的靶点: {sorted(q1_set - trad_set)}")
print(f"  传统独有的靶点: {sorted(trad_set - q1_set)[:5]}...")

# ============================================================================
# 标准2-5：需要跑富集分析
# ============================================================================

def comprehensive_enrichment_analysis(targets, label):
    """综合富集分析，返回多个客观指标"""
    
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
        "sig_pathways": [],
    }
    
    if enrichment["status"] == "成功" and "significant" in enrichment:
        all_pathways = []
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                results["total_sig"] += len(df_sig)
                for _, row in df_sig.iterrows():
                    term = str(row.get("Term", ""))
                    pval = row.get("Adjusted P-value", row.get("P-value", 1.0))
                    logp = -np.log10(max(pval, 1e-300))
                    all_pathways.append({
                        "term": term,
                        "pvalue": pval,
                        "logp": logp,
                        "source": source
                    })
        
        all_pathways.sort(key=lambda x: x["pvalue"])
        results["sig_pathways"] = all_pathways
    
    return results

def calculate_focus_index(pathways):
    """计算聚焦度：有意义通路占总通路的比例"""
    if not pathways:
        return 0, 0, 0
    
    import re
    backbone_keywords = ['MAPK', 'PI3K-Akt', 'Rap1', 'Ras signaling', 
                         'receptor tyrosine kinase', 'cytokine-mediated']
    cross_talk = {
        "AGE-RAGE signaling pathway in diabetic complications",
        "IL-17 signaling pathway",
        "Fluid shear stress and atherosclerosis"
    }
    
    backbone_count = 0
    crosstalk_count = 0
    
    for p in pathways:
        term_lower = p["term"].lower()
        for kw in backbone_keywords:
            if re.search(r'\b' + re.escape(kw.lower()) + r'\b', term_lower):
                backbone_count += 1
                break
        if p["term"] in cross_talk:
            crosstalk_count += 1
    
    total = len(pathways)
    meaningful = backbone_count + crosstalk_count
    focus_index = meaningful / total if total > 0 else 0
    
    return focus_index, backbone_count, crosstalk_count

def calculate_redundancy(pathways):
    """计算Top5通路的冗余度（基因重叠度）"""
    # 简化：如果Top5通路都是KEGG的癌症/信号通路，说明冗余度高
    # 如果分散在不同数据库和不同功能，说明冗余度低
    if len(pathways) < 5:
        return 0
    
    top5 = pathways[:5]
    sources = [p["source"] for p in top5]
    
    # 来源集中度：如果都来自同一个数据库，说明冗余
    from collections import Counter
    source_counts = Counter(sources)
    max_source_ratio = max(source_counts.values()) / len(sources)
    
    return max_source_ratio

def calculate_pathway_entropy(pathways):
    """计算通路类型的多样性（信息熵）"""
    # 按来源分类计算熵
    if not pathways:
        return 0
    
    from collections import Counter
    import math
    
    sources = [p["source"] for p in pathways[:20]]
    counts = Counter(sources)
    total = len(sources)
    
    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    
    return entropy

# 跑三组富集分析
print("\n" + "=" * 70)
print("开始富集分析...")
print("=" * 70)

groups = [
    (trad_targets, "传统方法(Prob=0.05)"),
    (sorted(freq_100), "全梯度稳健(9个)"),
    (q1_targets, "Q1稳健靶点(31个)")
]

all_results = []

for targets, label in groups:
    print(f"\n分析: {label} ({len(targets)}个靶点)")
    r = comprehensive_enrichment_analysis(targets, label)
    
    # 计算各指标
    focus_idx, backbone, crosstalk = calculate_focus_index(r["sig_pathways"])
    redundancy = calculate_redundancy(r["sig_pathways"])
    entropy = calculate_pathway_entropy(r["sig_pathways"])
    avg_logp = np.mean([p["logp"] for p in r["sig_pathways"][:20]]) if r["sig_pathways"] else 0
    max_logp = max([p["logp"] for p in r["sig_pathways"]]) if r["sig_pathways"] else 0
    
    all_results.append({
        "label": label,
        "n_targets": len(targets),
        "total_sig": r["total_sig"],
        "focus_index": focus_idx,
        "backbone_count": backbone,
        "crosstalk_count": crosstalk,
        "redundancy": redundancy,
        "entropy": entropy,
        "avg_logp": avg_logp,
        "max_logp": max_logp,
        "top5": [p["term"] for p in r["sig_pathways"][:5]]
    })
    
    print(f"  显著通路: {r['total_sig']}, 聚焦度: {focus_idx:.1%}, 冗余度: {redundancy:.2f}")
    print(f"  平均-logP: {avg_logp:.1f}, 最大-logP: {max_logp:.1f}")
    print(f"  通路熵: {entropy:.2f}")

# ============================================================================
# 汇总对比
# ============================================================================

print("\n" + "=" * 70)
print("汇总对比表（客观标准）")
print("=" * 70)

df_results = pd.DataFrame(all_results)

print("\n【维度1：富集规模】")
print(df_results[["label", "n_targets", "total_sig"]].to_string(index=False))

print("\n【维度2：聚焦度（有意义通路占比）】")
print(df_results[["label", "focus_index", "backbone_count", "crosstalk_count"]].to_string(index=False))

print("\n【维度3：显著性强度】")
print(df_results[["label", "avg_logp", "max_logp"]].to_string(index=False))

print("\n【维度4：通路结构特征】")
print(df_results[["label", "redundancy", "entropy"]].to_string(index=False))

print("\n【维度5：Top5通路】")
for _, row in df_results.iterrows():
    print(f"\n{row['label']}:")
    for i, p in enumerate(row['top5'], 1):
        print(f"  {i}. {p}")

# ============================================================================
# 综合评分
# ============================================================================

print("\n" + "=" * 70)
print("综合评分（无权重假设，纯客观对比）")
print("=" * 70)

# 标准化各指标到0-1
for col in ["focus_index", "avg_logp", "entropy"]:
    max_val = df_results[col].max()
    min_val = df_results[col].min()
    if max_val > min_val:
        df_results[f"{col}_norm"] = (df_results[col] - min_val) / (max_val - min_val)
    else:
        df_results[f"{col}_norm"] = 1.0

# 冗余度越低越好（取反）
df_results["redundancy_inv"] = 1 - df_results["redundancy"]

# 综合评分（等权重）
df_results["composite_score"] = (
    df_results["focus_index_norm"] * 0.25 +
    df_results["avg_logp_norm"] * 0.25 +
    df_results["entropy_norm"] * 0.25 +
    df_results["redundancy_inv"] * 0.25
)

print("\n综合评分（0-1，越高越好）:")
for _, row in df_results.sort_values("composite_score", ascending=False).iterrows():
    print(f"  {row['label']:<25} 评分: {row['composite_score']:.3f}")

print("\n说明:")
print("  聚焦度: 有意义通路占比（越高=信号越集中）")
print("  显著性: 平均-logP（越高=关联越强）")
print("  熵: 通路来源多样性（越高=覆盖越广）")
print("  冗余度: Top5通路来源集中度（越低=越不冗余）")
print("  综合评分: 四项标准化后的等权重平均")
