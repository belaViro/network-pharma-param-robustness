"""
参数网格敏感性分析

参数网格（72 组合）:
  - Prob (SwissTarget置信度): [0, 0.09, 0.1]  (0.2 已删除，因靶点过少导致富集失效)
  - STRING (PPI置信度): [0.30, 0.35, 0.40, 0.45, 0.50, 0.70]
  - TopK (Hub基因数): [25, 30, 35, 40]
Pvalue 固定为 0.01 用于网格扫描，另单独对比 Pvalue=0.05 对稳健靶点的影响。

用法:
    cd code && python scripts/analysis/run_sensitivity_analysis.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd
import numpy as np
import networkx as nx
from itertools import product
import time
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

from src.config import SWISS_TARGETS_DIR, BACKGROUND_GENES, MERGED_DISEASE_GENES_CSV
from src.sensitivity.enrichment_cache import cached_enrichment
from src.sensitivity.analysis import analyze_enrichment

# ============================================================================
# 参数网格
# ============================================================================
PROB_THRESHOLDS = [0, 0.09, 0.1]
STRING_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.70]
PVALUE_CUTOFF = 0.01   # 网格扫描固定值，后文单独对比 0.05
TOP_K_VALUES = [25, 30, 35, 40]

PARAM_COMBINATIONS = list(product(PROB_THRESHOLDS, STRING_THRESHOLDS, TOP_K_VALUES))

print("=" * 70)
print("参数网格敏感性分析")
print("=" * 70)
print(f"\n参数组合数: {len(PARAM_COMBINATIONS)}")
print(f"  Prob (SwissTarget): {PROB_THRESHOLDS}  (0.2 已排除，因进入富集失效区间)")
print(f"  STRING (PPI):       {STRING_THRESHOLDS}")
print(f"  Pvalue (固定):      {PVALUE_CUTOFF}，另单独对比 0.05")
print(f"  TopK (Hub基因):     {TOP_K_VALUES}")

# ============================================================================
# 加载数据
# ============================================================================
print("\n加载数据...")
pred_df = pd.read_csv(SWISS_TARGETS_DIR / "merged_targets.csv")
pred_df = pred_df[["Target_Gene", "Probability"]].rename(
    columns={"Target_Gene": "Gene"}
)
print(f"  SwissTarget 并集: {len(pred_df)} 个靶点")

# 加载疾病基因 & 预过滤
disease_genes_df = pd.read_csv(MERGED_DISEASE_GENES_CSV)
disease_genes_set = set(disease_genes_df["Gene_Symbol"].str.upper())
pred_df = pred_df[pred_df["Gene"].str.upper().isin(disease_genes_set)].copy()
print(f"  疾病交集后: {len(pred_df)} 个靶点 (固定)")

for p in PROB_THRESHOLDS:
    n = len(pred_df[pred_df["Probability"] >= p])
    print(f"    Prob>={p}: {n} 个靶点")

string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
print(f"  STRING: {len(string_edges)} 条边")

# ============================================================================
# 单次分析
# ============================================================================
def run_single(filtered_targets, string_edges_df, string_th, pval_th, background, top_k=30):
    n_targets = len(filtered_targets)
    if n_targets < 10:
        return None

    filtered_edges = string_edges_df[string_edges_df["combined_score"] >= string_th]
    if len(filtered_edges) == 0:
        return None

    G = nx.Graph()
    for _, row in filtered_edges.iterrows():
        G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

    valid_nodes = [n for n in G.nodes() if n in filtered_targets]
    G = G.subgraph(valid_nodes).copy()
    if G.number_of_nodes() < 10:
        return None

    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    topk = [g for g, d in degree_sorted[:top_k]]
    if len(topk) < 5:
        return None

    degree_dict = dict(G.degree())
    betweenness_dict = nx.betweenness_centrality(G)

    topk_info = []
    for i, (g, d) in enumerate(degree_sorted[:top_k]):
        topk_info.append({
            "gene": g,
            "degree": d,
            "degree_rank": i + 1,
            "betweenness": betweenness_dict.get(g, 0),
        })

    try:
        enrichment = cached_enrichment(topk, pval_th, background)
    except Exception as e:
        logger.warning(f"    富集分析异常: {e}")
        return None

    metrics = analyze_enrichment(enrichment, include_sources=["KEGG", "GO_BP"])

    return {
        "n_targets": n_targets,
        "n_ppi_nodes": G.number_of_nodes(),
        "n_ppi_edges": G.number_of_edges(),
        "top_k": len(topk),
        "topk_info": topk_info,
        "significant_pathways": metrics["total_significant"],
        "direct_ox_count": metrics["direct_ox"],
        "cross_talk_count": metrics["cross_talk"],
        "backbone_count": metrics["backbone"],
        "other_count": metrics["other"],
        "avg_ef": metrics["avg_ef"],
        "min_pvalue": metrics["min_pvalue"],
        "pathway_details": metrics["pathway_details"],
    }

# ============================================================================
# 批量运行
# ============================================================================
print("\n开始批量运行...")
results_records = []
all_target_records = []
all_pathway_records = []

start_time = time.time()

for idx, (prob_th, string_th, top_k) in enumerate(PARAM_COMBINATIONS):
    combo_str = f"P{prob_th}_S{string_th}_K{top_k}"
    filtered_targets = pred_df[pred_df["Probability"] > prob_th]["Gene"].tolist()
    result = run_single(filtered_targets, string_edges, string_th, PVALUE_CUTOFF, BACKGROUND_GENES, top_k)

    if result is None:
        print(f"  [{idx+1:2d}/{len(PARAM_COMBINATIONS)}] SKIP  {combo_str}")
        continue

    for tinfo in result["topk_info"]:
        all_target_records.append({
            **tinfo,
            "combo": combo_str,
            "prob_th": prob_th,
            "string_th": string_th,
            "pval_th": PVALUE_CUTOFF,
            "top_k": top_k,
        })

    for pdet in result["pathway_details"]:
        all_pathway_records.append({
            **pdet,
            "combo": combo_str,
            "prob_th": prob_th,
            "string_th": string_th,
            "pval_th": PVALUE_CUTOFF,
            "top_k": top_k,
            "total_significant": result["significant_pathways"],
        })

    results_records.append({
        "combo": combo_str,
        "prob_threshold": prob_th,
        "string_threshold": string_th,
        "pvalue_cutoff": PVALUE_CUTOFF,
        "top_k": top_k,
        "n_targets": result["n_targets"],
        "n_ppi_nodes": result["n_ppi_nodes"],
        "n_ppi_edges": result["n_ppi_edges"],
        "significant_pathways": result["significant_pathways"],
        "direct_ox_count": result["direct_ox_count"],
        "cross_talk_count": result["cross_talk_count"],
        "backbone_count": result["backbone_count"],
        "other_count": result["other_count"],
        "avg_ef": result["avg_ef"],
        "min_pvalue": result["min_pvalue"],
    })

    print(f"  [{idx+1:2d}/{len(PARAM_COMBINATIONS)}] OK    {combo_str} -> "
          f"靶点{result['n_targets']:3d} PPI{result['n_ppi_nodes']:3d}/{result['n_ppi_edges']:4d} "
          f"主干{result['backbone_count']:2d} 交叉{result['cross_talk_count']:2d} "
          f"直接{result['direct_ox_count']:2d} 总{result['significant_pathways']:3d}")

elapsed = time.time() - start_time
print(f"\n完成: {len(results_records)}/{len(PARAM_COMBINATIONS)} 个成功, 耗时 {elapsed:.1f}s")

# ============================================================================
# 保存结果
# ============================================================================
results_df = pd.DataFrame(results_records)
results_df.to_csv("src/outputs/sensitivity/param_grid_results.csv", index=False)
print(f"\n组合结果已保存: src/outputs/sensitivity/param_grid_results.csv")

# ============================================================================
# 靶点稳健性分析（修复版：Degree 标准化 + 成功组合分母）
# ============================================================================
print("\n" + "=" * 70)
print("靶点稳健性分析")
print("=" * 70)

target_df = pd.DataFrame(all_target_records)
# 修复：分母使用成功运行的组合数，而非总组合数
n_total = len(results_records)

# ---- 新增：计算每个组合的最大 Degree ----
combo_max_degree = {}
for combo, group in target_df.groupby("combo"):
    combo_max_degree[combo] = group["degree"].max()

# ---- 新增：计算标准化 Degree（相对 Degree = degree / max_degree_in_combo） ----
target_df["degree_norm"] = target_df.apply(
    lambda row: row["degree"] / combo_max_degree[row["combo"]], axis=1
)

target_stats = []
for gene, group in target_df.groupby("gene"):
    freq = len(group) / n_total
    avg_degree = group["degree"].mean()
    avg_degree_norm = group["degree_norm"].mean()  # 新增：标准化 Degree
    avg_degree_rank = group["degree_rank"].mean()
    avg_betweenness = group["betweenness"].mean()
    target_stats.append({
        "gene": gene,
        "frequency": freq,
        "count": len(group),
        "avg_degree": avg_degree,
        "avg_degree_norm": avg_degree_norm,  # 新增
        "avg_degree_rank": avg_degree_rank,
        "avg_betweenness": avg_betweenness,
        "min_rank": group["degree_rank"].min(),
        "max_rank": group["degree_rank"].max(),
    })

target_stats_df = pd.DataFrame(target_stats).sort_values("frequency", ascending=False)

freq_med = target_stats_df["frequency"].median()
# 修复：使用标准化 Degree 的中位数作为象限分割线
degree_med = target_stats_df["avg_degree_norm"].median()

def classify_target(row):
    # 修复：使用标准化 Degree 判断拓扑重要性
    if row["frequency"] >= freq_med and row["avg_degree_norm"] >= degree_med:
        return "Q1: 核心稳健"
    elif row["frequency"] < freq_med and row["avg_degree_norm"] >= degree_med:
        return "Q2: 低频高重要"
    elif row["frequency"] >= freq_med and row["avg_degree_norm"] < degree_med:
        return "Q3: 稳健但边缘"
    else:
        return "Q4: 可忽略"

target_stats_df["quadrant"] = target_stats_df.apply(classify_target, axis=1)

print(f"\n总组合数: {n_total} (成功运行)")
print(f"频率中位数: {freq_med:.2f}")
print(f"平均Degree中位数(原始): {target_stats_df['avg_degree'].median():.1f}")
print(f"平均Degree中位数(标准化): {degree_med:.3f}")

q1_targets = target_stats_df[target_stats_df["quadrant"].str.contains("Q1")]["gene"].tolist()
print(f"\nQ1 核心稳健靶点数: {len(q1_targets)}")
if q1_targets:
    print(f"  {', '.join(q1_targets)}")

target_stats_df.to_csv("src/outputs/sensitivity/target_robustness.csv", index=False)

# ============================================================================
# 通路排名稳定性分析
# ============================================================================
print("\n" + "=" * 70)
print("通路排名稳定性分析")
print("=" * 70)

pathway_df = pd.DataFrame(all_pathway_records)

pathway_ranks = []
for combo, group in pathway_df.groupby("combo"):
    sorted_group = group.sort_values("pvalue").reset_index(drop=True)
    for idx, row in sorted_group.iterrows():
        pathway_ranks.append({
            "pathway": row["pathway"],
            "source": row["source"],
            "combo": combo,
            "rank": idx + 1,
            "pvalue": row["pvalue"],
            "ef": row["ef"],
            "logp": -np.log10(max(row["pvalue"], 1e-300)),
            "pathway_category": row["pathway_category"],
        })

rank_df = pd.DataFrame(pathway_ranks)

pathway_rank_stats = []
for pathway, group in rank_df.groupby("pathway"):
    category_counts = group["pathway_category"].value_counts()
    dominant_category = category_counts.index[0]
    pathway_rank_stats.append({
        "pathway": pathway,
        "source": group["source"].iloc[0],
        "pathway_category": dominant_category,
        "frequency": len(group) / n_total,
        "count": len(group),
        "avg_rank": group["rank"].mean(),
        "min_rank": group["rank"].min(),
        "max_rank": group["rank"].max(),
        "rank_std": group["rank"].std(),
        "avg_logp": group["logp"].mean(),
        "avg_ef": group["ef"].mean(),
    })

rank_stats_df = pd.DataFrame(pathway_rank_stats).sort_values("avg_rank")
rank_stats_df.to_csv("src/outputs/sensitivity/pathway_rank_stability.csv", index=False)

CATEGORY_LABELS = {
    "direct_ox": "直接OX", "cross_talk": "交叉调控",
    "signaling_backbone": "信号主干", "other": "其他",
}

print(f"\n{'通路':<50} {'类别':<12} {'频率':<8} {'平均排名':<10} {'排名波动':<10} {'平均EF':<8}")
print("-" * 100)
for _, row in rank_stats_df.head(25).iterrows():
    cat_label = CATEGORY_LABELS.get(row["pathway_category"], row["pathway_category"])
    print(f"{row['pathway']:<48} {cat_label:<12} {row['frequency']:<8.2f} "
          f"{row['avg_rank']:<10.1f} {row['rank_std']:<10.1f} {row['avg_ef']:<8.1f}")

# ============================================================================
# Pvalue 阈值敏感性对比
# ============================================================================
print("\n" + "=" * 70)
print("Pvalue 阈值敏感性对比（Q1 稳健靶点）")
print("=" * 70)

q1_genes = sorted(target_stats_df[target_stats_df["quadrant"].str.contains("Q1")]["gene"].tolist())
print(f"\nQ1 核心稳健靶点: {len(q1_genes)} 个")

pval_results = {}
for pval in [0.01, 0.05]:
    print(f"\n>> Pvalue = {pval}")
    enrichment = cached_enrichment(q1_genes, pval, BACKGROUND_GENES)
    metrics = analyze_enrichment(enrichment, include_sources=["KEGG", "GO_BP"])
    pval_results[pval] = metrics
    print(f"  显著通路总数: {metrics['total_significant']}")
    print(f"  信号主干道:   {metrics['backbone']}")
    print(f"  交叉调控:     {metrics['cross_talk']}")
    print(f"  直接OX:       {metrics['direct_ox']}")
    print(f"  平均EF:       {metrics['avg_ef']:.2f}")
    print(f"  最小P-value:  {metrics['min_pvalue']:.2e}")

# 对比
m01 = pval_results[0.01]
m05 = pval_results[0.05]
print(f"\n--- 对比 ---")
print(f"  P=0.01 → 显著 {m01['total_significant']} 条通路")
print(f"  P=0.05 → 显著 {m05['total_significant']} 条通路")
print(f"  主干道:     P=0.01={m01['backbone']}, P=0.05={m05['backbone']}")
print(f"  交叉调控:   P=0.01={m01['cross_talk']}, P=0.05={m05['cross_talk']}")
print(f"  直接OX:     P=0.01={m01['direct_ox']}, P=0.05={m05['direct_ox']}")
shared = set(p["pathway"] for p in m01["pathway_details"]) & set(p["pathway"] for p in m05["pathway_details"])
if shared:
    print(f"  共有通路: {len(shared)} 条")
pval_stable = m01["backbone"] == m05["backbone"] and m01["cross_talk"] == m05["cross_talk"]
print(f"  分类分布一致: {'是 (Pvalue 稳健)' if pval_stable else '否 (Pvalue 敏感)'}")

# ============================================================================
# 最终判定
# ============================================================================
print("\n" + "=" * 70)
print("敏感性分析结果判定")
print("=" * 70)

checks = []
n_q1 = len(q1_targets)
checks.append(("Q1核心稳健靶点数>=5", n_q1 >= 5, f"{n_q1}"))

key_targets = ["TNF", "SRC", "EGFR", "STAT3", "BCL2", "MMP9", "HIF1A", "PTGS2"]
key_in_q1 = [t for t in key_targets if t in q1_targets]
checks.append(("关键氧化应激靶点进入Q1", len(key_in_q1) >= 3, f"{len(key_in_q1)}个({','.join(key_in_q1)})"))

stable_pathways = rank_stats_df[
    (rank_stats_df["frequency"] >= 0.5) & (rank_stats_df["avg_rank"] <= 10)
]
checks.append(("排名稳定通路数>=3", len(stable_pathways) >= 3, f"{len(stable_pathways)}"))

backbone_stable = stable_pathways[stable_pathways["pathway_category"] == "signaling_backbone"]
checks.append(("信号网络主干道排名稳定", len(backbone_stable) >= 1, f"{len(backbone_stable)}个"))

print(f"\n{'检查项':<40} {'结果':<30} {'判定'}")
print("-" * 90)
pass_count = sum(1 for _, passed, _ in checks if passed)
for check_name, passed, actual in checks:
    print(f"{check_name:<40} {actual:<30} {'PASS' if passed else 'FAIL'}")

print(f"\n通过检查: {pass_count}/{len(checks)}")
if pass_count >= 4:
    print("\n[判定] 敏感性分析通过! 结果支持故事成立。")
elif pass_count >= 2:
    print("\n[判定] 敏感性分析部分通过。")
else:
    print("\n[判定] 敏感性分析未通过。")

# ============================================================================
# 生成 MD 报告
# ============================================================================
print("\n生成详细报告...")
from datetime import datetime

CATEGORY_LABELS_MD = {
    "direct_ox": "直接氧化应激",
    "cross_talk": "交叉调控",
    "signaling_backbone": "信号主干通路",
    "other": "其他",
}

md_lines = []
md_lines.append(f"# 参数网格敏感性分析报告")
md_lines.append(f"")
md_lines.append(f"**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
md_lines.append(f"")
md_lines.append(f"## 1. 参数网格设置")
md_lines.append(f"")
md_lines.append(f"| 参数 | 取值 | 组合数 |")
md_lines.append(f"|------|------|--------|")
md_lines.append(f"| Prob (SwissTarget) | {PROB_THRESHOLDS} | 3 |")
md_lines.append(f"| STRING (PPI) | {STRING_THRESHOLDS} | 6 |")
md_lines.append(f"| TopK (Hub基因) | {TOP_K_VALUES} | 4 |")
md_lines.append(f"| Pvalue 固定 | {PVALUE_CUTOFF} | 1 |")
md_lines.append(f"| **总计** | | **{len(PARAM_COMBINATIONS)}** |")
md_lines.append(f"")

md_lines.append(f"## 2. 数据加载概况")
md_lines.append(f"")
md_lines.append(f"- SwissTarget 并集: {len(pd.read_csv(SWISS_TARGETS_DIR / 'merged_targets.csv'))} 靶点")
md_lines.append(f"- 疾病基因: {len(disease_genes_set)} 个")
md_lines.append(f"- 疾病交集后: {len(pred_df)} 靶点")
md_lines.append(f"- STRING PPI 边数: {len(string_edges)}")
md_lines.append(f"- 成功组合: {len(results_records)} / {len(PARAM_COMBINATIONS)}")
md_lines.append(f"")

md_lines.append(f"## 3. 各概率阈值靶点数")
md_lines.append(f"")
md_lines.append(f"| Prob 阈值 | 疾病相关靶点数 |")
md_lines.append(f"|-----------|---------------|")
for p in PROB_THRESHOLDS:
    n = len(pred_df[pred_df["Probability"] >= p])
    md_lines.append(f"| ≥ {p} | {n} |")
md_lines.append(f"")

# 每组合的富集详情
md_lines.append(f"## 4. 每组参数富集分析详情")
md_lines.append(f"")
md_lines.append(f"| 组合 | 靶点数 | PPI节点 | PPI边 | 显著通路 | 主干道 | 交叉调控 | 直接OX | 其他 | 平均EF | 最小P |")
md_lines.append(f"|------|--------|--------|-------|---------|-------|---------|--------|-------|--------|-------|")
for _, row in results_df.iterrows():
    md_lines.append(
        f"| {row['combo']} | {row['n_targets']} | {row['n_ppi_nodes']} | {row['n_ppi_edges']} | "
        f"{row['significant_pathways']} | {row['backbone_count']} | {row['cross_talk_count']} | "
        f"{row['direct_ox_count']} | {row['other_count']} | {row['avg_ef']:.1f} | {row['min_pvalue']:.2e} |"
    )
md_lines.append(f"")

# 各组合的富集通路明细
md_lines.append(f"### 4.1 各组合富集通路明细")
md_lines.append(f"")
pathway_details_by_combo = {}
for combo_rec in all_pathway_records:
    c = combo_rec["combo"]
    if c not in pathway_details_by_combo:
        pathway_details_by_combo[c] = []
    pathway_details_by_combo[c].append(combo_rec)

for combo, pathways in sorted(pathway_details_by_combo.items()):
    row = results_df[results_df["combo"] == combo].iloc[0]
    md_lines.append(f"<details>")
    md_lines.append(f"<summary><b>{combo}</b> — 显著{row['significant_pathways']}条, 主干{row['backbone_count']}, 交叉{row['cross_talk_count']}, 直接OX{row['direct_ox_count']}, EF={row['avg_ef']:.1f}</summary>")
    md_lines.append(f"")
    md_lines.append(f"| 通路名 | 类别 | P-value | EF |")
    md_lines.append(f"|--------|------|--------|----|")
    sorted_pathways = sorted(pathways, key=lambda x: x["pvalue"])
    for pdet in sorted_pathways[:20]:
        cat = CATEGORY_LABELS_MD.get(pdet["pathway_category"], pdet["pathway_category"])
        md_lines.append(f"| {pdet['pathway'][:55]} | {cat} | {pdet['pvalue']:.2e} | {pdet['ef']:.1f} |")
    if len(sorted_pathways) > 20:
        md_lines.append(f"| ... 还有 {len(sorted_pathways) - 20} 条通路 ... | | | |")
    md_lines.append(f"")
    md_lines.append(f"</details>")
    md_lines.append(f"")

# 靶点稳健性
md_lines.append(f"## 5. 靶点稳健性分析")
md_lines.append(f"")
md_lines.append(f"| 象限 | 说明 | 靶点数 |")
md_lines.append(f"|------|------|--------|")
for q, label in [("Q1: 核心稳健", "高频率 + 高度"), ("Q2: 低频高重要", "低频率 + 高度"),
                  ("Q3: 稳健但边缘", "高频率 + 低度"), ("Q4: 可忽略", "低频率 + 低度")]:
    count = len(target_stats_df[target_stats_df["quadrant"] == q])
    md_lines.append(f"| {q} | {label} | {count} |")
md_lines.append(f"")

md_lines.append(f"### 5.1 Q1 核心稳健靶点列表")
md_lines.append(f"")
q1_sorted = target_stats_df[target_stats_df["quadrant"].str.contains("Q1")].sort_values("frequency", ascending=False)
md_lines.append(f"| 排名 | 基因 | 频率 | 出现次数 | 平均Degree(原始) | 平均Degree(标准化) | 平均排名 | 平均中介中心性 |")
md_lines.append(f"|------|------|------|---------|----------------|-------------------|---------|--------------|")
for rank, (_, r) in enumerate(q1_sorted.iterrows(), 1):
    md_lines.append(f"| {rank} | {r['gene']} | {r['frequency']:.1%} | {int(r['count'])} | {r['avg_degree']:.1f} | {r['avg_degree_norm']:.3f} | {r['avg_degree_rank']:.1f} | {r['avg_betweenness']:.4f} |")
md_lines.append(f"")

# 通路排名稳定性
md_lines.append(f"## 6. 通路排名稳定性 (Top 25)")
md_lines.append(f"")
md_lines.append(f"| 排名 | 通路 | 类别 | 频率 | 平均排名 | 排名波动 | 平均logP | 平均EF |")
md_lines.append(f"|------|------|------|------|---------|---------|---------|-------|")
for _, r in rank_stats_df.head(25).iterrows():
    cat = CATEGORY_LABELS_MD.get(r["pathway_category"], r["pathway_category"])
    md_lines.append(f"| {int(r['avg_rank'])} | {r['pathway'][:50]} | {cat} | {r['frequency']:.1%} | {r['avg_rank']:.1f} | {r['rank_std']:.1f} | {r['avg_logp']:.1f} | {r['avg_ef']:.1f} |")
md_lines.append(f"")

# Pvalue 对比
md_lines.append(f"## 7. Pvalue 阈值敏感性对比")
md_lines.append(f"")
md_lines.append(f"| 指标 | P=0.01 | P=0.05 |")
md_lines.append(f"|------|--------|--------|")
md_lines.append(f"| 显著通路总数 | {m01['total_significant']} | {m05['total_significant']} |")
md_lines.append(f"| 主干道 | {m01['backbone']} | {m05['backbone']} |")
md_lines.append(f"| 交叉调控 | {m01['cross_talk']} | {m05['cross_talk']} |")
md_lines.append(f"| 直接OX | {m01['direct_ox']} | {m05['direct_ox']} |")
md_lines.append(f"| 其他 | {m01['other']} | {m05['other']} |")
md_lines.append(f"| 平均EF | {m01['avg_ef']:.2f} | {m05['avg_ef']:.2f} |")
md_lines.append(f"| 共有通路 | {len(shared)} | {len(shared)} |")
stable = "是（Pvalue 稳健）" if pval_stable else "否（Pvalue 敏感）"
md_lines.append(f"| 分类分布一致 | {stable} | {stable} |")
md_lines.append(f"")

# 最终判定
md_lines.append(f"## 8. 最终判定")
md_lines.append(f"")
md_lines.append(f"| 检查项 | 结果 | 判定 |")
md_lines.append(f"|--------|------|------|")
for check_name, passed, actual in checks:
    md_lines.append(f"| {check_name} | {actual} | {'PASS' if passed else 'FAIL'} |")
md_lines.append(f"")
md_lines.append(f"**通过: {pass_count}/{len(checks)}**")
if pass_count >= 4:
    md_lines.append(f"")
    md_lines.append(f"> ✅ **判定: 敏感性分析通过。结果支持故事成立。**")
elif pass_count >= 2:
    md_lines.append(f"")
    md_lines.append(f"> ⚠️ **判定: 敏感性分析部分通过。**")
else:
    md_lines.append(f"")
    md_lines.append(f"> ❌ **判定: 敏感性分析未通过。**")

# 写入文件
report_md_path = "src/outputs/sensitivity/sensitivity_analysis_report.md"
with open(report_md_path, "w", encoding="utf-8") as f:
    f.write("\n".join(md_lines))
print(f"\nMD 报告已保存: {report_md_path}")
