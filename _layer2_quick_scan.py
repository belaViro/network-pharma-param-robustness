"""
Layer 2 快速扫描：参数网格点估计（无 Bootstrap）
目的：
  1. 验证全参数网格扫描逻辑正确
  2. 快速得到参数敏感性热图（点估计版）
  3. 识别哪些参数组合会导致结果异常（如靶点过少）
  4. 估算完整运行的耗时

运行时间：预计 10-30 分钟（取决于网络）
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from itertools import product
from collections import Counter
import time
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

from src.config import SWISS_TARGETS_DIR, ANTIOXIDANT_KEYWORDS, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

print("=" * 70)
print("Layer 2: 快速参数扫描（点估计，无 Bootstrap）")
print("=" * 70)

# ============================================================================
# 0. 配置
# ============================================================================
# 精简参数网格（Layer 2 用较少的组合，10-15分钟跑完）
# 完整网格在 Layer 3 使用
PROB_THRESHOLDS = [0.0, 0.05, 0.1, 0.2, 0.3]
STRING_THRESHOLDS = [0.0, 0.3, 0.4, 0.5, 0.7]
PVALUE_CUTOFFS = [0.01, 0.05]
TOP_K_VALUES = [20, 25, 30, 35, 40]

PARAM_COMBINATIONS = list(product(PROB_THRESHOLDS, STRING_THRESHOLDS, PVALUE_CUTOFFS, TOP_K_VALUES))
print(f"\n参数组合数: {len(PARAM_COMBINATIONS)}")
print(f"  Prob: {PROB_THRESHOLDS}")
print(f"  STRING: {STRING_THRESHOLDS}")
print(f"  Pvalue: {PVALUE_CUTOFFS}")
print(f"  TopK: {TOP_K_VALUES}")

# ============================================================================
# 1. 加载数据
# ============================================================================
print("\n加载数据...")
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
print(f"  SwissTarget: {len(pred_df)} 个唯一靶点")

string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
print(f"  STRING: {len(string_edges)} 条边")

# ============================================================================
# 2. 单次分析函数（从主脚本提取，稍作优化）
# ============================================================================
def run_single_analysis(filtered_targets, string_edges_df, string_th, pval_th, background, top_k=30):
    """单次分析流程"""
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

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    if n_nodes < 10:
        return None

    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    topk = [g for g, d in degree_sorted[:top_k]]
    if len(topk) < 5:
        return None

    try:
        enrichment = run_enrichment_analysis(
            gene_list=topk,
            pvalue_cutoff=pval_th,
            top_n=20,
            background=background,
        )
    except Exception as e:
        logger.warning(f"    富集分析异常: {e}")
        return None

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

    return {
        "n_targets": n_targets,
        "n_ppi_nodes": n_nodes,
        "n_ppi_edges": n_edges,
        "top_k": len(topk),
        "top_genes": topk,
        "significant_pathways": total_significant,
        "oxidative_pathways": oxidative_count,
        "avg_ef": avg_ef,
        "min_pvalue": min_pval,
        "kegg_pathways": kegg_pathways,
    }

# ============================================================================
# 3. 批量运行（带进度和计时）
# ============================================================================
print("\n开始批量运行...")
results_records = []
failed_combos = []
start_time = time.time()

for idx, (prob_th, string_th, pval_th, top_k) in enumerate(PARAM_COMBINATIONS):
    combo_str = f"Prob>= {prob_th}, STRING>= {string_th}, Pvalue< {pval_th}, TopK= {top_k}"
    
    # Probability 筛选
    filtered_targets = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
    
    # 单次分析
    result = run_single_analysis(filtered_targets, string_edges, string_th, pval_th, BACKGROUND_GENES, top_k)
    
    if result is None:
        failed_combos.append(combo_str)
        print(f"  [{idx+1:3d}/{len(PARAM_COMBINATIONS)}] SKIP  {combo_str} -> 结果异常/靶点过少")
        continue
    
    results_records.append({
        "prob_threshold": prob_th,
        "string_threshold": string_th,
        "pvalue_cutoff": pval_th,
        "top_k": top_k,
        "n_targets": result["n_targets"],
        "n_ppi_nodes": result["n_ppi_nodes"],
        "n_ppi_edges": result["n_ppi_edges"],
        "significant_pathways": result["significant_pathways"],
        "oxidative_pathways": result["oxidative_pathways"],
        "avg_ef": result["avg_ef"],
        "min_pvalue": result["min_pvalue"],
    })
    
    print(f"  [{idx+1:3d}/{len(PARAM_COMBINATIONS)}] OK    {combo_str} -> "
          f"靶点{result['n_targets']:3d} PPI{result['n_ppi_nodes']:3d}/{result['n_ppi_edges']:4d} "
          f"氧化应激{result['oxidative_pathways']:2d} 总通路{result['significant_pathways']:3d}")

elapsed = time.time() - start_time
print(f"\n批量运行完成: {len(results_records)}/{len(PARAM_COMBINATIONS)} 个成功, 耗时 {elapsed:.1f}s")
if failed_combos:
    print(f"  跳过 {len(failed_combos)} 个异常组合")

# ============================================================================
# 4. 保存结果
# ============================================================================
results_df = pd.DataFrame(results_records)
if len(results_df) == 0:
    print("ERROR: 没有有效结果")
    sys.exit(1)

output_path = "src/outputs/sensitivity/quick_scan_results.csv"
results_df.to_csv(output_path, index=False)
print(f"\n结果已保存: {output_path}")

# ============================================================================
# 5. 参数敏感性分析（单因素分析）
# ============================================================================
print("\n" + "=" * 70)
print("参数敏感性分析")
print("=" * 70)

# 按每个参数分组，计算氧化应激通路数的标准差（敏感性指标）
for param_name, param_col in [
    ("SwissTarget Probability", "prob_threshold"),
    ("STRING Score", "string_threshold"),
    ("P-value cutoff", "pvalue_cutoff"),
    ("Top-K Hub", "top_k"),
]:
    grouped = results_df.groupby(param_col)["oxidative_pathways"].agg(["mean", "std", "min", "max"])
    sensitivity = grouped["std"].mean()  # 平均标准差作为敏感性指标
    print(f"\n{param_name} (敏感性: {sensitivity:.2f}):")
    for val, row in grouped.iterrows():
        print(f"  {param_name}= {val}: 氧化应激={row['mean']:.1f} (std={row['std']:.1f}, range=[{row['min']:.0f}, {row['max']:.0f}])")

# ============================================================================
# 6. 最优参数组合（点估计）
# ============================================================================
print("\n" + "=" * 70)
print("最优参数组合排名（按氧化应激通路数降序）")
print("=" * 70)

best = results_df.sort_values("oxidative_pathways", ascending=False)
print(f"{'Rank':<5} {'Prob':<6} {'STRING':<8} {'Pvalue':<8} {'TopK':<6} "
      f"{'氧化应激':<10} {'总通路':<8} {'PPI节点':<8} {'平均EF':<8}")
print("-" * 70)
for idx, (_, row) in enumerate(best.head(15).iterrows()):
    print(f"{idx+1:<5} {row['prob_threshold']:<6.2f} {row['string_threshold']:<8.1f} "
          f"{row['pvalue_cutoff']:<8.2f} {row['top_k']:<6.0f} {row['oxidative_pathways']:<10.0f} "
          f"{row['significant_pathways']:<8.0f} {row['n_ppi_nodes']:<8.0f} {row['avg_ef']:<8.2f}")

# ============================================================================
# 7. 生成快速热图
# ============================================================================
print("\n生成敏感性热图...")
try:
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 热图1: 氧化应激通路数（固定 Pvalue=0.01, TopK=30）
    ax1 = axes[0, 0]
    sub = results_df[(results_df["pvalue_cutoff"] == 0.01) & (results_df["top_k"] == 30)]
    if len(sub) > 0:
        pivot = sub.pivot_table(values="oxidative_pathways", index="prob_threshold", columns="string_threshold", aggfunc="mean")
        im = ax1.imshow(pivot.values, cmap="YlOrRd", aspect="auto", vmin=0)
        ax1.set_xticks(range(len(pivot.columns)))
        ax1.set_xticklabels([f"{x:.1f}" for x in pivot.columns])
        ax1.set_yticks(range(len(pivot.index)))
        ax1.set_yticklabels([f"{x:.2f}" for x in pivot.index])
        ax1.set_xlabel("STRING Threshold")
        ax1.set_ylabel("Probability Threshold")
        ax1.set_title("Oxidative Stress Pathways\n(Pvalue=0.01, TopK=30)")
        plt.colorbar(im, ax=ax1, fraction=0.046)
    
    # 热图2: 显著通路总数
    ax2 = axes[0, 1]
    if len(sub) > 0:
        pivot2 = sub.pivot_table(values="significant_pathways", index="prob_threshold", columns="string_threshold", aggfunc="mean")
        im2 = ax2.imshow(pivot2.values, cmap="Blues", aspect="auto")
        ax2.set_xticks(range(len(pivot2.columns)))
        ax2.set_xticklabels([f"{x:.1f}" for x in pivot2.columns])
        ax2.set_yticks(range(len(pivot2.index)))
        ax2.set_yticklabels([f"{x:.2f}" for x in pivot2.index])
        ax2.set_xlabel("STRING Threshold")
        ax2.set_ylabel("Probability Threshold")
        ax2.set_title("Total Significant Pathways\n(Pvalue=0.01, TopK=30)")
        plt.colorbar(im2, ax=ax2, fraction=0.046)
    
    # 热图3: Top-K 影响（固定 Pvalue=0.01）
    ax3 = axes[1, 0]
    sub3 = results_df[results_df["pvalue_cutoff"] == 0.01]
    if len(sub3) > 0:
        pivot3 = sub3.pivot_table(values="oxidative_pathways", index="prob_threshold", columns="top_k", aggfunc="mean")
        im3 = ax3.imshow(pivot3.values, cmap="YlOrRd", aspect="auto", vmin=0)
        ax3.set_xticks(range(len(pivot3.columns)))
        ax3.set_xticklabels([f"{x:.0f}" for x in pivot3.columns])
        ax3.set_yticks(range(len(pivot3.index)))
        ax3.set_yticklabels([f"{x:.2f}" for x in pivot3.index])
        ax3.set_xlabel("Top-K")
        ax3.set_ylabel("Probability Threshold")
        ax3.set_title("Oxidative Stress vs Top-K\n(Pvalue=0.01)")
        plt.colorbar(im3, ax=ax3, fraction=0.046)
    
    # 热图4: P-value 影响（固定 TopK=30）
    ax4 = axes[1, 1]
    sub4 = results_df[results_df["top_k"] == 30]
    if len(sub4) > 0:
        pivot4 = sub4.pivot_table(values="oxidative_pathways", index="prob_threshold", columns="pvalue_cutoff", aggfunc="mean")
        im4 = ax4.imshow(pivot4.values, cmap="YlOrRd", aspect="auto", vmin=0)
        ax4.set_xticks(range(len(pivot4.columns)))
        ax4.set_xticklabels([f"{x:.2f}" for x in pivot4.columns])
        ax4.set_yticks(range(len(pivot4.index)))
        ax4.set_yticklabels([f"{x:.2f}" for x in pivot4.index])
        ax4.set_xlabel("P-value cutoff")
        ax4.set_ylabel("Probability Threshold")
        ax4.set_title("Oxidative Stress vs P-value\n(TopK=30)")
        plt.colorbar(im4, ax=ax4, fraction=0.046)
    
    plt.tight_layout()
    fig_path = "src/outputs/sensitivity/quick_scan_heatmaps.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"  热图已保存: {fig_path}")
except Exception as e:
    print(f"  绘图失败（可忽略）: {e}")

print("\n" + "=" * 70)
print("Layer 2 快速扫描完成！")
print("=" * 70)
