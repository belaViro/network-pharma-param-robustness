"""
Layer 2 快速扫描（带本地缓存版）
核心优化：同一个基因列表 + 同一个Pvalue 只调一次Enrichr API
缓存位置: src/outputs/sensitivity/enrichment_cache/
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from itertools import product
import time
import json
import hashlib
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

from src.config import SWISS_TARGETS_DIR, ANTIOXIDANT_KEYWORDS, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ============================================================================
# 0. 缓存配置
# ============================================================================
CACHE_DIR = Path("src/outputs/sensitivity/enrichment_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

def get_cache_key(gene_list, pvalue_cutoff):
    """为基因列表生成唯一缓存键"""
    genes_str = ";".join(sorted(gene_list))
    hash_str = hashlib.md5(f"{genes_str}|{pvalue_cutoff}".encode()).hexdigest()[:12]
    return hash_str

def load_cached_result(cache_key):
    """从缓存读取富集结果"""
    cache_file = CACHE_DIR / f"{cache_key}.json"
    if cache_file.exists():
        with open(cache_file, "r", encoding="utf-8") as f:
            return json.load(f)
    return None

def save_cached_result(cache_key, result_dict):
    """保存富集结果到缓存"""
    cache_file = CACHE_DIR / f"{cache_key}.json"
    # 将 DataFrame 转为可序列化格式
    serializable = {}
    for source, data in result_dict.items():
        if isinstance(data, pd.DataFrame):
            serializable[source] = data.to_dict(orient="records")
        else:
            serializable[source] = data
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(serializable, f, indent=2)

def cached_enrichment_analysis(gene_list, pvalue_cutoff, top_n=20, background=20500):
    """带缓存的富集分析"""
    cache_key = get_cache_key(gene_list, pvalue_cutoff)
    cached = load_cached_result(cache_key)
    
    if cached is not None:
        # 从缓存重建 DataFrame
        result = {}
        for source, records in cached.items():
            result[source] = pd.DataFrame(records)
        return {"status": "成功", "significant": result}
    
    # 未缓存，调用真实 API
    enrichment = run_enrichment_analysis(
        gene_list=gene_list,
        pvalue_cutoff=pvalue_cutoff,
        top_n=top_n,
        background=background,
    )
    
    if enrichment["status"] == "成功" and "significant" in enrichment:
        save_cached_result(cache_key, enrichment["significant"])
    
    return enrichment

# ============================================================================
# 1. 参数网格（与论文一致）
# ============================================================================
print("=" * 70)
print("Layer 2: 参数扫描（带缓存加速）")
print("=" * 70)

PROB_THRESHOLDS = [0.0, 0.05, 0.1, 0.2, 0.3]
STRING_THRESHOLDS = [0.0, 0.3, 0.4, 0.5, 0.7]
PVALUE_CUTOFFS = [0.01, 0.05]
TOP_K_VALUES = [20, 25, 30, 35, 40]

PARAM_COMBINATIONS = list(product(PROB_THRESHOLDS, STRING_THRESHOLDS, PVALUE_CUTOFFS, TOP_K_VALUES))
print(f"\n参数组合数: {len(PARAM_COMBINATIONS)}")
print(f"  已有缓存: {len(list(CACHE_DIR.glob('*.json')))} 个")

# ============================================================================
# 2. 加载数据
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
# 3. 单次分析函数（使用缓存版富集分析）
# ============================================================================
def run_single_analysis(filtered_targets, string_edges_df, string_th, pval_th, background, top_k=30):
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

    try:
        enrichment = cached_enrichment_analysis(topk, pval_th, top_n=20, background=background)
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
        "n_ppi_nodes": G.number_of_nodes(),
        "n_ppi_edges": G.number_of_edges(),
        "top_k": len(topk),
        "top_genes": topk,
        "significant_pathways": total_significant,
        "oxidative_pathways": oxidative_count,
        "avg_ef": avg_ef,
        "min_pvalue": min_pval,
        "kegg_pathways": kegg_pathways,
    }

# ============================================================================
# 4. 批量运行
# ============================================================================
print("\n开始批量运行（按 Ctrl+C 可随时中断，已缓存结果不会丢失）...")
results_records = []
failed_combos = []
start_time = time.time()

for idx, (prob_th, string_th, pval_th, top_k) in enumerate(PARAM_COMBINATIONS):
    combo_str = f"Prob>={prob_th}, STRING>={string_th}, Pvalue<{pval_th}, TopK={top_k}"
    
    filtered_targets = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
    result = run_single_analysis(filtered_targets, string_edges, string_th, pval_th, BACKGROUND_GENES, top_k)
    
    if result is None:
        failed_combos.append(combo_str)
        print(f"  [{idx+1:3d}/{len(PARAM_COMBINATIONS)}] SKIP  {combo_str}")
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
        "top_genes": ";".join(result["top_genes"]),
        "kegg_pathways": ";".join(result["kegg_pathways"]),
    })
    
    print(f"  [{idx+1:3d}/{len(PARAM_COMBINATIONS)}] OK    {combo_str} -> "
          f"氧化应激{result['oxidative_pathways']:2d} 总通路{result['significant_pathways']:3d}")

elapsed = time.time() - start_time
print(f"\n完成: {len(results_records)}/{len(PARAM_COMBINATIONS)} 个成功, 耗时 {elapsed:.1f}s")

# ============================================================================
# 5. 保存结果
# ============================================================================
results_df = pd.DataFrame(results_records)
if len(results_df) == 0:
    print("ERROR: 没有有效结果")
    sys.exit(1)

output_path = "src/outputs/sensitivity/quick_scan_results.csv"
results_df.to_csv(output_path, index=False)
print(f"结果已保存: {output_path}")

# ============================================================================
# 6. 参数敏感性排名
# ============================================================================
print("\n" + "=" * 70)
print("参数敏感性分析：哪个参数对氧化应激通路数影响最大？")
print("=" * 70)

for param_name, param_col in [
    ("SwissTarget Probability", "prob_threshold"),
    ("STRING Score", "string_threshold"),
    ("P-value cutoff", "pvalue_cutoff"),
    ("Top-K Hub", "top_k"),
]:
    grouped = results_df.groupby(param_col)["oxidative_pathways"].agg(["mean", "std", "min", "max", "count"])
    sensitivity = grouped["std"].mean()
    print(f"\n{param_name} (敏感性指数: {sensitivity:.2f}):")
    for val, row in grouped.iterrows():
        print(f"  {param_name}={val}: 氧化应激={row['mean']:.1f} (std={row['std']:.1f}, n={row['count']})")

# ============================================================================
# 7. 最优参数组合
# ============================================================================
print("\n" + "=" * 70)
print("最优参数组合（按氧化应激通路数降序）")
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
# 8. ★★★ 稳健通路识别 ★★★
# ============================================================================
print("\n" + "=" * 70)
print("稳健通路识别：哪些疾病通路在参数变化下依然稳定出现？")
print("=" * 70)

from collections import Counter

pathway_counter = Counter()
pathway_param_details = {}

for _, row in results_df.iterrows():
    if pd.notna(row["kegg_pathways"]) and row["kegg_pathways"]:
        pathways = row["kegg_pathways"].split(";")
        param_label = f"P{row['prob_threshold']}_S{row['string_threshold']}_V{row['pvalue_cutoff']}_K{row['top_k']}"
        for p in pathways:
            if p:
                pathway_counter[p] += 1
                if p not in pathway_param_details:
                    pathway_param_details[p] = []
                pathway_param_details[p].append(param_label)

n_total = len(results_df)
print(f"\n有效参数组合数: {n_total}")
print(f"\n{'通路名称':<50} {'出现次数':<10} {'频率':<8} {'稳健性评级'}")
print("-" * 90)

robust_pathways = []
for pathway, count in pathway_counter.most_common(30):
    freq = count / n_total
    if freq >= 0.8:
        rating = "★★★ 核心稳健"
    elif freq >= 0.5:
        rating = "★★ 重要稳健"
    elif freq >= 0.3:
        rating = "★ 潜在稳健"
    else:
        rating = "  不稳定"
    
    print(f"{pathway:<48} {count:<10} {freq:<8.2f} {rating}")
    
    if freq >= 0.3:
        robust_pathways.append({
            "pathway": pathway,
            "count": count,
            "frequency": freq,
            "rating": rating.strip(),
        })

# 保存稳健通路
robust_df = pd.DataFrame(robust_pathways)
robust_path = "src/outputs/sensitivity/robust_pathways_layer2.csv"
robust_df.to_csv(robust_path, index=False)
print(f"\n稳健通路已保存: {robust_path}")

# ============================================================================
# 9. 稳健靶点识别
# ============================================================================
print("\n" + "=" * 70)
print("稳健靶点识别：哪些靶点在参数变化下依然稳定出现？")
print("=" * 70)

target_counter = Counter()
for _, row in results_df.iterrows():
    if pd.notna(row["top_genes"]) and row["top_genes"]:
        targets = row["top_genes"].split(";")
        for t in targets:
            if t:
                target_counter[t] += 1

print(f"\n{'靶点':<12} {'出现次数':<10} {'频率':<8} {'稳健性评级'}")
print("-" * 60)

robust_targets = []
for target, count in target_counter.most_common(30):
    freq = count / n_total
    if freq >= 0.8:
        rating = "★★★ 核心稳健"
    elif freq >= 0.5:
        rating = "★★ 重要稳健"
    elif freq >= 0.3:
        rating = "★ 潜在稳健"
    else:
        rating = "  不稳定"
    
    print(f"{target:<10} {count:<10} {freq:<8.2f} {rating}")
    
    if freq >= 0.3:
        robust_targets.append({
            "target": target,
            "count": count,
            "frequency": freq,
            "rating": rating.strip(),
        })

robust_target_df = pd.DataFrame(robust_targets)
robust_target_path = "src/outputs/sensitivity/robust_targets_layer2.csv"
robust_target_df.to_csv(robust_target_path, index=False)
print(f"\n稳健靶点已保存: {robust_target_path}")

# ============================================================================
# 10. 热图
# ============================================================================
print("\n生成热图...")
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
    
    # 热图3: Top-K 影响
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
    
    # 热图4: P-value 影响
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
    print(f"  绘图跳过: {e}")

print("\n" + "=" * 70)
print("Layer 2 完成！")
print("=" * 70)
print(f"\n产出文件:")
print(f"  1. {output_path} - 全部参数组合结果")
print(f"  2. {robust_path} - 稳健通路识别")
print(f"  3. {robust_target_path} - 稳健靶点识别")
print(f"  4. {fig_path} - 敏感性热图")
print(f"  5. {CACHE_DIR} - 富集分析缓存（Layer 3可复用）")
