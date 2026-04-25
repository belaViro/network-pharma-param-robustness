"""
MVP v2: 修正版
Prob阈值调整为 [0.0, 0.05, 0.1]，去掉0.2
解决Prob=0.2时靶点过少(13个)导致CA基因劫持的问题
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
import re

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

from src.config import SWISS_TARGETS_DIR, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ============================================================================
# 0. 缓存（复用v1的缓存）
# ============================================================================
CACHE_DIR = Path("src/outputs/sensitivity/enrichment_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

def get_cache_key(gene_list, pvalue_cutoff):
    genes_str = ";".join(sorted(gene_list))
    return hashlib.md5(f"{genes_str}|{pvalue_cutoff}".encode()).hexdigest()[:12]

def load_cached(cache_key):
    f = CACHE_DIR / f"{cache_key}.json"
    if f.exists():
        with open(f, "r") as fh:
            data = json.load(fh)
            return {k: pd.DataFrame(v) if isinstance(v, list) else v for k, v in data.items()}
    return None

def save_cached(cache_key, result):
    f = CACHE_DIR / f"{cache_key}.json"
    serializable = {k: (v.to_dict(orient="records") if isinstance(v, pd.DataFrame) else v) for k, v in result.items()}
    with open(f, "w") as fh:
        json.dump(serializable, fh)

def cached_enrichment(gene_list, pvalue_cutoff, background=20500):
    ck = get_cache_key(gene_list, pvalue_cutoff)
    cached = load_cached(ck)
    if cached is not None:
        return {"status": "成功", "significant": cached}
    enrichment = run_enrichment_analysis(gene_list=gene_list, pvalue_cutoff=pvalue_cutoff, top_n=20, background=background)
    if enrichment["status"] == "成功" and "significant" in enrichment:
        save_cached(ck, enrichment["significant"])
    return enrichment

# ============================================================================
# 1. 加载 GO 氧化应激集合 + 定义三级分类体系
# ============================================================================
GO_OX_FILE = Path("src/data/oxidative_stress_go_terms.txt")
GO_OX_IDS = set()
if GO_OX_FILE.exists():
    with open(GO_OX_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            go_id = line.strip().split('\t')[0]
            GO_OX_IDS.add(go_id)
    print(f"[INFO] 加载 GO 氧化应激集合: {len(GO_OX_IDS)} 个 terms")
else:
    print(f"[WARNING] GO 氧化应激集合不存在: {GO_OX_FILE}")

# Curated KEGG/Reactome 交叉调控通路（已知与氧化应激存在功能交叉）
CROSS_TALK_PATHWAYS = {
    "AGE-RAGE signaling pathway in diabetic complications",
    "IL-17 signaling pathway",
    "Ferroptosis",
    "Glutathione metabolism",
    "Peroxisome",
    "Fluid shear stress and atherosclerosis",
    "NOD-like receptor signaling pathway",
    "TNF signaling pathway",
    "Toll-like receptor signaling pathway",
    "NF-kappa B signaling pathway",
}

# 信号网络主干道关键词（核心信号转导通路）
BACKBONE_KEYWORDS = [
    'MAPK', 'PI3K-Akt', 'Rap1', 'Ras signaling',
    'receptor tyrosine kinase', 'cytokine-mediated',
    'calcium signaling', 'cAMP signaling',
    'chemokine signaling', 'JAK-STAT',
    'EGFR tyrosine kinase inhibitor resistance',
]

def classify_pathway(term: str) -> str:
    """
    三级分类体系：
    - direct_ox: GO层级验证，属于 GO:0006979 后代（最严格）
    - cross_talk: 已知与氧化应激存在交叉调控（curated列表）
    - signaling_backbone: 信号网络主干道（核心信号转导）
    - other: 其他通路
    """
    term_lower = term.lower()
    
    # 1. 直接氧化应激响应（GO层级验证）
    match = re.search(r'\(GO:(\d+)\)', term)
    if match:
        go_id = f"GO:{match.group(1)}"
        if go_id in GO_OX_IDS:
            return "direct_ox"
    
    # 2. 交叉调控通路（curated列表）
    if term in CROSS_TALK_PATHWAYS:
        return "cross_talk"
    
    # 3. 信号网络主干道（单词边界匹配）
    for kw in BACKBONE_KEYWORDS:
        if re.search(r'\b' + re.escape(kw.lower()) + r'\b', term_lower):
            return "signaling_backbone"
    
    return "other"


# ============================================================================
# 2. MVP v2参数网格（36个组合）
# ============================================================================
print("=" * 70)
print("MVP v2: 三级分类修正版（Prob=[0.0, 0.05, 0.1]）")
print("=" * 70)

PROB_THRESHOLDS = [0.0, 0.05, 0.1]
STRING_THRESHOLDS = [0.3, 0.4, 0.5, 0.7]
PVALUE_CUTOFFS = [0.01]
TOP_K_VALUES = [25, 30, 40]

PARAM_COMBINATIONS = list(product(PROB_THRESHOLDS, STRING_THRESHOLDS, PVALUE_CUTOFFS, TOP_K_VALUES))
print(f"\n参数组合数: {len(PARAM_COMBINATIONS)}")
print(f"  Prob: {PROB_THRESHOLDS}")
print(f"  STRING: {STRING_THRESHOLDS}")
print(f"  Pvalue: {PVALUE_CUTOFFS}")
print(f"  TopK: {TOP_K_VALUES}")
print(f"  已有缓存: {len(list(CACHE_DIR.glob('*.json')))} 个")

# ============================================================================
# 2. 加载数据
# ============================================================================
print("\n加载数据...")
all_predictions = []
for csv_file in SWISS_TARGETS_DIR.glob("*.csv"):
    df = pd.read_csv(csv_file)
    if "Target" in df.columns and "Common name" in df.columns and "Probability*" in df.columns:
        all_predictions.append(df[["Common name", "Probability*"]].rename(
            columns={"Common name": "Gene", "Probability*": "Probability"}
        ))

pred_df = pd.concat(all_predictions, ignore_index=True).drop_duplicates(subset=["Gene"])
print(f"  SwissTarget: {len(pred_df)} 个唯一靶点")

# 预览不同Prob阈值下的靶点数量
for p in PROB_THRESHOLDS:
    n = len(pred_df[pred_df["Probability"] >= p])
    print(f"    Prob>={p}: {n} 个靶点")

string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
print(f"  STRING: {len(string_edges)} 条边")

# ============================================================================
# 3. 单次分析
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

    # 三级分类计数器
    category_counts = {"direct_ox": 0, "cross_talk": 0, "signaling_backbone": 0, "other": 0}
    total_significant = 0
    avg_ef = 0.0
    min_pval = 1.0
    pathway_details = []

    if enrichment["status"] == "成功" and "significant" in enrichment:
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                total_significant += len(df_sig)
                for idx, row_sig in df_sig.iterrows():
                    term_original = str(row_sig.get("Term", ""))
                    pathway_category = classify_pathway(term_original)
                    category_counts[pathway_category] += 1
                    
                    ef = row_sig.get("EF", 0)
                    pval = row_sig.get("Adjusted P-value", 1.0)
                    logp = -np.log10(max(pval, 1e-300))
                    
                    pathway_details.append({
                        "pathway": term_original,
                        "source": source,
                        "pathway_category": pathway_category,
                        "pvalue": pval,
                        "ef": ef,
                        "logp": logp,
                    })
                    
                    if isinstance(ef, (int, float)) and isinstance(pval, (int, float)):
                        avg_ef += ef
                        if pval < min_pval:
                            min_pval = pval
        if total_significant > 0:
            avg_ef /= total_significant

    return {
        "n_targets": n_targets,
        "n_ppi_nodes": G.number_of_nodes(),
        "n_ppi_edges": G.number_of_edges(),
        "top_k": len(topk),
        "topk_info": topk_info,
        "significant_pathways": total_significant,
        "direct_ox_count": category_counts["direct_ox"],
        "cross_talk_count": category_counts["cross_talk"],
        "backbone_count": category_counts["signaling_backbone"],
        "other_count": category_counts["other"],
        "avg_ef": avg_ef,
        "min_pvalue": min_pval,
        "pathway_details": pathway_details,
    }

# ============================================================================
# 4. 批量运行
# ============================================================================
print("\n开始批量运行...")
results_records = []
all_target_records = []
all_pathway_records = []

start_time = time.time()

for idx, (prob_th, string_th, pval_th, top_k) in enumerate(PARAM_COMBINATIONS):
    combo_str = f"P{prob_th}_S{string_th}_V{pval_th}_K{top_k}"
    filtered_targets = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
    result = run_single(filtered_targets, string_edges, string_th, pval_th, BACKGROUND_GENES, top_k)
    
    if result is None:
        print(f"  [{idx+1:2d}/{len(PARAM_COMBINATIONS)}] SKIP  {combo_str}")
        continue
    
    for tinfo in result["topk_info"]:
        all_target_records.append({
            **tinfo,
            "combo": combo_str,
            "prob_th": prob_th,
            "string_th": string_th,
            "pval_th": pval_th,
            "top_k": top_k,
        })
    
    for pdet in result["pathway_details"]:
        all_pathway_records.append({
            **pdet,
            "combo": combo_str,
            "prob_th": prob_th,
            "pval_th": pval_th,
            "top_k": top_k,
            "total_significant": result["significant_pathways"],
        })
    
    results_records.append({
        "combo": combo_str,
        "prob_threshold": prob_th,
        "string_threshold": string_th,
        "pvalue_cutoff": pval_th,
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
          f"主干{result['backbone_count']:2d} 交叉{result['cross_talk_count']:2d} 直接{result['direct_ox_count']:2d} 总{result['significant_pathways']:3d}")

elapsed = time.time() - start_time
print(f"\n完成: {len(results_records)}/{len(PARAM_COMBINATIONS)} 个成功, 耗时 {elapsed:.1f}s")

# ============================================================================
# 5. 保存原始结果
# ============================================================================
results_df = pd.DataFrame(results_records)
results_df.to_csv("src/outputs/sensitivity/mvp_v2_combinations.csv", index=False)
print(f"\n组合结果已保存: src/outputs/sensitivity/mvp_v2_combinations.csv")

# ============================================================================
# 6. 靶点稳健性分析
# ============================================================================
print("\n" + "=" * 70)
print("靶点稳健性分析 (MVP v2)")
print("=" * 70)

target_df = pd.DataFrame(all_target_records)
n_total = len(PARAM_COMBINATIONS)

target_stats = []
for gene, group in target_df.groupby("gene"):
    freq = len(group) / n_total
    avg_degree = group["degree"].mean()
    avg_degree_rank = group["degree_rank"].mean()
    avg_betweenness = group["betweenness"].mean()
    min_rank = group["degree_rank"].min()
    max_rank = group["degree_rank"].max()
    
    target_stats.append({
        "gene": gene,
        "frequency": freq,
        "count": len(group),
        "avg_degree": avg_degree,
        "avg_degree_rank": avg_degree_rank,
        "avg_betweenness": avg_betweenness,
        "min_rank": min_rank,
        "max_rank": max_rank,
    })

target_stats_df = pd.DataFrame(target_stats).sort_values("frequency", ascending=False)

freq_med = target_stats_df["frequency"].median()
degree_med = target_stats_df["avg_degree"].median()

def classify_target(row):
    if row["frequency"] >= freq_med and row["avg_degree"] >= degree_med:
        return "Q1: 核心稳健"
    elif row["frequency"] < freq_med and row["avg_degree"] >= degree_med:
        return "Q2: 低频高重要"
    elif row["frequency"] >= freq_med and row["avg_degree"] < degree_med:
        return "Q3: 稳健但边缘"
    else:
        return "Q4: 可忽略"

target_stats_df["quadrant"] = target_stats_df.apply(classify_target, axis=1)

print(f"\n总组合数: {n_total}")
print(f"频率中位数: {freq_med:.2f}")
print(f"平均Degree中位数: {degree_med:.1f}")

print(f"\n{'靶点':<12} {'频率':<8} {'平均Degree':<12} {'平均排名':<10} {'象限'}")
print("-" * 60)

q1_targets = []
for _, row in target_stats_df.head(30).iterrows():
    print(f"{row['gene']:<10} {row['frequency']:<8.2f} {row['avg_degree']:<12.1f} {row['avg_degree_rank']:<10.1f} {row['quadrant']}")
    if "Q1" in row["quadrant"]:
        q1_targets.append(row["gene"])

print(f"\nQ1 核心稳健靶点数: {len(q1_targets)}")
if q1_targets:
    print(f"  {', '.join(q1_targets)}")

target_stats_df.to_csv("src/outputs/sensitivity/mvp_v2_target_robustness.csv", index=False)

# ============================================================================
# 7. 通路排名稳定性分析
# ============================================================================
print("\n" + "=" * 70)
print("通路排名稳定性分析 (MVP v2)")
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
            "logp": row["logp"],
            "pathway_category": row["pathway_category"],
        })

rank_df = pd.DataFrame(pathway_ranks)

pathway_rank_stats = []
for pathway, group in rank_df.groupby("pathway"):
    # 通路分类：取该通路在所有组合中出现过的分类（多数投票）
    category_counts = group["pathway_category"].value_counts()
    dominant_category = category_counts.index[0]
    
    avg_rank = group["rank"].mean()
    min_rank = group["rank"].min()
    max_rank = group["rank"].max()
    rank_std = group["rank"].std()
    count = len(group)
    freq = count / n_total
    avg_logp = group["logp"].mean()
    avg_ef = group["ef"].mean()
    
    pathway_rank_stats.append({
        "pathway": pathway,
        "source": group["source"].iloc[0],
        "pathway_category": dominant_category,
        "frequency": freq,
        "count": count,
        "avg_rank": avg_rank,
        "min_rank": min_rank,
        "max_rank": max_rank,
        "rank_std": rank_std,
        "avg_logp": avg_logp,
        "avg_ef": avg_ef,
    })

rank_stats_df = pd.DataFrame(pathway_rank_stats).sort_values("avg_rank")

print(f"\n{'通路':<50} {'来源':<10} {'类别':<12} {'频率':<8} {'平均排名':<10} {'排名波动':<10} {'平均EF':<8}")
print("-" * 120)

CATEGORY_LABELS = {
    "direct_ox": "直接OX",
    "cross_talk": "交叉调控",
    "signaling_backbone": "信号主干",
    "other": "其他",
}

for _, row in rank_stats_df.head(25).iterrows():
    cat_label = CATEGORY_LABELS.get(row["pathway_category"], row["pathway_category"])
    print(f"{row['pathway']:<48} {row['source']:<10} {cat_label:<12} {row['frequency']:<8.2f} "
          f"{row['avg_rank']:<10.1f} {row['rank_std']:<10.1f} {row['avg_ef']:<8.1f}")

rank_stats_df.to_csv("src/outputs/sensitivity/mvp_v2_pathway_rank_stability.csv", index=False)

# ============================================================================
# 8. v1 vs v2 关键对比
# ============================================================================
print("\n" + "=" * 70)
print("v1 vs v2 关键对比")
print("=" * 70)

print("\n【v1 问题】Prob=0.2 时靶点过少:")
print("  Prob=0.0: 282靶点 → 氧化应激通路12条")
print("  Prob=0.2:  13靶点 → 氧化应激通路 0条 (CA基因劫持)")

print("\n【v2 修正】Prob=[0.0, 0.05, 0.1]:")
for p in PROB_THRESHOLDS:
    n = len(pred_df[pred_df["Probability"] >= p])
    print(f"  Prob>={p}: {n}靶点")

# ============================================================================
# 9. 最终判定
# ============================================================================
print("\n" + "=" * 70)
print("MVP v2 测试结果判定")
print("=" * 70)

# 检查清单
checks = []

n_q1 = len(q1_targets)
checks.append(("Q1核心稳健靶点数>=5", n_q1 >= 5, f"{n_q1}"))

key_targets = ["TNF", "SRC", "EGFR", "STAT3", "BCL2", "MMP9", "HIF1A", "PTGS2"]
key_in_q1 = [t for t in key_targets if t in q1_targets]
checks.append(("关键氧化应激靶点进入Q1", len(key_in_q1) >= 3, f"{len(key_in_q1)}个({','.join(key_in_q1)})"))

ox_targets = ["TNF", "MAPK1", "MAPK3", "EGFR", "STAT3", "BCL2", "MMP9", "PTGS2", "HIF1A", "NOS3"]
ox_in_q1 = [t for t in ox_targets if t in q1_targets]
checks.append(("氧化应激靶点进入Q1", len(ox_in_q1) >= 2, f"{len(ox_in_q1)}个"))

stable_pathways = rank_stats_df[
    (rank_stats_df["frequency"] >= 0.5) & 
    (rank_stats_df["avg_rank"] <= 10)
]
n_stable = len(stable_pathways)
checks.append(("排名稳定通路数>=3", n_stable >= 3, f"{n_stable}"))

backbone_stable = stable_pathways[stable_pathways["pathway_category"] == "signaling_backbone"]
checks.append(("信号网络主干道排名稳定", len(backbone_stable) >= 1, f"{len(backbone_stable)}个"))

cross_stable = stable_pathways[stable_pathways["pathway_category"] == "cross_talk"]
checks.append(("交叉调控通路排名稳定", len(cross_stable) >= 1, f"{len(cross_stable)}个"))

print(f"\n{'检查项':<40} {'结果':<30} {'判定'}")
print("-" * 90)

pass_count = 0
for check_name, passed, actual in checks:
    status = "PASS" if passed else "FAIL"
    if passed:
        pass_count += 1
    print(f"{check_name:<40} {actual:<30} {status}")

print(f"\n通过检查: {pass_count}/{len(checks)}")

if pass_count >= 4:
    print("\n[判定] MVP v2 测试通过! 故事成立,建议全力推进A+E全量分析。")
elif pass_count >= 2:
    print("\n[判定] MVP v2 测试部分通过。故事基本成立,但可能需要调整参数网格或补充分析。")
else:
    print("\n[判定] MVP v2 测试未通过。建议重新评估策略。")

# 输出氧化应激相关靶点在Q1的详细情况
print("\n" + "=" * 70)
print("氧化应激相关靶点在Q1的情况")
print("=" * 70)
for t in ox_targets:
    row = target_stats_df[target_stats_df["gene"] == t]
    if len(row) > 0:
        r = row.iloc[0]
        in_q1 = "Q1" in r["quadrant"]
        print(f"  {t}: 频率={r['frequency']:.2f}, 平均Degree={r['avg_degree']:.1f}, {r['quadrant']}")

# 输出信号网络主干道通路排名
print("\n" + "=" * 70)
print("信号网络主干道通路排名 (Top 20)")
print("=" * 70)
backbone_pathways = rank_stats_df[rank_stats_df["pathway_category"] == "signaling_backbone"].sort_values("avg_rank")
for _, row in backbone_pathways.head(20).iterrows():
    print(f"  {row['pathway']:<50} 平均排名={row['avg_rank']:.1f}, 频率={row['frequency']:.2f}, EF={row['avg_ef']:.1f}")

# 输出交叉调控通路排名
print("\n" + "=" * 70)
print("交叉调控通路排名 (Top 20)")
print("=" * 70)
cross_pathways = rank_stats_df[rank_stats_df["pathway_category"] == "cross_talk"].sort_values("avg_rank")
for _, row in cross_pathways.head(20).iterrows():
    print(f"  {row['pathway']:<50} 平均排名={row['avg_rank']:.1f}, 频率={row['frequency']:.2f}, EF={row['avg_ef']:.1f}")

# 输出直接氧化应激响应通路（预期极少或空白）
print("\n" + "=" * 70)
print("直接氧化应激响应通路 (GO层级验证)")
print("=" * 70)
direct_ox_pathways = rank_stats_df[rank_stats_df["pathway_category"] == "direct_ox"].sort_values("avg_rank")
if len(direct_ox_pathways) > 0:
    for _, row in direct_ox_pathways.iterrows():
        print(f"  {row['pathway']:<50} 平均排名={row['avg_rank']:.1f}, 频率={row['frequency']:.2f}")
else:
    print("  [无] 在本数据集中，没有富集到 GO:0006979 后代的直接氧化应激响应通路")
    print("  说明：草莓叶活性成分可能通过间接调控信号网络来影响氧化应激平衡")

print("\n" + "=" * 70)
