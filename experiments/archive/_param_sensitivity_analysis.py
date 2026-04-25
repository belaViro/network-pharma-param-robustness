"""
参数敏感性分析 (A) + 不确定性量化 (E) = A+E 联合稳健性分析框架

核心思路：
1. A部分（参数敏感性）：系统扫描 SwissTarget Probability / STRING Score / P-value cutoff
   的参数空间，评估结果随参数的变化
2. E部分（不确定性量化）：对每个参数组合做 Bootstrap 重采样，评估结果的统计可靠性
3. A+E 结合：不仅回答"换参数结果怎么变"，还要回答"这个结果我有多大把握相信"

输出：
- 带 95% 置信区间的稳定性热图
- 靶点/通路出现频率的置信区间
- 综合稳健性评分（参数敏感性 × 统计可靠性）
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from itertools import product
from collections import Counter, defaultdict
import logging
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

from src.config import SWISS_TARGETS_DIR, ANTIOXIDANT_KEYWORDS, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ============================================================================
# 0. A+E 全局配置
# ============================================================================
N_BOOTSTRAP = 100          # Bootstrap 重采样次数（论文建议 ≥100，算力足可用 1000）
CI_ALPHA = 0.05            # 95% 置信区间
BOOTSTRAP_SEED = 42
np.random.seed(BOOTSTRAP_SEED)

RANDOM_STATE = np.random.RandomState(BOOTSTRAP_SEED)

# ============================================================================
# 1. 加载原始数据
# ============================================================================
print("=" * 70)
print("A+E 联合稳健性分析：参数敏感性 + 不确定性量化")
print("=" * 70)

# 加载所有 SwissTarget 预测结果（带 probability）
all_predictions = []
for csv_file in SWISS_TARGETS_DIR.glob("*.csv"):
    df = pd.read_csv(csv_file)
    if "Target" in df.columns and "Common name" in df.columns:
        if "Probability*" in df.columns:
            all_predictions.append(df[["Common name", "Probability*"]].rename(
                columns={"Common name": "Gene", "Probability*": "Probability"}
            ))
        elif "Probability" in df.columns:
            all_predictions.append(df[["Common name", "Probability"]].rename(
                columns={"Common name": "Gene"}
            ))
        elif "probability" in df.columns:
            all_predictions.append(df[["Common name", "probability"]].rename(
                columns={"Common name": "Gene", "probability": "Probability"}
            ))

if all_predictions:
    pred_df = pd.concat(all_predictions, ignore_index=True)
    pred_df = pred_df.drop_duplicates(subset=["Gene"])
    logger.info(f"加载 SwissTarget 预测结果: {len(pred_df)} 个唯一靶点 (基因符号)")
else:
    logger.error("未找到 SwissTarget 预测文件")
    sys.exit(1)

# 加载 STRING 边数据
string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
logger.info(f"加载 STRING 边数据: {len(string_edges)} 条边")

# ============================================================================
# 2. 定义参数网格
# ============================================================================
print("\n" + "=" * 70)
print("Step 1: 定义参数扫描空间 (A部分)")
print("=" * 70)

prob_thresholds = [0.0, 0.05, 0.1, 0.2, 0.3]
string_thresholds = [0.0, 0.3, 0.4, 0.5, 0.7]
pvalue_cutoffs = [0.01, 0.05]
top_k_values = [10, 20, 30, 50]

param_combinations = list(product(prob_thresholds, string_thresholds, pvalue_cutoffs, top_k_values))
logger.info(f"参数组合数: {len(param_combinations)} = "
            f"{len(prob_thresholds)} × {len(string_thresholds)} × {len(pvalue_cutoffs)} × {len(top_k_values)}")
logger.info(f"每个组合 Bootstrap {N_BOOTSTRAP} 次 → 总运行次数约: "
            f"{len(param_combinations)} × {N_BOOTSTRAP} = {len(param_combinations) * N_BOOTSTRAP}")
logger.info("注意: 200 组合 × 100 Bootstrap ≈ 20000 次分析，预计运行 2-6 小时")

# ============================================================================
# 3. 辅助函数：单次分析流程（用于 A 点估计和 E 每次 Bootstrap）
# ============================================================================
def run_single_analysis(filtered_targets, string_edges_df, string_th, pval_th, background, top_k=30):
    """
    对给定靶点列表执行：PPI构建 → Degree排序 → 富集分析 → 指标统计
    返回字典，包含所有关键指标。失败返回 None。
    """
    n_targets = len(filtered_targets)
    if n_targets < 10:
        return None

    # STRING score 过滤
    filtered_edges = string_edges_df[string_edges_df["combined_score"] >= string_th]
    G = nx.Graph()
    for _, row in filtered_edges.iterrows():
        G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

    valid_nodes = [n for n in G.nodes() if n in filtered_targets]
    G = G.subgraph(valid_nodes).copy()

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    if n_nodes < 10:
        return None

    # Degree 排序取 Top-k
    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    topk = [g for g, d in degree_sorted[:top_k]]
    if len(topk) < top_k:
        topk = [g for g, d in degree_sorted]

    # 富集分析
    try:
        enrichment = run_enrichment_analysis(
            gene_list=topk,
            pvalue_cutoff=pval_th,
            top_n=20,
            background=background,
        )
    except Exception:
        return None

    # 统计指标
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
# 4. 批量运行 A+E 分析
# ============================================================================
print("\n" + "=" * 70)
print("Step 2: A+E 联合分析（点估计 + Bootstrap 置信区间）")
print("=" * 70)

results_records = []
all_targets_by_param = {}      # 点估计的 Top-30 靶点
all_pathways_by_param = {}     # 点估计的 KEGG 通路
bootstrap_all = {}             # 每个参数的 Bootstrap 原始结果

for idx, (prob_th, string_th, pval_th, top_k) in enumerate(tqdm(param_combinations, desc="参数组合")):
    # ---- 4.1 A部分：Probability 阈值筛选靶点（点估计） ----
    filtered_targets = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
    n_targets = len(filtered_targets)

    if n_targets < 10:
        continue

    # ---- 4.2 A部分：单次分析（点估计） ----
    point_est = run_single_analysis(filtered_targets, string_edges, string_th, pval_th, BACKGROUND_GENES, top_k)
    if point_est is None:
        continue

    logger.info(f"\n[{idx+1}/{len(param_combinations)}] "
                f"Prob≥{prob_th}, STRING≥{string_th}, Pvalue<{pval_th}, TopK={top_k} | "
                f"靶点{n_targets} → PPI节点{point_est['n_ppi_nodes']} → "
                f"氧化应激通路{point_est['oxidative_pathways']}")

    # ---- 4.3 E部分：Bootstrap 重采样（不确定性量化） ----
    boot_records = []
    boot_target_presence = defaultdict(list)   # target -> [bool, bool, ...] 每次是否出现
    boot_pathway_presence = defaultdict(list)  # pathway -> [bool, bool, ...]

    # 预过滤 edges 加速
    filtered_edges = string_edges[string_edges["combined_score"] >= string_th]

    for b in range(N_BOOTSTRAP):
        # 从 filtered_targets 中有放回抽样
        boot_sample = RANDOM_STATE.choice(filtered_targets, size=len(filtered_targets), replace=True)
        boot_sample = list(set(boot_sample))  # 去重

        boot_result = run_single_analysis(boot_sample, filtered_edges, string_th, pval_th, BACKGROUND_GENES, top_k)

        if boot_result is None:
            continue

        boot_records.append({
            "oxidative_pathways": boot_result["oxidative_pathways"],
            "significant_pathways": boot_result["significant_pathways"],
            "avg_ef": boot_result["avg_ef"],
            "n_ppi_nodes": boot_result["n_ppi_nodes"],
        })

        # 记录靶点出现情况
        for t in boot_result["top_genes"]:
            boot_target_presence[t].append(1)
        # 记录未出现的靶点为 0（为了后续计算 CI，需要补齐）
        # 但更高效的方式是最后统一处理

        # 记录通路出现情况
        for p in boot_result["kegg_pathways"]:
            boot_pathway_presence[p].append(1)

    n_valid_boot = len(boot_records)
    if n_valid_boot < N_BOOTSTRAP * 0.5:
        logger.warning(f"  Bootstrap 成功率过低 ({n_valid_boot}/{N_BOOTSTRAP})，跳过该组合")
        continue

    # ---- 4.4 E部分：计算 95% 置信区间 ----
    ox_vals = [r["oxidative_pathways"] for r in boot_records]
    total_vals = [r["significant_pathways"] for r in boot_records]
    ef_vals = [r["avg_ef"] for r in boot_records]
    node_vals = [r["n_ppi_nodes"] for r in boot_records]

    ox_mean = np.mean(ox_vals)
    ox_ci_low, ox_ci_high = np.percentile(ox_vals, [2.5, 97.5])
    ox_std = np.std(ox_vals)

    total_mean = np.mean(total_vals)
    total_ci_low, total_ci_high = np.percentile(total_vals, [2.5, 97.5])

    ef_mean = np.mean(ef_vals)
    ef_ci_low, ef_ci_high = np.percentile(ef_vals, [2.5, 97.5])

    node_mean = np.mean(node_vals)
    node_ci_low, node_ci_high = np.percentile(node_vals, [2.5, 97.5])

    # CI 宽度（用于衡量不确定性大小）
    ox_ci_width = ox_ci_high - ox_ci_low

    logger.info(f"  Bootstrap: 氧化应激 {ox_mean:.1f} [{ox_ci_low:.1f}, {ox_ci_high:.1f}] "
                f"(CI宽度={ox_ci_width:.1f})")

    # ---- 4.5 记录 A+E 结果 ----
    record = {
        # 参数
        "prob_threshold": prob_th,
        "string_threshold": string_th,
        "pvalue_cutoff": pval_th,
        "top_k": top_k,
        # A部分：点估计
        "n_targets": point_est["n_targets"],
        "n_ppi_nodes": point_est["n_ppi_nodes"],
        "n_ppi_edges": point_est["n_ppi_edges"],
        "top_k": point_est["top_k"],
        "significant_pathways": point_est["significant_pathways"],
        "oxidative_pathways": point_est["oxidative_pathways"],
        "avg_ef": point_est["avg_ef"],
        "min_pvalue": point_est["min_pvalue"],
        "top_genes": ";".join(point_est["top_genes"]),
        # E部分：Bootstrap 汇总
        "n_valid_bootstrap": n_valid_boot,
        "ox_mean": ox_mean,
        "ox_ci_low": ox_ci_low,
        "ox_ci_high": ox_ci_high,
        "ox_ci_width": ox_ci_width,
        "ox_std": ox_std,
        "total_mean": total_mean,
        "total_ci_low": total_ci_low,
        "total_ci_high": total_ci_high,
        "ef_mean": ef_mean,
        "ef_ci_low": ef_ci_low,
        "ef_ci_high": ef_ci_high,
        "node_mean": node_mean,
        "node_ci_low": node_ci_low,
        "node_ci_high": node_ci_high,
    }
    results_records.append(record)
    param_key = (prob_th, string_th, pval_th, top_k)
    all_targets_by_param[param_key] = point_est["top_genes"]
    all_pathways_by_param[param_key] = point_est["kegg_pathways"]
    bootstrap_all[param_key] = {
        "boot_records": boot_records,
        "boot_target_presence": dict(boot_target_presence),
        "boot_pathway_presence": dict(boot_pathway_presence),
        "n_valid_boot": n_valid_boot,
    }

# ============================================================================
# 5. 结果汇总与分析（A+E）
# ============================================================================
print("\n" + "=" * 70)
print("Step 3: A+E 结果汇总分析")
print("=" * 70)

results_df = pd.DataFrame(results_records)
if len(results_df) == 0:
    logger.error("没有有效结果，请检查参数范围或数据")
    sys.exit(1)

# 保存原始结果
results_df.to_csv("src/outputs/sensitivity/param_sensitivity_results.csv", index=False)
logger.info(f"结果已保存: src/outputs/sensitivity/param_sensitivity_results.csv")

# ---- 5.1 参数组合效果排序（A+E 综合） ----
print("\n【A+E 参数组合排名】（按氧化应激均值降序，标注 CI）")
print("-" * 110)
print(f"{'Rank':<5} {'Prob':<6} {'STRING':<8} {'Pvalue':<8} {'TopK':<6} "
      f"{'氧化应激':<10} {'95% CI':<18} {'CI宽度':<8} "
      f"{'总通路':<10} {'PPI节点':<8} {'可靠性评级'}")
print("-" * 110)

sorted_results = results_df.sort_values("ox_mean", ascending=False)
for idx, (_, row) in enumerate(sorted_results.head(15).iterrows()):
    # 可靠性评级：均值高 + CI窄 = 高可靠性
    reliability = ""
    if row["ox_mean"] >= sorted_results["ox_mean"].quantile(0.75) and row["ox_ci_width"] <= 3:
        reliability = "★★★ 高可靠"
    elif row["ox_mean"] >= sorted_results["ox_mean"].median() and row["ox_ci_width"] <= 5:
        reliability = "★★ 中可靠"
    else:
        reliability = "★ 低可靠/高不确定"

    print(f"{idx+1:<5} {row['prob_threshold']:<6.2f} {row['string_threshold']:<8.1f} "
          f"{row['pvalue_cutoff']:<8.2f} {row['top_k']:<6.0f} {row['ox_mean']:<10.1f} "
          f"[{row['ox_ci_low']:.1f}, {row['ox_ci_high']:.1f}]{''.ljust(6)} "
          f"{row['ox_ci_width']:<8.1f} {row['total_mean']:<10.1f} "
          f"{row['node_mean']:<8.0f} {reliability}")

# ---- 5.2 最优参数窗口（A+E） ----
print("\n【A+E 最优参数窗口】")
best_ox = results_df.loc[results_df["ox_mean"].idxmax()]
print(f"氧化应激通路最多的参数组合（点估计）:")
print(f"  Probability ≥ {best_ox['prob_threshold']}, "
      f"STRING ≥ {best_ox['string_threshold']}, P-value < {best_ox['pvalue_cutoff']}, TopK={best_ox['top_k']:.0f}")
print(f"  Bootstrap 均值: {best_ox['ox_mean']:.1f} [{best_ox['ox_ci_low']:.1f}, {best_ox['ox_ci_high']:.1f}]")

# 最可靠参数（氧化应激高 + CI窄）
results_df["reliability_score"] = (
    results_df["ox_mean"] / (results_df["ox_mean"].max() + 1e-8) * 0.6 +
    (1 - results_df["ox_ci_width"] / (results_df["ox_ci_width"].max() + 1e-8)) * 0.4
)
most_reliable = results_df.loc[results_df["reliability_score"].idxmax()]
print(f"\n最可靠参数组合（高氧化应激 + 窄置信区间）:")
print(f"  Probability ≥ {most_reliable['prob_threshold']}, "
      f"STRING ≥ {most_reliable['string_threshold']}, P-value < {most_reliable['pvalue_cutoff']}, TopK={most_reliable['top_k']:.0f}")
print(f"  氧化应激: {most_reliable['ox_mean']:.1f} [{most_reliable['ox_ci_low']:.1f}, {most_reliable['ox_ci_high']:.1f}], "
      f"CI宽度={most_reliable['ox_ci_width']:.1f}")

# ============================================================================
# 6. 稳健靶点识别（含 Bootstrap 置信区间）
# ============================================================================
print("\n" + "=" * 70)
print("Step 4: 稳健靶点识别（A+E）")
print("=" * 70)

target_counter = Counter()
target_oxidative_count = defaultdict(int)

for param, targets in all_targets_by_param.items():
    ox_count = results_df[
        (results_df["prob_threshold"] == param[0]) &
        (results_df["string_threshold"] == param[1]) &
        (results_df["pvalue_cutoff"] == param[2]) &
        (results_df["top_k"] == param[3])
    ]["oxidative_pathways"].values
    has_ox = len(ox_count) > 0 and ox_count[0] > 0
    for t in targets:
        target_counter[t] += 1
        if has_ox:
            target_oxidative_count[t] += 1

n_valid_params = len(all_targets_by_param)

# 对每个靶点，利用 Bootstrap 结果计算出现频率的 CI
target_bootstrap_stats = []
for target, count in target_counter.most_common(50):
    freq = count / n_valid_params
    ox_freq = target_oxidative_count[target] / n_valid_params if target in target_oxidative_count else 0

    # 从所有参数的 Bootstrap 中收集该靶点的出现情况
    all_presence = []
    for param, boot_data in bootstrap_all.items():
        presence_dict = boot_data["boot_target_presence"]
        if target in presence_dict:
            all_presence.extend(presence_dict[target])
        else:
            # 该次 Bootstrap 未出现，补 0
            all_presence.extend([0] * boot_data["n_valid_boot"])

    if len(all_presence) > 10:
        boot_freq = np.mean(all_presence)
        boot_ci_low, boot_ci_high = np.percentile(all_presence, [2.5, 97.5])
    else:
        boot_freq = freq
        boot_ci_low, boot_ci_high = 0, 1

    robust_score = freq * 0.4 + ox_freq * 0.4 + (1 - (boot_ci_high - boot_ci_low)) * 0.2
    target_bootstrap_stats.append({
        "target": target,
        "count": count,
        "freq": freq,
        "ox_freq": ox_freq,
        "boot_freq": boot_freq,
        "boot_ci_low": boot_ci_low,
        "boot_ci_high": boot_ci_high,
        "robust_score": robust_score,
    })

target_bootstrap_df = pd.DataFrame(target_bootstrap_stats)
target_bootstrap_df = target_bootstrap_df.sort_values("robust_score", ascending=False)

print(f"\n有效参数组合数: {n_valid_params}")
print(f"\n{'靶点':<12} {'点估计频率':<10} {'Bootstrap频率':<12} {'95% CI':<16} {'稳健性评分'}")
print("-" * 65)

for _, row in target_bootstrap_df.head(20).iterrows():
    if row["freq"] >= 0.3:  # 只显示出现频率>30%的
        print(f"{row['target']:<12} {row['freq']:<10.2f} {row['boot_freq']:<12.2f} "
              f"[{row['boot_ci_low']:.2f}, {row['boot_ci_high']:.2f}]{''.ljust(6)} "
              f"{row['robust_score']:.3f}")

# 保存靶点稳健性结果
target_bootstrap_df.to_csv("src/outputs/sensitivity/robust_targets_with_ci.csv", index=False)

# ============================================================================
# 7. 核心通路识别（含 Bootstrap 置信区间）
# ============================================================================
print("\n" + "=" * 70)
print("Step 5: 核心通路识别（A+E）")
print("=" * 70)

pathway_counter = Counter()
for param, pathways in all_pathways_by_param.items():
    for p in pathways:
        pathway_counter[p] += 1

# 通路 Bootstrap CI
pathway_bootstrap_stats = []
for pathway, count in pathway_counter.most_common(30):
    freq = count / n_valid_params
    all_presence = []
    for param, boot_data in bootstrap_all.items():
        presence_dict = boot_data["boot_pathway_presence"]
        if pathway in presence_dict:
            all_presence.extend(presence_dict[pathway])
        else:
            all_presence.extend([0] * boot_data["n_valid_boot"])

    if len(all_presence) > 10:
        boot_freq = np.mean(all_presence)
        boot_ci_low, boot_ci_high = np.percentile(all_presence, [2.5, 97.5])
    else:
        boot_freq = freq
        boot_ci_low, boot_ci_high = 0, 1

    pathway_bootstrap_stats.append({
        "pathway": pathway,
        "count": count,
        "freq": freq,
        "boot_freq": boot_freq,
        "boot_ci_low": boot_ci_low,
        "boot_ci_high": boot_ci_high,
    })

pathway_bootstrap_df = pd.DataFrame(pathway_bootstrap_stats)
pathway_bootstrap_df = pathway_bootstrap_df.sort_values("boot_freq", ascending=False)

print(f"\n{'通路':<45} {'频率':<8} {'Bootstrap':<10} {'95% CI':<16} {'标记'}")
print("-" * 95)
for _, row in pathway_bootstrap_df.iterrows():
    marker = ""
    if row["boot_ci_high"] >= 0.8:
        marker = "★★★ 核心"
    elif row["boot_ci_high"] >= 0.5:
        marker = "★★ 重要"
    elif row["boot_ci_high"] >= 0.3:
        marker = "★ 潜在"
    else:
        marker = "  不确定"

    print(f"{row['pathway']:<43} {row['freq']:<8.2f} {row['boot_freq']:<10.2f} "
          f"[{row['boot_ci_low']:.2f}, {row['boot_ci_high']:.2f}]{''.ljust(6)} {marker}")

pathway_bootstrap_df.to_csv("src/outputs/sensitivity/core_pathways_with_ci.csv", index=False)

# ============================================================================
# 8. 可视化（A+E）
# ============================================================================
print("\n" + "=" * 70)
print("Step 6: A+E 可视化")
print("=" * 70)

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# 为热图选择参考子集: top_k=30（默认值），对 pvalue 取平均
ref_df = results_df[results_df["top_k"] == 30]

fig = plt.figure(figsize=(20, 16))

# ---- 8.1 热图1: 氧化应激通路数（点估计）----
ax1 = plt.subplot(3, 3, 1)
pivot_ox = ref_df.pivot_table(
    values="oxidative_pathways",
    index="prob_threshold",
    columns="string_threshold",
    aggfunc="mean"
)
im1 = ax1.imshow(pivot_ox.values, cmap="YlOrRd", aspect="auto", vmin=0)
ax1.set_xticks(range(len(pivot_ox.columns)))
ax1.set_xticklabels([f"{x:.1f}" for x in pivot_ox.columns])
ax1.set_yticks(range(len(pivot_ox.index)))
ax1.set_yticklabels([f"{x:.2f}" for x in pivot_ox.index])
ax1.set_xlabel("STRING Score Threshold")
ax1.set_ylabel("SwissTarget Probability Threshold")
ax1.set_title("A: Oxidative Stress Pathways\n(Point Estimate, TopK=30)")
plt.colorbar(im1, ax=ax1, fraction=0.046)

# ---- 8.2 热图2: Bootstrap 均值 ----
ax2 = plt.subplot(3, 3, 2)
pivot_ox_mean = ref_df.pivot_table(
    values="ox_mean",
    index="prob_threshold",
    columns="string_threshold",
    aggfunc="mean"
)
im2 = ax2.imshow(pivot_ox_mean.values, cmap="YlOrRd", aspect="auto", vmin=0)
ax2.set_xticks(range(len(pivot_ox_mean.columns)))
ax2.set_xticklabels([f"{x:.1f}" for x in pivot_ox_mean.columns])
ax2.set_yticks(range(len(pivot_ox_mean.index)))
ax2.set_yticklabels([f"{x:.2f}" for x in pivot_ox_mean.index])
ax2.set_xlabel("STRING Score Threshold")
ax2.set_ylabel("SwissTarget Probability Threshold")
ax2.set_title("E: Oxidative Stress Mean\n(Bootstrap, TopK=30)")
plt.colorbar(im2, ax=ax2, fraction=0.046)

# ---- 8.3 热图3: CI 宽度（不确定性大小） ----
ax3 = plt.subplot(3, 3, 3)
pivot_ox_width = ref_df.pivot_table(
    values="ox_ci_width",
    index="prob_threshold",
    columns="string_threshold",
    aggfunc="mean"
)
# CI 宽度越小越好，所以用反转色图
cmap_r = plt.cm.get_cmap("YlOrRd").reversed()
im3 = ax3.imshow(pivot_ox_width.values, cmap=cmap_r, aspect="auto", vmin=0)
ax3.set_xticks(range(len(pivot_ox_width.columns)))
ax3.set_xticklabels([f"{x:.1f}" for x in pivot_ox_width.columns])
ax3.set_yticks(range(len(pivot_ox_width.index)))
ax3.set_yticklabels([f"{x:.2f}" for x in pivot_ox_width.index])
ax3.set_xlabel("STRING Score Threshold")
ax3.set_ylabel("SwissTarget Probability Threshold")
ax3.set_title("E: Uncertainty (CI Width)\n(Smaller = More Reliable, TopK=30)")
plt.colorbar(im3, ax=ax3, fraction=0.046)

# ---- 8.4 热图4: 显著通路总数 ----
ax4 = plt.subplot(3, 3, 4)
pivot_total = ref_df.pivot_table(
    values="significant_pathways",
    index="prob_threshold",
    columns="string_threshold",
    aggfunc="mean"
)
im4 = ax4.imshow(pivot_total.values, cmap="Blues", aspect="auto")
ax4.set_xticks(range(len(pivot_total.columns)))
ax4.set_xticklabels([f"{x:.1f}" for x in pivot_total.columns])
ax4.set_yticks(range(len(pivot_total.index)))
ax4.set_yticklabels([f"{x:.2f}" for x in pivot_total.index])
ax4.set_xlabel("STRING Score Threshold")
ax4.set_ylabel("SwissTarget Probability Threshold")
ax4.set_title("A: Total Significant Pathways\n(TopK=30)")
plt.colorbar(im4, ax=ax4, fraction=0.046)

# ---- 8.5 热图5: Bootstrap 总通路均值 ----
ax5 = plt.subplot(3, 3, 5)
pivot_total_mean = ref_df.pivot_table(
    values="total_mean",
    index="prob_threshold",
    columns="string_threshold",
    aggfunc="mean"
)
im5 = ax5.imshow(pivot_total_mean.values, cmap="Blues", aspect="auto")
ax5.set_xticks(range(len(pivot_total_mean.columns)))
ax5.set_xticklabels([f"{x:.1f}" for x in pivot_total_mean.columns])
ax5.set_yticks(range(len(pivot_total_mean.index)))
ax5.set_yticklabels([f"{x:.2f}" for x in pivot_total_mean.index])
ax5.set_xlabel("STRING Score Threshold")
ax5.set_ylabel("SwissTarget Probability Threshold")
ax5.set_title("E: Total Pathways Mean\n(Bootstrap, TopK=30)")
plt.colorbar(im5, ax=ax5, fraction=0.046)

# ---- 8.6 热图6: 综合可靠性评分 ----
ax6 = plt.subplot(3, 3, 6)
pivot_rel = ref_df.pivot_table(
    values="reliability_score",
    index="prob_threshold",
    columns="string_threshold",
    aggfunc="mean"
)
im6 = ax6.imshow(pivot_rel.values, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
ax6.set_xticks(range(len(pivot_rel.columns)))
ax6.set_xticklabels([f"{x:.1f}" for x in pivot_rel.columns])
ax6.set_yticks(range(len(pivot_rel.index)))
ax6.set_yticklabels([f"{x:.2f}" for x in pivot_rel.index])
ax6.set_xlabel("STRING Score Threshold")
ax6.set_ylabel("SwissTarget Probability Threshold")
ax6.set_title("A+E: Reliability Score\n(High Ox + Low Uncertainty, TopK=30)")
plt.colorbar(im6, ax=ax6, fraction=0.046)

# ---- 8.7 Top-K 影响分析: 氧化应激均值随 TopK 变化 ----
ax7 = plt.subplot(3, 3, 7)
for pval in sorted(results_df["pvalue_cutoff"].unique()):
    sub = results_df[results_df["pvalue_cutoff"] == pval]
    topk_group = sub.groupby("top_k")["ox_mean"].agg(["mean", "std"]).reset_index()
    ax7.plot(topk_group["top_k"], topk_group["mean"], marker='o', linewidth=2,
             label=f"P-value < {pval}")
    ax7.fill_between(topk_group["top_k"],
                     topk_group["mean"] - topk_group["std"],
                     topk_group["mean"] + topk_group["std"],
                     alpha=0.2)
ax7.set_xlabel("Top-K Hub Genes")
ax7.set_ylabel("Bootstrap Mean: Oxidative Stress Pathways")
ax7.set_title("A+E: Top-K Impact on Results\n(Mean ± Std across all Prob/STRING)")
ax7.legend(loc="upper left")
ax7.grid(True, alpha=0.3)

# ---- 8.8 Top-K 影响分析: CI 宽度随 TopK 变化 ----
ax8 = plt.subplot(3, 3, 8)
for pval in sorted(results_df["pvalue_cutoff"].unique()):
    sub = results_df[results_df["pvalue_cutoff"] == pval]
    topk_group = sub.groupby("top_k")["ox_ci_width"].agg(["mean", "std"]).reset_index()
    ax8.plot(topk_group["top_k"], topk_group["mean"], marker='s', linewidth=2,
             label=f"P-value < {pval}")
    ax8.fill_between(topk_group["top_k"],
                     topk_group["mean"] - topk_group["std"],
                     topk_group["mean"] + topk_group["std"],
                     alpha=0.2)
ax8.set_xlabel("Top-K Hub Genes")
ax8.set_ylabel("Mean CI Width (Uncertainty)")
ax8.set_title("A+E: Top-K Impact on Uncertainty\n(Smaller = More Reliable)")
ax8.legend(loc="upper right")
ax8.grid(True, alpha=0.3)

# ---- 8.9 散点图: 氧化应激均值 vs CI 宽度（参数组合定位） ----
ax9 = plt.subplot(3, 3, 9)
# 用不同颜色区分不同 top_k
colors = {10: "#e74c3c", 20: "#3498db", 30: "#2ecc71", 50: "#9b59b6"}
for tk in sorted(results_df["top_k"].unique()):
    sub = results_df[results_df["top_k"] == tk]
    ax9.scatter(
        sub["ox_mean"],
        sub["ox_ci_width"],
        c=colors.get(tk, "gray"),
        s=80,
        edgecolors="black",
        linewidth=0.3,
        alpha=0.7,
        label=f"TopK={tk:.0f}"
    )
# 标注最可靠点
best_idx = results_df["reliability_score"].idxmax()
best_row = results_df.loc[best_idx]
ax9.annotate(
    f"Most Reliable\n({best_row['prob_threshold']:.2f}, {best_row['string_threshold']:.1f}, K={best_row['top_k']:.0f})",
    xy=(best_row["ox_mean"], best_row["ox_ci_width"]),
    xytext=(best_row["ox_mean"] + 1, best_row["ox_ci_width"] + 0.5),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="red"),
    color="red",
    fontweight="bold"
)
ax9.set_xlabel("Bootstrap Mean: Oxidative Stress Pathways")
ax9.set_ylabel("CI Width (Uncertainty)")
ax9.set_title("A+E: Parameter Landscape by Top-K")
ax9.legend(loc="upper right", fontsize=8)
# 添加理想区域标注
ax9.axhspan(0, 2, xmin=0.6, xmax=1.0, alpha=0.1, color="green")
ax9.text(
    results_df["ox_mean"].max() * 0.85, 1.0,
    "Ideal Zone\n(High Ox + Low Uncertainty)",
    fontsize=9, color="darkgreen", ha="center", style="italic"
)

plt.tight_layout()
plt.savefig("src/outputs/sensitivity/param_sensitivity_A_plus_E.png", dpi=200, bbox_inches="tight")
print(f"  A+E 综合热图已保存: src/outputs/sensitivity/param_sensitivity_A_plus_E.png")

# 额外保存一个高质量版本用于论文
plt.savefig("src/outputs/sensitivity/param_sensitivity_A_plus_E.pdf", bbox_inches="tight")
print(f"  PDF 版本已保存: src/outputs/sensitivity/param_sensitivity_A_plus_E.pdf")

# ============================================================================
# 9. 最终推荐（A+E）
# ============================================================================
print("\n" + "=" * 70)
print("Step 7: A+E 参数推荐")
print("=" * 70)

# 综合评分: 氧化应激均值(35%) + 显著通路均值(15%) + EF(15%) + 可靠性(35%)
results_df["ox_norm"] = (results_df["ox_mean"] - results_df["ox_mean"].min()) / \
                        (results_df["ox_mean"].max() - results_df["ox_mean"].min() + 1e-8)
results_df["total_norm"] = (results_df["total_mean"] - results_df["total_mean"].min()) / \
                           (results_df["total_mean"].max() - results_df["total_mean"].min() + 1e-8)
results_df["ef_norm"] = (results_df["ef_mean"] - results_df["ef_mean"].min()) / \
                        (results_df["ef_mean"].max() - results_df["ef_mean"].min() + 1e-8)
# 可靠性评分已计算

results_df["composite_score"] = (
    results_df["ox_norm"] * 0.35 +
    results_df["total_norm"] * 0.15 +
    results_df["ef_norm"] * 0.15 +
    results_df["reliability_score"] * 0.35
)

best_overall = results_df.loc[results_df["composite_score"].idxmax()]
print(f"\n【A+E 综合最优参数组合】")
print(f"  SwissTarget Probability ≥ {best_overall['prob_threshold']}")
print(f"  STRING Combined Score ≥ {best_overall['string_threshold']}")
print(f"  富集分析 P-value < {best_overall['pvalue_cutoff']}")
print(f"  Hub基因数 TopK = {best_overall['top_k']:.0f}")
print(f"  综合评分: {best_overall['composite_score']:.3f}")
print(f"  该组合下:")
print(f"    - 氧化应激通路: {best_overall['ox_mean']:.1f} [{best_overall['ox_ci_low']:.1f}, {best_overall['ox_ci_high']:.1f}] (点估计: {best_overall['oxidative_pathways']})")
print(f"    - 总通路: {best_overall['total_mean']:.1f} [{best_overall['total_ci_low']:.1f}, {best_overall['total_ci_high']:.1f}]")
print(f"    - PPI网络: {best_overall['node_mean']:.0f} [{best_overall['node_ci_low']:.0f}, {best_overall['node_ci_high']:.0f}] 节点")
print(f"    - CI宽度: {best_overall['ox_ci_width']:.1f} (不确定性指标)")

print(f"\n【A+E 稳健靶点 Top-10】（高频率 + 窄 CI）")
for _, row in target_bootstrap_df.head(10).iterrows():
    ci_marker = "✓" if (row["boot_ci_high"] - row["boot_ci_low"]) < 0.3 else "~"
    print(f"  {ci_marker} {row['target']}: Bootstrap频率 {row['boot_freq']:.1%} "
          f"[{row['boot_ci_low']:.2f}, {row['boot_ci_high']:.2f}], 评分 {row['robust_score']:.3f}")

print(f"\n【A+E 核心通路】（高频率 + 窄 CI）")
for _, row in pathway_bootstrap_df.head(10).iterrows():
    if row["boot_ci_high"] >= 0.3:
        ci_marker = "★" if row["boot_ci_high"] >= 0.7 else "☆"
        print(f"  {ci_marker} {row['pathway']}: Bootstrap频率 {row['boot_freq']:.1%} "
              f"[{row['boot_ci_low']:.2f}, {row['boot_ci_high']:.2f}]")

# 保存最终结果
results_df.to_csv("src/outputs/sensitivity/param_sensitivity_A_plus_E_results.csv", index=False)
print(f"\n完整结果已保存: src/outputs/sensitivity/param_sensitivity_A_plus_E_results.csv")

print("\n" + "=" * 70)
print("A+E 分析完成！")
print("=" * 70)
print("\n一句话总结：")
print("不仅回答'换参数结果怎么变'，还要回答'这个结果我有多大把握相信'")
print("=" * 70)
