"""
A+E 策略有效性快速验证

目标：用最小计算成本验证三件事：
1. 换参数结果是否真的变？（A 有效性）
2. Bootstrap 重采样是否真的能反映不确定性？（E 有效性）
3. A+E 结合是否能识别出"高均值+低不确定性"的最优参数？（A+E 价值）

策略：
- 只跑 6 个代表性参数组合（而非全部 25 个）
- 每个组合只做 20 次 Bootstrap（而非 100 次）
- 总计算量：6 × 20 = 120 次富集分析，预计 10-20 分钟
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter, defaultdict
import logging
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

from src.config import SWISS_TARGETS_DIR, ANTIOXIDANT_KEYWORDS, BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ========================================================================
# 配置：小规模快速验证
# ========================================================================
N_BOOTSTRAP = 5           # 快速验证用 5 次
RANDOM_STATE = np.random.RandomState(42)

print("=" * 70)
print("A+E Strategy Validity Check")
print("=" * 70)
print(f"Bootstrap: {N_BOOTSTRAP} (quick mode)")
print(f"Param combos: 6 (representative)")
print("=" * 70)

# ========================================================================
# 1. 加载数据
# ========================================================================
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

pred_df = pd.concat(all_predictions, ignore_index=True).drop_duplicates(subset=["Gene"])
logger.info(f"加载 SwissTarget: {len(pred_df)} 个唯一靶点 (基因符号)")

string_edges = pd.read_csv("src/outputs/shared/ppi_network/string_interactions_interactions.csv")
logger.info(f"加载 STRING: {len(string_edges)} 条边")

# ========================================================================
# 2. 选取 6 个代表性参数组合（覆盖参数空间的关键角落）
# ========================================================================
# 组合设计：
# - 宽松 vs 严格（Prob + STRING 的极端组合）
# - 保守 vs 宽松（Pvalue 0.01 vs 0.05）
param_combinations = [
    (0.0, 0.0, 0.05),   # 最宽松：所有靶点都进，低 STRING，宽松 Pvalue
    (0.0, 0.7, 0.05),   # 宽松 Prob + 高 STRING
    (0.1, 0.4, 0.05),   # 中等参数（接近经验默认值）
    (0.2, 0.5, 0.05),   # 较严格
    (0.3, 0.7, 0.05),   # 最严格
    (0.1, 0.4, 0.01),   # 同 (0.1, 0.4) 但 Pvalue 更严格
]

# ========================================================================
# 3. 单次分析函数
# ========================================================================
def run_single(filtered_targets, string_edges_df, string_th, pval_th):
    """单次分析：PPI → Degree Top-30 → 富集"""
    if len(filtered_targets) < 10:
        return None

    edges = string_edges_df[string_edges_df["combined_score"] >= string_th]
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

    valid_nodes = [n for n in G.nodes() if n in filtered_targets]
    G = G.subgraph(valid_nodes).copy()

    if G.number_of_nodes() < 10:
        return None

    degree_sorted = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    top30 = [g for g, d in degree_sorted[:30]]
    if len(top30) < 30:
        top30 = [g for g, d in degree_sorted]

    try:
        enrichment = run_enrichment_analysis(
            gene_list=top30,
            pvalue_cutoff=pval_th,
            top_n=20,
            background=BACKGROUND_GENES,
        )
    except Exception:
        return None

    ox_count = 0
    total_sig = 0
    if enrichment["status"] == "成功" and "significant" in enrichment:
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                total_sig += len(df_sig)
                for _, row_sig in df_sig.iterrows():
                    term = str(row_sig.get("Term", "")).lower()
                    if any(kw.lower() in term for kw in ANTIOXIDANT_KEYWORDS):
                        ox_count += 1

    return {
        "n_targets": len(filtered_targets),
        "n_ppi_nodes": G.number_of_nodes(),
        "top_k": len(top30),
        "oxidative_pathways": ox_count,
        "significant_pathways": total_sig,
        "top_genes": top30,
    }

# ========================================================================
# 4. 主循环：对每个参数组合做点估计 + Bootstrap
# ========================================================================
results = []

for prob_th, string_th, pval_th in param_combinations:
    print(f"\n{'-' * 60}")
    print(f"Params: Prob>={prob_th}, STRING>={string_th}, Pvalue<{pval_th}")
    print(f"{'-' * 60}")

    # A 部分：点估计
    filtered = pred_df[pred_df["Probability"] >= prob_th]["Gene"].tolist()
    point = run_single(filtered, string_edges, string_th, pval_th)

    if point is None:
        print(f"  [WARN] Point estimate failed (too few targets/network), skip")
        continue

    print(f"  [A] Point estimate: targets={point['n_targets']} -> PPI_nodes={point['n_ppi_nodes']} -> "
          f"ox_pathways={point['oxidative_pathways']}, total_pathways={point['significant_pathways']}")

    # E 部分：Bootstrap
    print(f"  [E] Bootstrapping {N_BOOTSTRAP} times...", end=" ", flush=True)
    boot_ox = []
    boot_total = []
    boot_genes = []

    for b in range(N_BOOTSTRAP):
        boot_sample = list(set(RANDOM_STATE.choice(filtered, size=len(filtered), replace=True)))
        boot_result = run_single(boot_sample, string_edges, string_th, pval_th)
        if boot_result:
            boot_ox.append(boot_result["oxidative_pathways"])
            boot_total.append(boot_result["significant_pathways"])
            boot_genes.append(boot_result["top_genes"])

    print(f"success {len(boot_ox)}/{N_BOOTSTRAP}")

    if len(boot_ox) < N_BOOTSTRAP * 0.5:
        print(f"  [WARN] Bootstrap success rate too low, skip")
        continue

    ox_mean = np.mean(boot_ox)
    ox_std = np.std(boot_ox)
    ox_ci_low, ox_ci_high = np.percentile(boot_ox, [2.5, 97.5])
    ox_ci_width = ox_ci_high - ox_ci_low

    total_mean = np.mean(boot_total)
    total_ci_low, total_ci_high = np.percentile(boot_total, [2.5, 97.5])

    # 靶点频率 Bootstrap
    gene_counter = Counter()
    for genes in boot_genes:
        for g in genes:
            gene_counter[g] += 1
    max_freq = max(gene_counter.values()) / len(boot_genes) if gene_counter else 0
    stable_genes = [g for g, c in gene_counter.items() if c / len(boot_genes) >= 0.7]

    print(f"  [E] Ox pathways Bootstrap: mean={ox_mean:.1f}, std={ox_std:.1f}, "
          f"95% CI=[{ox_ci_low:.1f}, {ox_ci_high:.1f}], CI_width={ox_ci_width:.1f}")
    print(f"  [E] Stable targets (freq>=70%): {len(stable_genes)}")
    if stable_genes:
        print(f"       {', '.join(stable_genes[:8])}{'...' if len(stable_genes) > 8 else ''}")

    results.append({
        "prob": prob_th,
        "string": string_th,
        "pval": pval_th,
        "point_ox": point["oxidative_pathways"],
        "point_total": point["significant_pathways"],
        "boot_ox_mean": ox_mean,
        "boot_ox_std": ox_std,
        "boot_ox_ci_low": ox_ci_low,
        "boot_ox_ci_high": ox_ci_high,
        "boot_ox_ci_width": ox_ci_width,
        "boot_total_mean": total_mean,
        "boot_total_ci_low": total_ci_low,
        "boot_total_ci_high": total_ci_high,
        "n_boot": len(boot_ox),
        "n_stable_genes": len(stable_genes),
        "stable_genes": ";".join(stable_genes),
    })

# ========================================================================
# 5. 综合评估：A+E 策略是否有效？
# ========================================================================
print("\n" + "=" * 70)
print("A+E Strategy Validity Report")
print("=" * 70)

if len(results) < 3:
    print("[FAIL] Valid param combos < 3, A+E cannot be validated on current data")
    print("   Possible reasons: param range too strict, or data too small")
    sys.exit(1)

df = pd.DataFrame(results)

# -- Metric 1: A validity (param sensitivity) --
print("\n[Metric 1] Param Sensitivity (A validity)")
ox_range = df["point_ox"].max() - df["point_ox"].min()
ox_cv = df["point_ox"].std() / (df["point_ox"].mean() + 1e-8)
print(f"  Ox pathway range: {df['point_ox'].min():.0f} ~ {df['point_ox'].max():.0f} (span={ox_range:.0f})")
print(f"  Coefficient of Variation: {ox_cv:.2%}")

if ox_range >= 2:
    print(f"  [PASS] Conclusion: param changes cause obvious result changes, A is VALID")
elif ox_range >= 1:
    print(f"  [WARN] Conclusion: param changes cause slight changes, A is WEAKLY valid")
else:
    print(f"  [FAIL] Conclusion: param changes have almost no effect, A is INVALID")

# -- Metric 2: E validity (uncertainty quantification) --
print("\n[Metric 2] Uncertainty Quantification (E validity)")
avg_ci_width = df["boot_ox_ci_width"].mean()
max_ci_width = df["boot_ox_ci_width"].max()
print(f"  Average CI width: {avg_ci_width:.1f}")
print(f"  Max CI width: {max_ci_width:.1f}")

high_uncertainty = df[df["boot_ox_ci_width"] >= 3]
if len(high_uncertainty) > 0:
    print(f"  High uncertainty combos: {len(high_uncertainty)}/{len(df)}")
    print(f"     These combos have large ox pathway fluctuations, interpret with caution")

if avg_ci_width >= 1.5:
    print(f"  [PASS] Conclusion: Bootstrap CI width reasonable, E is VALID")
elif avg_ci_width >= 0.5:
    print(f"  [WARN] Conclusion: Bootstrap CI width small, E is WEAKLY valid")
else:
    print(f"  [FAIL] Conclusion: Bootstrap CI width near zero, E is INVALID")

# -- Metric 3: A+E combined value --
print("\n[Metric 3] A+E Combined Value (can identify best params?)")
df["reliability"] = (df["boot_ox_mean"] / (df["boot_ox_mean"].max() + 1e-8)) * 0.6 + \
                     (1 - df["boot_ox_ci_width"] / (df["boot_ox_ci_width"].max() + 1e-8)) * 0.4

best = df.loc[df["reliability"].idxmax()]
worst = df.loc[df["reliability"].idxmin()]

print(f"  Most reliable params: Prob>={best['prob']}, STRING>={best['string']}, Pvalue<{best['pval']}")
print(f"     -> ox mean={best['boot_ox_mean']:.1f}, CI=[{best['boot_ox_ci_low']:.1f}, {best['boot_ox_ci_high']:.1f}], width={best['boot_ox_ci_width']:.1f}")
print(f"  Least reliable params: Prob>={worst['prob']}, STRING>={worst['string']}, Pvalue<{worst['pval']}")
print(f"     -> ox mean={worst['boot_ox_mean']:.1f}, CI=[{worst['boot_ox_ci_low']:.1f}, {worst['boot_ox_ci_high']:.1f}], width={worst['boot_ox_ci_width']:.1f}")

if best.name != worst.name and best["boot_ox_mean"] > worst["boot_ox_mean"] * 1.1 or best["boot_ox_ci_width"] < worst["boot_ox_ci_width"] * 0.8:
    print(f"  [PASS] Conclusion: A+E can distinguish high/low reliability, COMBO is VALUABLE")
else:
    print(f"  [WARN] Conclusion: A+E distinction limited, try wider param range")

# -- Metric 4: Pvalue cutoff impact --
print("\n[Metric 4] Pvalue cutoff impact")
pval05 = df[df["pval"] == 0.05]
pval01 = df[df["pval"] == 0.01]
if len(pval05) > 0 and len(pval01) > 0:
    # 找相同 prob+string 不同 pval 的组合
    merged = pval05.merge(pval01, on=["prob", "string"], suffixes=("_05", "_01"))
    if len(merged) > 0:
        row = merged.iloc[0]
        print(f"  Same params (Prob>={row['prob']}, STRING>={row['string']}):")
        print(f"     P<0.05: ox={row['point_ox_05']:.0f}, total={row['point_total_05']:.0f}")
        print(f"     P<0.01: ox={row['point_ox_01']:.0f}, total={row['point_total_01']:.0f}")
        if row['point_ox_05'] != row['point_ox_01']:
            print(f"  [PASS] Pvalue cutoff has obvious impact")
        else:
            print(f"  [WARN] Pvalue cutoff has no impact (same pathway count)")

# ========================================================================
# 6. 最终裁决
# ========================================================================
print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)

score_a = 1 if ox_range >= 2 else (0.5 if ox_range >= 1 else 0)
score_e = 1 if avg_ci_width >= 1.5 else (0.5 if avg_ci_width >= 0.5 else 0)
score_combo = 1 if (best.name != worst.name and (best["boot_ox_mean"] > worst["boot_ox_mean"] * 1.1 or best["boot_ox_ci_width"] < worst["boot_ox_ci_width"] * 0.8)) else 0.5
total_score = score_a + score_e + score_combo

print(f"\nA validity score: {score_a}/1.0")
print(f"E validity score: {score_e}/1.0")
print(f"A+E combo score: {score_combo}/1.0")
print(f"Total score: {total_score}/3.0")

if total_score >= 2.5:
    verdict = "[PASS] A+E is HIGHLY EFFECTIVE, proceed to full run immediately"
elif total_score >= 1.5:
    verdict = "[WARN] A+E is MODERATELY effective, tune param range then full run"
else:
    verdict = "[FAIL] A+E is NOT effective under current settings, redesign params"

print(f"\n{verdict}")

# Save results
output_path = Path("src/outputs/gnn/A_plus_E_validity_check.csv")
output_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_path, index=False)
print(f"\nDetailed results saved: {output_path}")
print("=" * 70)
