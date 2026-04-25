"""
稳健靶点 vs 各参数组合的富集效果对比

比较Q1核心稳健靶点(33个)和Tier 1绝对稳健靶点(9个)的富集效果
是否优于每个参数组合的富集效果。

用法:
    cd code && python scripts/analysis/compare_robust_vs_parametric.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd
import json

from src.config import BACKGROUND_GENES
from src.sensitivity.enrichment_cache import cached_enrichment
from src.sensitivity.analysis import analyze_enrichment

# ============================================================================
# 加载各参数组合的富集结果
# ============================================================================
print("=" * 80)
print("稳健靶点 vs 各参数组合 富集效果对比")
print("=" * 80)

combo_df = pd.read_csv("src/outputs/sensitivity/param_grid_results.csv")
print(f"\n加载 {len(combo_df)} 个参数组合的结果")

print(f"\n--- 各参数组合富集指标统计 ---")
for idx, row in combo_df.iterrows():
    print(f"  {row['combo']:<25} "
          f"显著={row['significant_pathways']:3d} "
          f"主干={row['backbone_count']:2d} "
          f"交叉={row['cross_talk_count']:2d} "
          f"直接OX={row['direct_ox_count']:2d} "
          f"平均EF={row['avg_ef']:6.1f} "
          f"minP={row['min_pvalue']:.2e}")

print(f"\n--- 各参数组合指标汇总 ---")
for metric, label in [("significant_pathways", "总显著通路数"),
                       ("backbone_count", "信号主干道数"),
                       ("cross_talk_count", "交叉调控通路数"),
                       ("avg_ef", "平均EF"),
                       ("min_pvalue", "最小P-value")]:
    vals = combo_df[metric]
    print(f"  {label}: min={vals.min():.2f}, max={vals.max():.2f}, "
          f"mean={vals.mean():.2f}, median={vals.median():.2f}")

# ============================================================================
# 运行Q1/Tier1稳健靶点的富集分析
# ============================================================================
target_df = pd.read_csv("src/outputs/sensitivity/target_robustness.csv")
q1_genes = sorted(target_df[target_df["quadrant"] == "Q1: 核心稳健"]["gene"].tolist())
tier1_genes = sorted(target_df[target_df["frequency"] == 1.0]["gene"].tolist())

print(f"\n{'='*80}")
print(f"Q1 核心稳健靶点: {len(q1_genes)} 个")
print(f"  {', '.join(q1_genes)}")
print(f"Tier 1 绝对稳健靶点: {len(tier1_genes)} 个")
print(f"  {', '.join(tier1_genes)}")

PVALUE_CUTOFF = 0.01
results = {}

for label, genes in [("Q1_核心稳健(33个)", q1_genes), ("Tier1_绝对稳健(9个)", tier1_genes)]:
    print(f"\n>> 运行 {label} 的富集分析...")
    enrichment = cached_enrichment(genes, PVALUE_CUTOFF, BACKGROUND_GENES)
    metrics = analyze_enrichment(enrichment, include_sources=["KEGG", "GO_BP"])
    results[label] = metrics

    print(f"  显著通路总数: {metrics['total_significant']}")
    print(f"  信号主干道: {metrics['backbone']}")
    print(f"  交叉调控: {metrics['cross_talk']}")
    print(f"  直接OX: {metrics['direct_ox']}")
    print(f"  平均EF: {metrics['avg_ef']}")
    print(f"  最小P-value: {metrics['min_pvalue']:.2e}")

    print(f"\n  富集到的具体通路 (按分类):")
    for cat_name in ["signaling_backbone", "cross_talk", "direct_ox", "other"]:
        cat_pathways = [p for p in metrics["pathway_details"] if p["pathway_category"] == cat_name]
        if cat_pathways:
            cat_label = {"signaling_backbone": "信号主干道", "cross_talk": "交叉调控",
                         "direct_ox": "直接OX", "other": "其他"}[cat_name]
            print(f"  [{cat_label}] {len(cat_pathways)} 条:")
            for p in sorted(cat_pathways, key=lambda x: x["pvalue"])[:10]:
                print(f"    - {p['pathway'][:60]:<60} EF={p['ef']:.1f} P={p['pvalue']:.2e}")

# ============================================================================
# 对比分析
# ============================================================================
print(f"\n{'='*80}")
print("对比分析：稳健靶点富集效果是否优于各参数组合？")
print("=" * 80)

comparisons = [
    ("总显著通路数 (越多越好)", "significant_pathways", "higher"),
    ("信号主干道通路数 (越多越好)", "backbone_count", "higher"),
    ("交叉调控通路数 (越多越好)", "cross_talk_count", "higher"),
    ("平均Enrichment Factor (越高越好)", "avg_ef", "higher"),
    ("最小P-value (越小越好)", "min_pvalue", "lower"),
]

metric_map = {
    "significant_pathways": "total_significant",
    "backbone_count": "backbone",
    "cross_talk_count": "cross_talk",
    "avg_ef": "avg_ef",
    "min_pvalue": "min_pvalue",
}

for label, metric, direction in comparisons:
    combo_vals = combo_df[metric]
    print(f"\n--- {label} ---")
    print(f"  参数组合范围: {combo_vals.min():.4f} ~ {combo_vals.max():.4f} "
          f"(均值 {combo_vals.mean():.4f})")

    for rlabel, rmetrics in results.items():
        rval = rmetrics.get(metric_map[metric], 0)
        if direction == "higher":
            rank = int((combo_vals < rval).sum())
            better = rval >= combo_vals.max()
        else:
            rank = int((combo_vals > rval).sum())
            better = rval <= combo_vals.min()

        percentile = rank / len(combo_df) * 100
        status = "优于全部组合" if better else "未优于全部组合"
        print(f"  {rlabel}: {rval:.4f} "
              f"(优于 {rank}/{len(combo_df)} 个组合 = {percentile:.0f}% 分位) [{status}]")

# ============================================================================
# 综合判定
# ============================================================================
print(f"\n{'='*80}")
print("综合判定")
print("=" * 80)

q1_metrics = results["Q1_核心稳健(33个)"]
t1_metrics = results["Tier1_绝对稳健(9个)"]

def check(val, better, combo_vals, metric_name, higher_better=True):
    if higher_better:
        max_val = combo_vals.max()
        flag = "+" if better else "-"
        return f"  [{flag}] {metric_name}: {val:.2f} vs 最佳组合 {max_val:.2f}"
    else:
        min_val = combo_vals.min()
        flag = "+" if better else "-"
        return f"  [{flag}] {metric_name}: {val:.2e} vs 最佳组合 {min_val:.2e}"

# Q1判定
q1_total_better = q1_metrics["total_significant"] >= combo_df["significant_pathways"].max()
q1_backbone_better = q1_metrics["backbone"] >= combo_df["backbone_count"].max()
q1_ef_better = q1_metrics["avg_ef"] >= combo_df["avg_ef"].max()
q1_pval_better = q1_metrics["min_pvalue"] <= combo_df["min_pvalue"].min()
q1_score = sum([q1_total_better, q1_backbone_better, q1_ef_better, q1_pval_better])

print(f"\nQ1 核心稳健靶点 ({len(q1_genes)}个):")
print(check(q1_metrics['total_significant'], q1_total_better, combo_df['significant_pathways'], '总显著通路', True))
print(check(q1_metrics['backbone'], q1_backbone_better, combo_df['backbone_count'], '主干道数', True))
print(check(q1_metrics['avg_ef'], q1_ef_better, combo_df['avg_ef'], '平均EF', True))
print(check(q1_metrics['min_pvalue'], q1_pval_better, combo_df['min_pvalue'], '最小P-value', False))
print(f"  => 通过 {q1_score}/4 项")

# Tier1判定
t1_total_better = t1_metrics["total_significant"] >= combo_df["significant_pathways"].max()
t1_backbone_better = t1_metrics["backbone"] >= combo_df["backbone_count"].max()
t1_ef_better = t1_metrics["avg_ef"] >= combo_df["avg_ef"].max()
t1_pval_better = t1_metrics["min_pvalue"] <= combo_df["min_pvalue"].min()
t1_score = sum([t1_total_better, t1_backbone_better, t1_ef_better, t1_pval_better])

print(f"\nTier 1 绝对稳健靶点 ({len(tier1_genes)}个):")
print(check(t1_metrics['total_significant'], t1_total_better, combo_df['significant_pathways'], '总显著通路', True))
print(check(t1_metrics['backbone'], t1_backbone_better, combo_df['backbone_count'], '主干道数', True))
print(check(t1_metrics['avg_ef'], t1_ef_better, combo_df['avg_ef'], '平均EF', True))
print(check(t1_metrics['min_pvalue'], t1_pval_better, combo_df['min_pvalue'], '最小P-value', False))
print(f"  => 通过 {t1_score}/4 项")

# 整体结论
print(f"\n{'='*80}")
if t1_score >= 3 or q1_score >= 3:
    print("[结论] ++ 稳健靶点的富集效果确实优于或等于绝大部分参数组合")
    print("  说明：通过稳健性筛选的靶点集合，其生物学信号比任何单一参数")
    print("  组合的Top-K筛选结果都更强、更可靠。")
else:
    print("[结论] -- 稳健靶点的富集效果并非全面优于各参数组合")
    print("  这说明稳健性筛选在保持富集质量的同时，可能牺牲了部分信号强度。")
    print("  但稳健靶点的核心优势在于可重复性，而非极值表现。")
print("=" * 80)

# 保存对比结果
comparison_summary = {
    "dataset": "稳健靶点 vs 各参数组合 富集效果对比",
    "n_param_combos": len(combo_df),
    "n_q1_targets": len(q1_genes),
    "n_tier1_targets": len(tier1_genes),
    "combo_stats": {
        "significant_pathways": {"min": int(combo_df["significant_pathways"].min()),
                                 "max": int(combo_df["significant_pathways"].max()),
                                 "mean": round(float(combo_df["significant_pathways"].mean()), 1)},
        "backbone_count": {"min": int(combo_df["backbone_count"].min()),
                           "max": int(combo_df["backbone_count"].max()),
                           "mean": round(float(combo_df["backbone_count"].mean()), 1)},
        "avg_ef": {"min": round(float(combo_df["avg_ef"].min()), 1),
                   "max": round(float(combo_df["avg_ef"].max()), 1),
                   "mean": round(float(combo_df["avg_ef"].mean()), 1)},
    },
    "q1_metrics": {k: v for k, v in q1_metrics.items() if k != "pathway_details"},
    "tier1_metrics": {k: v for k, v in t1_metrics.items() if k != "pathway_details"},
    "q1_beats_all": {
        "total_significant": q1_total_better, "backbone": q1_backbone_better,
        "avg_ef": q1_ef_better, "min_pvalue": q1_pval_better,
    },
    "tier1_beats_all": {
        "total_significant": t1_total_better, "backbone": t1_backbone_better,
        "avg_ef": t1_ef_better, "min_pvalue": t1_pval_better,
    },
}

output_path = "src/outputs/sensitivity/robust_vs_param_comparison.json"
with open(output_path, "w") as f:
    json.dump(comparison_summary, f, indent=2, default=str)
print(f"\n对比结果已保存: {output_path}")
