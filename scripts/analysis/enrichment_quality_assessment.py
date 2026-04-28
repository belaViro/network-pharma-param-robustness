#!/usr/bin/env python3
"""
富集质量综合评估：比较 Q1 靶点富集通路 vs 全通路池的质量差异

用法:
    cd code && python scripts/analysis/enrichment_quality_assessment.py
"""

import sys
from pathlib import Path
import json

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd
import numpy as np
from scipy import stats

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
PATHWAY_CSV = OUTPUTS_DIR / "pathway_rank_stability_augmented.csv"
AUDIT_JSON = OUTPUTS_DIR / "q1_enrichment_audit.json"
OUTPUT_JSON = OUTPUTS_DIR / "enrichment_quality_assessment.json"

# ── Load data ──────────────────────────────────────────────────────────
pathway_df = pd.read_csv(PATHWAY_CSV)
audit_data = json.load(open(AUDIT_JSON, 'r', encoding='utf-8'))

all_pathways = pathway_df.copy()

# Q1靶点富集到的系统通路 (KEGG + GO_BP = 40条)
enriched_pathways = []
for source in ("KEGG", "GO_BP"):
    if source in audit_data['pathway_details']:
        for p in audit_data['pathway_details'][source]:
            enriched_pathways.append(p['pathway'])

# 在 pathway_df 中标记哪些是被富集到的
all_pathways['enriched'] = all_pathways['pathway'].isin(enriched_pathways)

enriched_df = all_pathways[all_pathways['enriched'] == True]
non_enriched_df = all_pathways[all_pathways['enriched'] == False]

print("=" * 70)
print("富集质量综合评估")
print("=" * 70)
print(f"\n全通路池: {len(all_pathways)} 条")
print(f"Q1靶点富集通路: {len(enriched_df)} 条")
print(f"未被富集通路: {len(non_enriched_df)} 条")

# ── 核心质量指标对比 ───────────────────────────────────────────────────
metrics = ['frequency', 'avg_logp', 'avg_rank', 'rank_std']

print(f"\n{'='*70}")
print("核心质量指标对比 (平均值 ± 标准差)")
print(f"{'='*70}")
all_label = f"全通路池 ({len(all_pathways)}条)"
enriched_label = f"富集通路 ({len(enriched_df)}条)"
print(f"{'指标':<20} {all_label:<22} {enriched_label:<22} {'差异':<10}")
print("-" * 70)

comparison_results = {}

for metric in metrics:
    all_mean = all_pathways[metric].mean()
    all_std = all_pathways[metric].std()
    enriched_mean = enriched_df[metric].mean()
    enriched_std = enriched_df[metric].std()
    diff_pct = (enriched_mean - all_mean) / all_mean * 100 if all_mean != 0 else 0
    
    # 统计检验 (Mann-Whitney U，因为数据可能不正态)
    stat, pvalue = stats.mannwhitneyu(
        enriched_df[metric].dropna(), 
        non_enriched_df[metric].dropna(),
        alternative='two-sided'
    )
    
    comparison_results[metric] = {
        'all_mean': float(all_mean),
        'all_std': float(all_std),
        'enriched_mean': float(enriched_mean),
        'enriched_std': float(enriched_std),
        'diff_pct': float(diff_pct),
        'mann_whitney_p': float(pvalue),
    }
    
    sig = "***" if pvalue < 0.001 else "**" if pvalue < 0.01 else "*" if pvalue < 0.05 else "ns"
    print(f"{metric:<20} {all_mean:>8.3f} ± {all_std:<6.3f}    {enriched_mean:>8.3f} ± {enriched_std:<6.3f}    {diff_pct:>+6.1f}% {sig}")

print("\n注: *** p<0.001, ** p<0.01, * p<0.05, ns 不显著")

# ── 中位数对比 ─────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("中位数对比")
print(f"{'='*70}")
print(f"{'指标':<20} {'全通路池':<15} {'富集通路':<15} {'差异':<10}")
print("-" * 55)

for metric in metrics:
    all_median = all_pathways[metric].median()
    enriched_median = enriched_df[metric].median()
    diff_pct = (enriched_median - all_median) / all_median * 100 if all_median != 0 else 0
    print(f"{metric:<20} {all_median:>10.3f}    {enriched_median:>10.3f}    {diff_pct:>+6.1f}%")

# ── 加权稳健率 ─────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("加权稳健率 (替代简单的计数比例)")
print(f"{'='*70}")

# 简单计数法
q1_count = sum(1 for p in enriched_df['robustness_quadrant'] if 'Q1' in p)
simple_q1_rate = q1_count / len(enriched_df) * 100

# 按 frequency 加权
freq_weighted = enriched_df.copy()
freq_weighted['q1_weight'] = freq_weighted['frequency'] * freq_weighted['robustness_quadrant'].apply(lambda x: 1.0 if 'Q1' in x else 0.0)
freq_weighted_q1_rate = freq_weighted['q1_weight'].sum() / freq_weighted['frequency'].sum() * 100

# 按 avg_logp 加权
logp_weighted = enriched_df.copy()
logp_weighted['q1_weight'] = logp_weighted['avg_logp'] * logp_weighted['robustness_quadrant'].apply(lambda x: 1.0 if 'Q1' in x else 0.0)
logp_weighted_q1_rate = logp_weighted['q1_weight'].sum() / logp_weighted['avg_logp'].sum() * 100

# 按 frequency × avg_logp 加权
combined_weighted = enriched_df.copy()
combined_weighted['weight'] = combined_weighted['frequency'] * combined_weighted['avg_logp']
combined_weighted['q1_weight'] = combined_weighted['weight'] * combined_weighted['robustness_quadrant'].apply(lambda x: 1.0 if 'Q1' in x else 0.0)
combined_weighted_q1_rate = combined_weighted['q1_weight'].sum() / combined_weighted['weight'].sum() * 100

print(f"  简单计数 Q1 率:            {simple_q1_rate:.1f}%")
print(f"  按 frequency 加权 Q1 率:   {freq_weighted_q1_rate:.1f}%")
print(f"  按 avg_logP 加权 Q1 率:    {logp_weighted_q1_rate:.1f}%")
print(f"  按 freq×logP 加权 Q1 率:   {combined_weighted_q1_rate:.1f}%")

# ── 分层质量分析 ───────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("富集通路的分层质量")
print(f"{'='*70}")

for quadrant in ['Q1: 候选稳健通路', 'Q2: 低频高显著', 'Q3: 稳健但效应弱', 'Q4: 可忽略']:
    subset = enriched_df[enriched_df['robustness_quadrant'] == quadrant]
    if len(subset) == 0:
        continue
    avg_freq = subset['frequency'].mean()
    avg_logp = subset['avg_logp'].mean()
    avg_rank = subset['avg_rank'].mean()
    print(f"\n  {quadrant} ({len(subset)}条):")
    print(f"    平均频率: {avg_freq:.3f} | 平均 avg_logP: {avg_logp:.3f} | 平均排名: {avg_rank:.1f}")

# ── 关键结论 ───────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("核心结论")
print(f"{'='*70}")

freq_diff = comparison_results['frequency']['diff_pct']
logp_diff = comparison_results['avg_logp']['diff_pct']
rank_diff = comparison_results['avg_rank']['diff_pct']

print(f"\n  Q1 靶点富集通路的平均频率比全池高 {freq_diff:+.1f}%")
print(f"  Q1 靶点富集通路的平均 avg_logP 比全池高 {logp_diff:+.1f}%")
print(f"  Q1 靶点富集通路的平均排名比全池低 {rank_diff:+.1f}% (排名越低越好)")

if freq_diff > 50 and logp_diff > 20:
    print(f"\n  [结论] Q1 靶点富集到的通路在参数空间中具有显著更高的稳健性和")
    print(f"         显著性，说明靶点-通路的映射是高质量、有方向性的。")
    print(f"         即使 Q1 计数率为 {simple_q1_rate:.0f}%，富集通路的整体质量")
    print(f"         显著优于背景通路池。")
else:
    print(f"\n  [结论] Q1 靶点富集通路与全池质量差异有限，需进一步分析。")

# ── Save ───────────────────────────────────────────────────────────────
output = {
    "total_pathways": len(all_pathways),
    "enriched_pathways": len(enriched_df),
    "metric_comparison": comparison_results,
    "weighted_q1_rates": {
        "simple": simple_q1_rate,
        "frequency_weighted": freq_weighted_q1_rate,
        "logp_weighted": logp_weighted_q1_rate,
        "combined_weighted": combined_weighted_q1_rate,
    },
    "enriched_pathway_ids": enriched_pathways,
}

with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

print(f"\n[OK] 评估结果已保存: {OUTPUT_JSON}")
