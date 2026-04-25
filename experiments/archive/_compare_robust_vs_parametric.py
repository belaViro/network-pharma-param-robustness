"""
稳健靶点 vs 各参数组合的富集效果对比
比较Q1核心稳健靶点(33个)和Tier 1绝对稳健靶点(9个)的富集效果
是否优于每个参数组合的富集效果
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import re
import json
import hashlib
from dataclasses import dataclass
from typing import List, Dict

from src.config import BACKGROUND_GENES
from src.enrichment.enrichment_analysis import run_enrichment_analysis

# ============================================================================
# 1. 加载分类器（复用 MVP v2 代码）
# ============================================================================
GO_OX_FILE = Path("src/data/oxidative_stress_go_terms.txt")
GO_OX_IDS = set()
if GO_OX_FILE.exists():
    with open(GO_OX_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            go_id = line.strip().split('\t')[0]
            GO_OX_IDS.add(go_id)

CROSS_TALK_PATHWAYS = {
    "AGE-RAGE signaling pathway in diabetic complications",
    "IL-17 signaling pathway", "Ferroptosis",
    "Glutathione metabolism", "Peroxisome",
    "Fluid shear stress and atherosclerosis",
    "NOD-like receptor signaling pathway",
    "TNF signaling pathway",
    "Toll-like receptor signaling pathway",
    "NF-kappa B signaling pathway",
}

BACKBONE_KEYWORDS = [
    'MAPK', 'PI3K-Akt', 'Rap1', 'Ras signaling',
    'receptor tyrosine kinase', 'protein tyrosine kinase',
    'cytokine-mediated', 'calcium signaling', 'cAMP signaling',
    'chemokine signaling', 'JAK-STAT',
    'EGFR tyrosine kinase inhibitor resistance',
]

def classify_pathway(term: str) -> str:
    term_lower = term.lower()
    match = re.search(r'\(GO:(\d+)\)', term)
    if match:
        go_id = f"GO:{match.group(1)}"
        if go_id in GO_OX_IDS:
            return "direct_ox"
    if term in CROSS_TALK_PATHWAYS:
        return "cross_talk"
    for kw in BACKBONE_KEYWORDS:
        if re.search(r'\b' + re.escape(kw.lower()) + r'\b', term_lower):
            return "signaling_backbone"
    return "other"

# ============================================================================
# 2. 缓存（复用已有缓存）
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

def analyze_enrichment(enrichment):
    """从富集结果中提取分类指标"""
    category_counts = {"direct_ox": 0, "cross_talk": 0, "signaling_backbone": 0, "other": 0}
    total_significant = 0
    avg_ef = 0.0
    min_pval = 1.0
    pathway_details = []

    if enrichment["status"] == "成功" and "significant" in enrichment:
        for source, df_sig in enrichment["significant"].items():
            if isinstance(df_sig, pd.DataFrame) and len(df_sig) > 0:
                total_significant += len(df_sig)
                for _, row_sig in df_sig.iterrows():
                    term_original = str(row_sig.get("Term", ""))
                    pathway_category = classify_pathway(term_original)
                    category_counts[pathway_category] += 1
                    ef = row_sig.get("EF", 0)
                    pval = row_sig.get("Adjusted P-value", 1.0)
                    pathway_details.append({
                        "pathway": term_original,
                        "source": source,
                        "pathway_category": pathway_category,
                        "pvalue": pval,
                        "ef": ef,
                    })
                    if isinstance(ef, (int, float)) and isinstance(pval, (int, float)):
                        avg_ef += ef
                        if pval < min_pval:
                            min_pval = pval
        if total_significant > 0:
            avg_ef /= total_significant

    return {
        "total_significant": total_significant,
        "direct_ox": category_counts["direct_ox"],
        "cross_talk": category_counts["cross_talk"],
        "backbone": category_counts["signaling_backbone"],
        "other": category_counts["other"],
        "avg_ef": round(avg_ef, 2),
        "min_pvalue": min_pval,
        "pathway_details": pathway_details,
    }

# ============================================================================
# 3. 加载各参数组合的富集结果
# ============================================================================
print("=" * 80)
print("稳健靶点 vs 各参数组合 富集效果对比")
print("=" * 80)

# 加载MVP v2组合结果
combo_df = pd.read_csv("src/outputs/sensitivity/mvp_v2_combinations.csv")
print(f"\n加载 {len(combo_df)} 个参数组合的结果")

# 各组合指标统计
print(f"\n--- 各参数组合富集指标统计 ---")
for idx, row in combo_df.iterrows():
    print(f"  {row['combo']:<25} "
          f"显著={row['significant_pathways']:3d} "
          f"主干={row['backbone_count']:2d} "
          f"交叉={row['cross_talk_count']:2d} "
          f"直接OX={row['direct_ox_count']:2d} "
          f"平均EF={row['avg_ef']:6.1f} "
          f"minP={row['min_pvalue']:.2e}")

# 汇总统计
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
# 4. 运行Q1稳健靶点的富集分析
# ============================================================================
target_df = pd.read_csv("src/outputs/sensitivity/mvp_v2_target_robustness.csv")
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
    metrics = analyze_enrichment(enrichment)
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
# 5. 对比分析：稳健靶点 vs 各参数组合
# ============================================================================
print(f"\n{'='*80}")
print("对比分析：稳健靶点富集效果是否优于各参数组合？")
print("=" * 80)

# 定义对比指标和方向
comparisons = [
    ("总显著通路数 (越多越好)", "significant_pathways", "higher"),
    ("信号主干道通路数 (越多越好)", "backbone_count", "higher"),
    ("交叉调控通路数 (越多越好)", "cross_talk_count", "higher"),
    ("平均Enrichment Factor (越高越好)", "avg_ef", "higher"),
    ("最小P-value (越小越好)", "min_pvalue", "lower"),
]

for label, metric, direction in comparisons:
    combo_vals = combo_df[metric]
    print(f"\n--- {label} ---")
    print(f"  参数组合范围: {combo_vals.min():.4f} ~ {combo_vals.max():.4f} (均值 {combo_vals.mean():.4f})")

    for rlabel, rmetrics in results.items():
        if metric == "min_pvalue":
            rval = rmetrics.get("min_pvalue", 1.0)
        else:
            rval = rmetrics.get(metric.split("_")[0] if "_" in metric else metric, 0)
            # map metric names
            metric_map = {
                "significant_pathways": "total_significant",
                "backbone_count": "backbone",
                "cross_talk_count": "cross_talk",
                "avg_ef": "avg_ef",
            }
            rval = rmetrics.get(metric_map.get(metric, metric), 0)

        if direction == "higher":
            rank = int((combo_vals < rval).sum())
            better = rval >= combo_vals.max()
        else:
            rank = int((combo_vals > rval).sum())
            better = rval <= combo_vals.min()

        percentile = rank / len(combo_df) * 100
        status = "优于全部组合" if better else "未优于全部组合"
        print(f"  {rlabel}: {rval:.4f} (优于 {rank}/{len(combo_df)} 个组合 = {percentile:.0f}% 分位) [{status}]")


# ============================================================================
# 6. 综合判定
# ============================================================================
print(f"\n{'='*80}")
print("综合判定")
print("=" * 80)

q1_metrics = results["Q1_核心稳健(33个)"]
t1_metrics = results["Tier1_绝对稳健(9个)"]

# 定义判定规则
judgments = []

# Q1判定
q1_total_better = q1_metrics["total_significant"] >= combo_df["significant_pathways"].max()
q1_backbone_better = q1_metrics["backbone"] >= combo_df["backbone_count"].max()
q1_ef_better = q1_metrics["avg_ef"] >= combo_df["avg_ef"].max()

# 对于min_pvalue，越小越好
q1_pval_better = q1_metrics["min_pvalue"] <= combo_df["min_pvalue"].min()

q1_score = sum([q1_total_better, q1_backbone_better, q1_ef_better, q1_pval_better])
def check(val, better, combo_vals, metric_name, higher_better=True):
    if higher_better:
        max_val = combo_vals.max()
        flag = "+" if better else "-"
        return f"  [{flag}] {metric_name}: {val:.2f} vs 最佳组合 {max_val:.2f}"
    else:
        min_val = combo_vals.min()
        flag = "+" if better else "-"
        return f"  [{flag}] {metric_name}: {val:.2e} vs 最佳组合 {min_val:.2e}"

print(f"\nQ1 核心稳健靶点 ({len(q1_genes)}个):")
print(check(q1_metrics['total_significant'], q1_total_better, combo_df['significant_pathways'], '总显著通路', True))
print(check(q1_metrics['backbone'], q1_backbone_better, combo_df['backbone_count'], '主干道数', True))
print(check(q1_metrics['avg_ef'], q1_ef_better, combo_df['avg_ef'], '平均EF', True))
print(check(q1_metrics['min_pvalue'], q1_pval_better, combo_df['min_pvalue'], '最小P-value', False))
print(f"  => 通过 {q1_score}/4 项")

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
        "total_significant": q1_total_better,
        "backbone": q1_backbone_better,
        "avg_ef": q1_ef_better,
        "min_pvalue": q1_pval_better,
    },
    "tier1_beats_all": {
        "total_significant": t1_total_better,
        "backbone": t1_backbone_better,
        "avg_ef": t1_ef_better,
        "min_pvalue": t1_pval_better,
    },
}

import json
with open("src/outputs/sensitivity/robust_vs_param_comparison.json", "w") as f:
    json.dump(comparison_summary, f, indent=2, default=str)
print(f"\n对比结果已保存: src/outputs/sensitivity/robust_vs_param_comparison.json")
