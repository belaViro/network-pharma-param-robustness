#!/usr/bin/env python3
"""
Q1靶点富集结果审计：看看Q1靶点驱动的 GO BP/KEGG 功能富集条目落在哪些稳健性象限

用法:
    cd code && python scripts/analysis/q1_targets_enrichment_audit.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd
import json

from src.config import BACKGROUND_GENES
from src.sensitivity.enrichment_cache import cached_enrichment
from src.sensitivity.analysis import analyze_enrichment

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
TARGET_CSV = OUTPUTS_DIR / "target_robustness.csv"
PATHWAY_CSV = OUTPUTS_DIR / "pathway_rank_stability_augmented.csv"
OUTPUT_JSON = OUTPUTS_DIR / "q1_enrichment_audit.json"

PVALUE_CUTOFF = 0.01
Q1_TERM_LABEL = "Q1: 候选稳健功能条目"
Q1_TERM_OLD_LABEL = "Q1: 候选稳健通路"

# ── Load Q1 targets ────────────────────────────────────────────────────
print("=" * 70)
print("Q1靶点富集结果审计：GO BP/KEGG 功能富集条目落点分析")
print("=" * 70)

target_df = pd.read_csv(TARGET_CSV)
q1_genes = sorted(target_df[target_df["quadrant"] == "Q1: 候选稳健"]["gene"].tolist())
print(f"\nQ1 候选稳健靶点: {len(q1_genes)} 个")

# ── Run enrichment on Q1 targets ───────────────────────────────────────
print(f"\n运行 Q1 靶点富集分析 (Pvalue={PVALUE_CUTOFF})...")
enrichment = cached_enrichment(q1_genes, PVALUE_CUTOFF, BACKGROUND_GENES)
metrics = analyze_enrichment(enrichment, include_sources=["KEGG", "GO_BP"])

print(f"  显著功能条目总数: {metrics['total_significant']}")
print(f"  直接OX:       {metrics['direct_ox']}")
print(f"  交叉调控:     {metrics['cross_talk']}")
print(f"  信号主干:     {metrics['backbone']}")
print(f"  其他:         {metrics['other']}")

# ── Load augmented pathway table ───────────────────────────────────────
pathway_df = pd.read_csv(PATHWAY_CSV)
pathway_df["robustness_quadrant"] = pathway_df["robustness_quadrant"].replace({Q1_TERM_OLD_LABEL: Q1_TERM_LABEL})
# 建立功能条目名称→象限的映射
pathway_quadrant_map = dict(zip(pathway_df["pathway"], pathway_df["robustness_quadrant"]))
pathway_freq_map = dict(zip(pathway_df["pathway"], pathway_df["frequency"]))
pathway_rank_map = dict(zip(pathway_df["pathway"], pathway_df["avg_rank"]))

# ── Audit each enriched pathway ────────────────────────────────────────
print(f"\n{'='*70}")
print("Q1靶点富集功能条目落点分析")
print("=" * 70)

results_by_source = {}
quadrant_summary = {Q1_TERM_LABEL: 0, "Q2: 低频高显著": 0,
                    "Q3: 稳健但效应弱": 0, "Q4: 可忽略": 0, "未在扫描中出现": 0}
category_quadrant_breakdown = {}

for source, df_sig in enrichment["significant"].items():
    if not isinstance(df_sig, pd.DataFrame) or len(df_sig) == 0:
        continue
    
    results_by_source[source] = []
    print(f"\n--- {source} ({len(df_sig)} 条) ---")
    print(f"{'功能条目名称':<55} {'五层分类':<12} {'象限':<18} {'频率':<8} {'排名':<6}")
    print("-" * 105)
    
    for _, row in df_sig.iterrows():
        term = row.get("Term", "")
        category = "N/A"
        
        # 匹配五层分类
        for p in metrics["pathway_details"]:
            if p["pathway"] == term:
                category = p["pathway_category"]
                break
        
        # 匹配稳健性象限
        quadrant = pathway_quadrant_map.get(term, "未在扫描中出现")
        freq = pathway_freq_map.get(term, 0.0)
        rank = pathway_rank_map.get(term, 999)
        
        quadrant_summary[quadrant] = quadrant_summary.get(quadrant, 0) + 1
        
        cat_key = f"{category} → {quadrant}"
        category_quadrant_breakdown[cat_key] = category_quadrant_breakdown.get(cat_key, 0) + 1
        
        q_short = quadrant.replace(Q1_TERM_LABEL, "Q1").replace(Q1_TERM_OLD_LABEL, "Q1").replace("Q2: 低频高显著", "Q2")\
                         .replace("Q3: 稳健但效应弱", "Q3").replace("Q4: 可忽略", "Q4")\
                         .replace("未在扫描中出现", "NEW")
        
        name = term[:52] + "..." if len(term) > 55 else term
        print(f"{name:<55} {category:<12} {q_short:<18} {freq:<7.1%} {rank:<6.1f}")
        
        results_by_source[source].append({
            "pathway": term,
            "category": category,
            "quadrant": quadrant,
            "frequency": float(freq),
            "avg_rank": float(rank),
            "pvalue": float(row.get("Adjusted P-value", 1.0)),
            "ef": float(row.get("EF", 0)),
        })

# ── Summary by source ──────────────────────────────────────────────────
print(f"\n{'='*70}")
print("按来源分层统计")
print("=" * 70)

total_enriched = metrics['total_significant']

# 按来源统计
source_stats = {}
for source, pathways in results_by_source.items():
    total = len(pathways)
    q1 = sum(1 for p in pathways if p["quadrant"] == Q1_TERM_LABEL)
    q3 = sum(1 for p in pathways if p["quadrant"] == "Q3: 稳健但效应弱")
    q4 = sum(1 for p in pathways if p["quadrant"] == "Q4: 可忽略")
    q2 = sum(1 for p in pathways if p["quadrant"] == "Q2: 低频高显著")
    q1_rate = q1 / total * 100 if total else 0
    source_stats[source] = {"total": total, "q1": q1, "q2": q2, "q3": q3, "q4": q4, "q1_rate": q1_rate}
    print(f"  {source:12s}: 总数={total:2d} | Q1={q1:2d}({q1_rate:5.1f}%) | Q2={q2} | Q3={q3:2d} | Q4={q4}")

# 核心功能条目（KEGG + GO_BP）
kegg_bp_pathways = []
for source in ("KEGG", "GO_BP"):
    if source in results_by_source:
        kegg_bp_pathways.extend(results_by_source[source])
kegg_bp_total = len(kegg_bp_pathways)
kegg_bp_q1 = sum(1 for p in kegg_bp_pathways if p["quadrant"] == Q1_TERM_LABEL)
kegg_bp_rate = kegg_bp_q1 / kegg_bp_total * 100 if kegg_bp_total else 0

print(f"\n  {'='*50}")
print(f"  KEGG + GO_BP (系统级功能条目): 总数={kegg_bp_total} | Q1={kegg_bp_q1}({kegg_bp_rate:.1f}%)")
print(f"  {'='*50}")

print(f"\n按稳健性象限分布（全部）:")
for q, cnt in quadrant_summary.items():
    if cnt > 0:
        pct = cnt / total_enriched * 100
        print(f"  {q:<20}: {cnt:3d} 条 ({pct:5.1f}%)")

print(f"\n按五层分类 × 象限交叉:")
for cat_quad, cnt in sorted(category_quadrant_breakdown.items(), key=lambda x: -x[1]):
    print(f"  {cat_quad}: {cnt} 条")

# ── Key insight ────────────────────────────────────────────────────────
q1_count = quadrant_summary.get(Q1_TERM_LABEL, 0)
q1_rate = q1_count / total_enriched * 100 if total_enriched else 0

print(f"\n{'='*70}")
print("核心结论")
print("=" * 70)
print(f"\n  [全部功能条目] Q1靶点富集功能条目中，落在 Q1 稳健象限的比例: {q1_rate:.1f}%")
print(f"  [系统功能条目] KEGG + GO_BP 中，落在 Q1 稳健象限的比例: {kegg_bp_rate:.1f}%")
print(f"  未在参数扫描中出现的'新'功能条目: {quadrant_summary.get('未在扫描中出现', 0)} 条")

if kegg_bp_rate >= 90:
    print(f"\n  [结论] 系统级功能条目（KEGG/GO_BP）的稳健传导率极高（{kegg_bp_rate:.0f}%），")
    print(f"         靶点→功能条目的传导在信号网络层面几乎完美。")
    print(f"         GO_MF/CC 比例较低，因其本质是分子注释而非系统功能条目。")
elif kegg_bp_rate >= 70:
    print(f"\n  [结论] 系统级功能条目稳健传导率良好，靶点→功能条目的传导有效。")
else:
    print(f"\n  [结论] 系统级功能条目稳健传导率不足，需讨论参数敏感性。")

# ── Save JSON ──────────────────────────────────────────────────────────
output = {
    "q1_targets": q1_genes,
    "q1_targets_count": len(q1_genes),
    "enrichment_summary": {
        "total_significant": metrics["total_significant"],
        "direct_ox": metrics["direct_ox"],
        "cross_talk": metrics["cross_talk"],
        "backbone": metrics["backbone"],
        "other": metrics["other"],
    },
    "quadrant_distribution": quadrant_summary,
    "category_quadrant_breakdown": category_quadrant_breakdown,
    "key_metrics": {
        "overall_q1_rate": round(q1_rate, 2),
        "system_functional_term_q1_rate": round(kegg_bp_rate, 2),
    },
    "pathway_details": results_by_source,
}

with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"\n[OK] 审计结果已保存: {OUTPUT_JSON}")
