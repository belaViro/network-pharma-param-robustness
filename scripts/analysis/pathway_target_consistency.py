#!/usr/bin/env python3
"""
靶点-功能富集条目一致性验证脚本（P0）

双重验证策略：
1. 视角A（靶点驱动）: Q1候选稳健靶点 → 富集分析 → 具体功能条目列表
2. 视角B（独立扫描）: 参数扫描中高频高显著功能条目
3. 一致性: 两个视角的交集比例

用法:
    cd code && python scripts/analysis/pathway_target_consistency.py
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
PATHWAY_CSV = OUTPUTS_DIR / "pathway_rank_stability.csv"
OUTPUT_JSON = OUTPUTS_DIR / "pathway_target_consistency.json"

PVALUE_CUTOFF = 0.01

# ── Load Q1 targets ────────────────────────────────────────────────────
print("=" * 70)
print("靶点-功能富集条目一致性验证（双重验证策略）")
print("=" * 70)

target_df = pd.read_csv(TARGET_CSV)
q1_genes = sorted(target_df[target_df["quadrant"] == "Q1: 候选稳健"]["gene"].tolist())
print(f"\nQ1 候选稳健靶点: {len(q1_genes)} 个")
print(f"  {', '.join(q1_genes)}")

# ── Perspective A: Target-driven validation ────────────────────────────
print(f"\n{'='*70}")
print("[视角A] 靶点驱动验证：Q1靶点富集分析")
print("=" * 70)

enrichment = cached_enrichment(q1_genes, PVALUE_CUTOFF, BACKGROUND_GENES)
metrics = analyze_enrichment(enrichment, include_sources=["KEGG", "GO_BP"])

print(f"  显著功能条目总数: {metrics['total_significant']}")
print(f"  信号主干:     {metrics['backbone']}")
print(f"  交叉调控:     {metrics['cross_talk']}")
print(f"  直接OX:       {metrics['direct_ox']}")
print(f"  其他:         {metrics['other']}")
print(f"  平均EF:       {metrics['avg_ef']:.2f}")

# 提取相关功能条目（cross_talk + signaling_backbone）
a_driver_pathways = set(
    p["pathway"] for p in metrics["pathway_details"]
    if p["pathway_category"] in ("cross_talk", "signaling_backbone")
)
print(f"\n  靶点驱动的相关功能条目 (cross_talk + signaling_backbone): {len(a_driver_pathways)} 条")

# ── Perspective B: Independent scan validation ─────────────────────────
print(f"\n{'='*70}")
print("[视角B] 独立扫描验证：参数空间高频高显著功能条目")
print("=" * 70)

pathway_df = pd.read_csv(PATHWAY_CSV)

# 高频高显著功能条目：frequency >= 0.5 且 avg_rank <= 15，且分类为相关
scan_robust = pathway_df[
    (pathway_df["frequency"] >= 0.5) &
    (pathway_df["avg_rank"] <= 15) &
    (pathway_df["pathway_category"].isin(["cross_talk", "signaling_backbone"]))
].sort_values("avg_rank")

b_scan_pathways = set(scan_robust["pathway"])
print(f"  参数扫描高频稳健功能条目 (freq>=50%, rank<=15): {len(b_scan_pathways)} 条")

# 也统计一个更严格的版本（freq>=0.6）
scan_strict = pathway_df[
    (pathway_df["frequency"] >= 0.6) &
    (pathway_df["avg_rank"] <= 15) &
    (pathway_df["pathway_category"].isin(["cross_talk", "signaling_backbone"]))
]
b_scan_strict = set(scan_strict["pathway"])
print(f"  严格版 (freq>=60%): {len(b_scan_strict)} 条")

# ── Consistency analysis ───────────────────────────────────────────────
print(f"\n{'='*70}")
print("一致性分析")
print("=" * 70)

overlap = a_driver_pathways & b_scan_pathways
only_in_driver = a_driver_pathways - b_scan_pathways
only_in_scan = b_scan_pathways - a_driver_pathways

consistency_rate = len(overlap) / len(a_driver_pathways) if a_driver_pathways else 0

print(f"\n  靶点驱动功能条目总数: {len(a_driver_pathways)}")
print(f"  扫描高频功能条目总数: {len(b_scan_pathways)}")
print(f"  交集功能条目数:       {len(overlap)}")
print(f"  仅靶点驱动:           {len(only_in_driver)}")
print(f"  仅扫描高频:           {len(only_in_scan)}")
print(f"\n  >>> 双重验证一致性: {consistency_rate:.1%} <<<")

# 严格版一致性
overlap_strict = a_driver_pathways & b_scan_strict
consistency_strict = len(overlap_strict) / len(a_driver_pathways) if a_driver_pathways else 0
print(f"  >>> 严格版一致性 (freq>=60%): {consistency_strict:.1%} <<<")

# ── Detailed overlap list ──────────────────────────────────────────────
if overlap:
    print(f"\n  [PASS] 双重验证通过 - 既是Q1靶点驱动、又在参数扫描中高频出现:")
    for pw in sorted(overlap):
        info = pathway_df[pathway_df["pathway"] == pw].iloc[0]
        print(f"    - {pw[:55]:<55} freq={info['frequency']:<6.1%} rank={info['avg_rank']:<5.1f} logP={info['avg_logp']:<6.2f}")

if only_in_driver:
    print(f"\n  [WARN] 靶点特有 - Q1靶点驱动但扫描中未达高频阈值 ({len(only_in_driver)}条):")
    for pw in sorted(only_in_driver)[:10]:
        print(f"    - {pw}")
    if len(only_in_driver) > 10:
        print(f"    ... 还有 {len(only_in_driver)-10} 条")

if only_in_scan:
    print(f"\n  [INFO] 扫描特有 - 参数扫描高频但非Q1靶点直接驱动 ({len(only_in_scan)}条):")
    for pw in sorted(only_in_scan)[:10]:
        info = pathway_df[pathway_df["pathway"] == pw].iloc[0]
        print(f"    - {pw[:55]:<55} freq={info['frequency']:<6.1%}")
    if len(only_in_scan) > 10:
        print(f"    ... 还有 {len(only_in_scan)-10} 条")

# ── Save JSON ──────────────────────────────────────────────────────────
output = {
    "method": "双重验证策略（靶点驱动 + 独立扫描）",
    "q1_targets_count": len(q1_genes),
    "q1_targets": q1_genes,
    "perspective_a_target_driven": {
        "total_significant_pathways": metrics["total_significant"],
        "relevant_pathways_count": len(a_driver_pathways),
        "relevant_pathways": sorted(a_driver_pathways),
    },
    "perspective_b_independent_scan": {
        "threshold": "frequency >= 0.5 AND avg_rank <= 15",
        "robust_pathways_count": len(b_scan_pathways),
        "robust_pathways": sorted(b_scan_pathways),
        "strict_threshold": "frequency >= 0.6 AND avg_rank <= 15",
        "strict_pathways_count": len(b_scan_strict),
    },
    "consistency": {
        "overlap_count": len(overlap),
        "overlap_pathways": sorted(overlap),
        "only_in_driver_count": len(only_in_driver),
        "only_in_scan_count": len(only_in_scan),
        "consistency_rate": round(consistency_rate, 4),
        "consistency_strict_rate": round(consistency_strict, 4),
    },
    "interpretation": {
        "high_consistency": consistency_rate >= 0.7,
        "narrative": (
            "双重验证高度一致，功能条目稳健性是靶点稳健性的可靠推论。"
            if consistency_rate >= 0.7 else
            "一致性中等，部分靶点驱动功能条目在参数扫描中波动较大，需单独讨论。"
        ),
    },
}

with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"\n[OK] 一致性验证结果已保存: {OUTPUT_JSON}")
