#!/usr/bin/env python3
"""
通路稳健性象限分类增强脚本（P0）

直接读取现有的 pathway_rank_stability.csv，
追加二维象限分类（频率 × 显著性强度），输出增强版结果。

用法:
    cd code && python scripts/analysis/augment_pathway_robustness.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd
import json

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = BASE_DIR / "src" / "outputs" / "sensitivity"
INPUT_CSV = OUTPUTS_DIR / "pathway_rank_stability.csv"
OUTPUT_CSV = OUTPUTS_DIR / "pathway_rank_stability_augmented.csv"
OUTPUT_JSON = OUTPUTS_DIR / "pathway_robustness_quadrant.json"

# ── Load ───────────────────────────────────────────────────────────────
print("=" * 70)
print("通路稳健性象限分类增强")
print("=" * 70)

df = pd.read_csv(INPUT_CSV)
print(f"加载 {len(df)} 条通路记录")

# ── Data-driven median splits ──────────────────────────────────────────
freq_median = df["frequency"].median()
logp_median = df["avg_logp"].median()

print(f"\n数据驱动中位数:")
print(f"  频率中位数: {freq_median:.4f}")
print(f"  显著性中位数 (-log10P): {logp_median:.4f}")

# ── Quadrant classification ────────────────────────────────────────────
def classify_pathway_quadrant(row):
    freq_high = row["frequency"] >= freq_median
    sig_high = row["avg_logp"] >= logp_median
    if freq_high and sig_high:
        return "Q1: 候选稳健通路"
    elif not freq_high and sig_high:
        return "Q2: 低频高显著"
    elif freq_high and not sig_high:
        return "Q3: 稳健但效应弱"
    else:
        return "Q4: 可忽略"

df["robustness_quadrant"] = df.apply(classify_pathway_quadrant, axis=1)

# ── Summary by quadrant ────────────────────────────────────────────────
print(f"\n{'象限':<20} {'数量':<6} {'占比':<8} {'平均频率':<10} {'平均-log10P':<12}")
print("-" * 60)

quadrant_summary = {}
for q in ["Q1: 候选稳健通路", "Q2: 低频高显著", "Q3: 稳健但效应弱", "Q4: 可忽略"]:
    sub = df[df["robustness_quadrant"] == q]
    if len(sub) > 0:
        pct = len(sub) / len(df) * 100
        avg_freq = sub["frequency"].mean()
        avg_logp = sub["avg_logp"].mean()
        print(f"{q:<20} {len(sub):<6} {pct:<7.1f}% {avg_freq:<10.2f} {avg_logp:<12.2f}")
        quadrant_summary[q] = {
            "count": int(len(sub)),
            "percentage": round(pct, 2),
            "avg_frequency": round(float(avg_freq), 4),
            "avg_logp": round(float(avg_logp), 4),
        }

# ── Q1 candidate robust pathways detail ────────────────────────────────
q1_df = df[df["robustness_quadrant"] == "Q1: 候选稳健通路"].sort_values("avg_rank")

print(f"\nQ1 候选稳健通路详情 (n={len(q1_df)}):")
print(f"{'通路名称':<55} {'来源':<8} {'类别':<18} {'频率':<8} {'排名':<6} {'-log10P':<8}")
print("-" * 105)

CATEGORY_LABELS = {
    "direct_ox": "直接氧化应激",
    "cross_talk": "交叉调控",
    "signaling_backbone": "信号主干",
    "other": "其他",
}

for _, r in q1_df.iterrows():
    cat = CATEGORY_LABELS.get(r["pathway_category"], r["pathway_category"])
    name = r["pathway"][:52] + "..." if len(r["pathway"]) > 55 else r["pathway"]
    print(f"{name:<55} {r['source']:<8} {cat:<18} {r['frequency']:<7.1%} {r['avg_rank']:<6.1f} {r['avg_logp']:<8.2f}")

# ── Save augmented CSV ─────────────────────────────────────────────────
df.to_csv(OUTPUT_CSV, index=False)
print(f"\n[OK] 增强版CSV已保存: {OUTPUT_CSV}")

# ── Save JSON summary ──────────────────────────────────────────────────
output_json = {
    "method": "二维象限分类（频率 × 显著性强度）",
    "thresholds": {
        "frequency_median": round(float(freq_median), 4),
        "avg_logp_median": round(float(logp_median), 4),
    },
    "total_pathways": int(len(df)),
    "quadrant_summary": quadrant_summary,
    "q1_candidate_robust_pathways": [
        {
            "pathway": r["pathway"],
            "source": r["source"],
            "pathway_category": r["pathway_category"],
            "frequency": round(float(r["frequency"]), 4),
            "avg_rank": round(float(r["avg_rank"]), 2),
            "avg_logp": round(float(r["avg_logp"]), 4),
            "avg_ef": round(float(r["avg_ef"]), 4),
        }
        for _, r in q1_df.iterrows()
    ],
}

with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
    json.dump(output_json, f, indent=2, ensure_ascii=False)
print(f"[OK] 象限分类JSON已保存: {OUTPUT_JSON}")
