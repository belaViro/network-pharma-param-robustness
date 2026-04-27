#!/usr/bin/env python3
"""
PubMed 文献共现验证脚本

用 PubMed 文献数据作为外部"金标准"，验证当前 KEGG 通路分类的客观性。

核心逻辑：
1. 对每条 KEGG 通路名称，在 PubMed 中搜索该通路与 "oxidative stress" 的共现文章
2. 文章数 >= threshold 的通路，被 PubMed 金标准判定为"氧化应激相关"
3. 与当前框架的 cross_talk + signaling_backbone 做一致性对比
4. 计算 Jaccard、Precision、Recall

API: NCBI E-utilities (免费，限制 ~3 req/sec)
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Set

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd

from src.sensitivity.classification import classify_pathway
from src.connectors.pubmed_connector import PubMedConnector, PubMedConfig

# 共现文章数阈值
PMID_THRESHOLD = 10


def load_kegg_pathways(stability_file: Path) -> tuple[List[str], Dict[str, str]]:
    """Load KEGG pathway list and their framework classifications."""
    df = pd.read_csv(stability_file)
    kegg_df = df[df["source"] == "KEGG"]
    pathways = kegg_df["pathway"].unique().tolist()
    classifications = {
        row["pathway"]: row["pathway_category"]
        for _, row in kegg_df.iterrows()
    }
    return pathways, classifications


def main():
    print("=" * 70)
    print("PubMed 文献共现验证：KEGG 通路氧化应激相关性外部金标准")
    print("=" * 70)

    # =====================================================================
    # Step 1: Load KEGG pathways
    # =====================================================================
    print("\n[Step 1] 加载当前框架的 KEGG 通路...")
    stability_file = (
        Path(__file__).resolve().parent.parent.parent
        / "src" / "outputs" / "sensitivity" / "pathway_rank_stability.csv"
    )

    if stability_file.exists():
        kegg_pathways, current_classifications = load_kegg_pathways(stability_file)
        print(f"  从 stability 文件加载: {len(kegg_pathways)} 条 KEGG 通路")
    else:
        # Fallback test set
        kegg_pathways = [
            "MAPK signaling pathway",
            "PI3K-Akt signaling pathway",
            "Pathways in cancer",
            "TNF signaling pathway",
            "AGE-RAGE signaling pathway in diabetic complications",
            "Nitrogen metabolism",
            "Glutathione metabolism",
        ]
        current_classifications = {
            p: classify_pathway(p, "KEGG") for p in kegg_pathways
        }
        print(f"  使用内置测试通路: {len(kegg_pathways)} 条")

    if not kegg_pathways:
        print("  [ERROR] 无 KEGG 通路可验证")
        return

    # =====================================================================
    # Step 2: PubMed co-occurrence queries
    # =====================================================================
    print(f"\n[Step 2] PubMed 共现查询（阈值 >={PMID_THRESHOLD} 篇文章）...")
    print(f"  预计耗时: ~{len(kegg_pathways) * 0.35:.0f} 秒")

    connector = PubMedConnector(config=PubMedConfig(delay_seconds=0.35))
    pubmed_counts = connector.batch_cooccurrence(
        kegg_pathways,
        context_term="oxidative stress",
        use_cache=True,
        show_progress=True,
    )

    pubmed_gold_standard: Set[str] = {
        p for p, c in pubmed_counts.items() if c >= PMID_THRESHOLD
    }
    print(f"\n  PubMed 金标准判定为相关: {len(pubmed_gold_standard)}/{len(kegg_pathways)}")

    # =====================================================================
    # Step 3: Consistency metrics
    # =====================================================================
    print("\n[Step 3] 计算 PubMed 金标准 vs 当前框架的一致性...")

    predicted_positive = {
        k for k, v in current_classifications.items()
        if v in ("cross_talk", "signaling_backbone")
    }
    all_pathways = set(kegg_pathways)

    metrics = PubMedConnector.compute_metrics(
        pubmed_gold_standard, predicted_positive, all_pathways
    )

    print(f"\n  {'Metric':<20} {'Value':<10}")
    print(f"  {'-'*30}")
    print(f"  {'True Positives':<20} {metrics['tp']}")
    print(f"  {'False Positives':<20} {metrics['fp']}")
    print(f"  {'False Negatives':<20} {metrics['fn']}")
    print(f"  {'True Negatives':<20} {metrics['tn']}")
    print(f"  {'Precision':<20} {metrics['precision']:.2%}")
    print(f"  {'Recall':<20} {metrics['recall']:.2%}")
    print(f"  {'F1 Score':<20} {metrics['f1']:.2%}")
    print(f"  {'Jaccard Overlap':<20} {metrics['jaccard']:.2%}")
    print(f"  {'Accuracy':<20} {metrics['accuracy']:.2%}")

    # =====================================================================
    # Step 4: Detailed comparison table
    # =====================================================================
    print("\n[Step 4] 详细通路对比:")
    print(f"\n  {'通路名称':<50} {'PubMed n':<10} {'当前分类':<15} {'PubMed':<8} {'框架':<8}")
    print(f"  {'-'*90}")

    for pathway in sorted(kegg_pathways):
        count = pubmed_counts.get(pathway, -1)
        current_cat = current_classifications.get(pathway, "N/A")
        pubmed_label = "相关" if pathway in pubmed_gold_standard else "无关"
        framework_label = "相关" if current_cat in ("cross_talk", "signaling_backbone") else "无关"
        agree = "[OK]" if pubmed_label == framework_label else "[XX]"
        print(f"  {pathway[:48]:<50} {count:>8}   {current_cat:<15} {pubmed_label:<8} {framework_label:<8} {agree}")

    # =====================================================================
    # Step 5: Missed high-count pathways
    # =====================================================================
    missed = sorted(
        [(p, pubmed_counts[p]) for p in pubmed_gold_standard
         if current_classifications.get(p) == "other"],
        key=lambda x: x[1], reverse=True
    )

    if missed:
        print(f"\n[Step 5] 高共现但被当前框架归为 'other' 的通路（潜在漏网，共 {len(missed)} 条）:")
        for p, count in missed[:20]:
            print(f"  n={count:4d} | {p}")
    else:
        print("\n[Step 5] 无漏网通路（PubMed 相关的通路都被框架捕获）")

    # =====================================================================
    # Step 6: False alarms
    # =====================================================================
    false_alarms = sorted(
        [(p, pubmed_counts[p]) for p in predicted_positive
         if p not in pubmed_gold_standard and pubmed_counts.get(p, -1) >= 0],
        key=lambda x: x[1]
    )

    if false_alarms:
        print(f"\n[Step 6] 低共现但被当前框架归为相关的通路（潜在误报，共 {len(false_alarms)} 条）:")
        for p, count in false_alarms[:20]:
            print(f"  n={count:4d} | {p}")
    else:
        print("\n[Step 6] 无误报通路")

    # =====================================================================
    # Step 7: Save results
    # =====================================================================
    output = {
        "pmid_threshold": PMID_THRESHOLD,
        "ncbi_delay": connector.config.delay_seconds,
        "total_kegg_pathways": len(kegg_pathways),
        "pubmed_positive": len(pubmed_gold_standard),
        "framework_positive": len(predicted_positive),
        "metrics": {
            k: round(v, 4) if isinstance(v, float) else v
            for k, v in metrics.items()
        },
        "pathway_details": [
            {
                "pathway": p,
                "pubmed_count": pubmed_counts.get(p, -1),
                "pubmed_positive": p in pubmed_gold_standard,
                "current_category": current_classifications.get(p, "N/A"),
                "framework_positive": p in predicted_positive,
            }
            for p in sorted(kegg_pathways)
        ],
        "missed_high_count": [
            {"pathway": p, "pubmed_count": c}
            for p, c in missed[:20]
        ],
        "false_alarms_low_count": [
            {"pathway": p, "pubmed_count": c}
            for p, c in false_alarms[:20]
        ],
    }

    out_path = (
        Path(__file__).resolve().parent.parent.parent
        / "src" / "outputs" / "sensitivity" / "pubmed_validation.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    print(f"\n[OK] 验证结果已保存: {out_path}")


if __name__ == "__main__":
    main()
