#!/usr/bin/env python3
"""
方案 B 概念验证：基因集重叠法

用 KEGG 通路基因集与氧化应激基因集的重叠度来判定通路关联性。
"""

import re
import time
from pathlib import Path

import pandas as pd
import requests

# KEGG API endpoints
KEGG_GET_URL = "https://rest.kegg.jp/get"
KEGG_LINK_URL = "https://rest.kegg.jp/link"

# Load oxidative stress gene set
OX_GENES_PATH = Path("src/outputs/shared/disease_genes/merged_disease_genes.csv")
ox_df = pd.read_csv(OX_GENES_PATH)
ox_gene_set = set(g.strip().upper() for g in ox_df["Gene_Symbol"].dropna().unique())
print(f"Loaded oxidative stress gene set: {len(ox_gene_set)} genes")


def get_kegg_pathway_genes(kegg_id: str) -> set:
    """Fetch gene symbols for a KEGG pathway via KEGG API."""
    url = f"{KEGG_GET_URL}/{kegg_id}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        print(f"  [ERROR] Failed to fetch {kegg_id}: {e}")
        return set()

    genes = set()
    in_gene_section = False
    for line in resp.text.split("\n"):
        if line.startswith("GENE"):
            in_gene_section = True
            # Parse first gene on same line
            m = re.search(r"GENE\s+\d+\s+([A-Za-z0-9_-]+)", line)
            if m:
                genes.add(m.group(1).strip().upper())
        elif in_gene_section:
            if line.startswith(" ") or line.startswith("\t"):
                # Continuation of GENE section
                m = re.search(r"\s+\d+\s+([A-Za-z0-9_-]+)", line)
                if m:
                    genes.add(m.group(1).strip().upper())
            else:
                in_gene_section = False
    return genes


def compute_overlap(pathway_genes: set, ox_genes: set) -> dict:
    """Compute Jaccard and overlap statistics."""
    intersection = pathway_genes & ox_genes
    union = pathway_genes | ox_genes
    jaccard = len(intersection) / len(union) if union else 0
    overlap_ratio = len(intersection) / len(pathway_genes) if pathway_genes else 0
    return {
        "pathway_genes": len(pathway_genes),
        "ox_genes": len(intersection),
        "jaccard": jaccard,
        "overlap_ratio": overlap_ratio,
    }


# Test pathways
test_pathways = [
    ("hsa04010", "MAPK signaling pathway", "signaling_backbone"),
    ("hsa04151", "PI3K-Akt signaling pathway", "signaling_backbone"),
    ("hsa04668", "TNF signaling pathway", "cross_talk"),
    ("hsa04216", "Ferroptosis", "cross_talk"),
    ("hsa00910", "Nitrogen metabolism", "other"),
    ("hsa04962", "AGE-RAGE signaling pathway in diabetic complications", "cross_talk"),
    ("hsa04974", "Protein digestion and absorption", "other"),  # likely other
    ("hsa04670", "Leukocyte transendothelial migration", "cross_talk"),
]

print("\n" + "=" * 80)
print("方案 B 概念验证：KEGG 通路基因集 vs 氧化应激基因集重叠分析")
print("=" * 80)
print(f"\n{'通路名称':<50} {'框架分类':<18} {'通路基因数':>10} {'OX基因数':>10} {'重叠率':>8} {'Jaccard':>8}")
print("-" * 110)

results = []
for kegg_id, name, fw_cat in test_pathways:
    genes = get_kegg_pathway_genes(kegg_id)
    stats = compute_overlap(genes, ox_gene_set)
    results.append({
        "kegg_id": kegg_id,
        "name": name,
        "fw_category": fw_cat,
        **stats,
    })
    print(f"{name[:48]:<50} {fw_cat:<18} {stats['pathway_genes']:>10} {stats['ox_genes']:>10} {stats['overlap_ratio']:>8.2%} {stats['jaccard']:>8.4f}")
    time.sleep(0.5)  # Be nice to KEGG API

print("\n" + "=" * 80)
print("解读")
print("=" * 80)

# Sort by overlap ratio
results_sorted = sorted(results, key=lambda x: x["overlap_ratio"], reverse=True)
print("\n按重叠率排序：")
for r in results_sorted:
    marker = "✓" if (r["fw_category"] in ("cross_talk", "signaling_backbone") and r["overlap_ratio"] > 0.05) or (r["fw_category"] == "other" and r["overlap_ratio"] < 0.05) else "?"
    print(f"  {marker} {r['name'][:45]:<45} | 重叠率={r['overlap_ratio']:>6.2%} | 框架={r['fw_category']}")
