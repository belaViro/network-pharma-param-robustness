#!/usr/bin/env python3
"""
Reactome 外部验证脚本

用 Reactome curated 数据作为"金标准"，验证当前 KEGG 通路分类的可靠性。

流程：
1. Reactome ContentService API 搜索氧化应激相关通路
2. 获取 Reactome 通路的 KEGG cross-reference 映射
3. 建立 KEGG ID -> 名称 映射（通过 KEGG API）
4. 与当前框架的 cross_talk + signaling_backbone 做一致性对比
"""

import html
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

# ---------------------------------------------------------------------------
# API 配置
# ---------------------------------------------------------------------------
REACTOME_BASE = "https://reactome.org/ContentService"
KEGG_BASE = "https://rest.kegg.jp"
TIMEOUT = 60

SEARCH_TERMS = [
    "oxidative stress",
    "reactive oxygen species",
    "glutathione",
    "ferroptosis",
    "Nrf2",
]

# ---------------------------------------------------------------------------
# Reactome API
# ---------------------------------------------------------------------------

def reactome_search(query: str, species: str = "Homo sapiens", page_size: int = 50) -> List[Dict]:
    """在 Reactome 中搜索通路。返回通路条目列表。"""
    url = f"{REACTOME_BASE}/search/query"
    params = {
        "query": query,
        "species": species,
        "types": "Pathway",
        "pageSize": page_size,
    }
    try:
        resp = requests.get(url, params=params, headers={"accept": "application/json"}, timeout=TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        # Reactome 返回结构: {"results": [{"entries": [...]}]}
        entries = []
        for group in data.get("results", []):
            entries.extend(group.get("entries", []))
        return entries
    except Exception as e:
        print(f"[WARN] Reactome search failed for '{query}': {e}")
        return []


def reactome_cross_refs(st_id: str) -> List[str]:
    """
    获取 Reactome 通路的 KEGG cross-reference。
    返回 KEGG ID 列表（如 hsa04010）。
    Reactome crossReference 端点返回 JSON 列表，每项含 database + identifier。
    """
    url = f"{REACTOME_BASE}/data/query/{st_id}/crossReference"
    try:
        resp = requests.get(url, headers={"accept": "application/json"}, timeout=TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        kegg_ids = []
        for item in data:
            db_name = item.get("database", {}).get("displayName", "")
            if "KEGG" in db_name.upper():
                ident = item.get("identifier", "")
                if ident and ident.startswith("hsa"):
                    kegg_ids.append(ident)
        return kegg_ids
    except Exception as e:
        print(f"[WARN] CrossRef failed for {st_id}: {e}")
        return []


# ---------------------------------------------------------------------------
# KEGG API：ID -> 名称映射
# ---------------------------------------------------------------------------

_kegg_id_to_name_cache: Dict[str, Optional[str]] = {}


def kegg_id_to_name(kegg_id: str) -> Optional[str]:
    """将 KEGG ID（如 hsa04010）转换为通路名称（如 MAPK signaling pathway）。"""
    if kegg_id in _kegg_id_to_name_cache:
        return _kegg_id_to_name_cache[kegg_id]

    url = f"{KEGG_BASE}/get/{kegg_id}"
    try:
        resp = requests.get(url, timeout=TIMEOUT)
        resp.raise_for_status()
        # 解析 NAME 行：NAME        MAPK signaling pathway - Homo sapiens (human)
        for line in resp.text.split("\n"):
            if line.startswith("NAME"):
                # 去掉 "NAME" 前缀，去掉物种后缀
                raw = line.replace("NAME", "", 1).strip()
                name = raw.split(" - ")[0].strip()
                _kegg_id_to_name_cache[kegg_id] = name
                return name
    except Exception as e:
        print(f"[WARN] KEGG ID lookup failed for {kegg_id}: {e}")

    _kegg_id_to_name_cache[kegg_id] = None
    return None


def clean_reactome_name(name: str) -> str:
    """去除 Reactome 返回的 HTML 高亮标签。"""
    return html.unescape(re.sub(r'<[^>]+>', '', name))


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Reactome External Validation: Oxidative Stress Pathway Classification")
    print("=" * 70)

    # =====================================================================
    # Step 1: Reactome 搜索
    # =====================================================================
    print("\n[Step 1] Searching Reactome for oxidative stress pathways...")
    reactome_map: Dict[str, str] = {}  # stId -> clean name

    for term in SEARCH_TERMS:
        entries = reactome_search(term)
        print(f"  '{term}' -> {len(entries)} pathways")
        for entry in entries:
            st_id = entry.get("stId")
            raw_name = entry.get("name", "")
            clean_name = clean_reactome_name(raw_name)
            if st_id and st_id not in reactome_map:
                reactome_map[st_id] = clean_name
        time.sleep(0.5)

    print(f"\n  Total unique Reactome pathways: {len(reactome_map)}")
    if not reactome_map:
        print("  [ERROR] No Reactome pathways found. Check API connectivity.")
        return

    # =====================================================================
    # Step 2: KEGG cross-reference mapping
    # =====================================================================
    print("\n[Step 2] Fetching KEGG cross-references...")
    reactome_to_kegg_ids: Dict[str, List[str]] = {}
    kegg_id_sources: Dict[str, Set[str]] = defaultdict(set)
    all_kegg_ids: Set[str] = set()

    for i, (st_id, name) in enumerate(reactome_map.items(), 1):
        kegg_ids = reactome_cross_refs(st_id)
        reactome_to_kegg_ids[st_id] = kegg_ids
        for kid in kegg_ids:
            kegg_id_sources[kid].add(st_id)
            all_kegg_ids.add(kid)
        if kegg_ids:
            print(f"  [{i}/{len(reactome_map)}] {name[:40]:<40} -> {kegg_ids}")
        time.sleep(0.2)

    mapped_count = sum(1 for v in reactome_to_kegg_ids.values() if v)
    print(f"\n  Reactome pathways with KEGG mapping: {mapped_count}/{len(reactome_map)}")
    print(f"  Unique KEGG IDs found: {len(all_kegg_ids)}")

    if not all_kegg_ids:
        print("\n  [ERROR] No KEGG mappings found.")
        return

    # =====================================================================
    # Step 3: KEGG ID -> Name conversion
    # =====================================================================
    print("\n[Step 3] Converting KEGG IDs to pathway names...")
    kegg_id_to_name_map: Dict[str, str] = {}
    for kid in sorted(all_kegg_ids):
        name = kegg_id_to_name(kid)
        if name:
            kegg_id_to_name_map[kid] = name
            print(f"  {kid} -> {name}")
        else:
            print(f"  {kid} -> [lookup failed]")
        time.sleep(0.3)

    gold_standard_names: Set[str] = set(kegg_id_to_name_map.values())
    print(f"\n  Successfully mapped: {len(gold_standard_names)}/{len(all_kegg_ids)}")

    if not gold_standard_names:
        print("  [ERROR] Could not map any KEGG IDs to names.")
        return

    # =====================================================================
    # Step 4: Load current framework classifications
    # =====================================================================
    print("\n[Step 4] Loading current framework KEGG classifications...")
    stability_file = (
        Path(__file__).resolve().parent.parent.parent
        / "src" / "outputs" / "sensitivity" / "pathway_rank_stability.csv"
    )

    current_classifications: Dict[str, str] = {}
    if stability_file.exists():
        import pandas as pd
        df = pd.read_csv(stability_file)
        kegg_df = df[df["source"] == "KEGG"]
        for _, row in kegg_df.iterrows():
            current_classifications[row["pathway"]] = row["pathway_category"]
        print(f"  Loaded {len(current_classifications)} KEGG pathways from {stability_file.name}")
    else:
        print(f"  [WARN] {stability_file} not found. Will use on-the-fly classification instead.")
        # 如果文件不存在，用 classification.py 实时分类 gold_standard_names 中的通路
        for name in gold_standard_names:
            current_classifications[name] = classify_pathway(name, "KEGG")
        print(f"  Classified {len(current_classifications)} pathways on-the-fly.")

    # =====================================================================
    # Step 5: Consistency metrics
    # =====================================================================
    print("\n[Step 5] Computing consistency metrics...")

    gold = gold_standard_names
    predicted_positive = {
        k for k, v in current_classifications.items()
        if v in ("cross_talk", "signaling_backbone")
    }
    common = gold & set(current_classifications.keys())

    print(f"\n  Reactome gold standard (KEGG names): {len(gold)}")
    print(f"  Current framework total KEGG: {len(current_classifications)}")
    print(f"  Overlap (common): {len(common)}")

    if not common:
        print("\n  [WARN] No overlap between Reactome and current framework.")
        print("  Sample Reactome KEGG names:")
        for name in list(gold)[:10]:
            print(f"    - {name}")
        print("\n  Sample current framework KEGG names:")
        for name in list(current_classifications.keys())[:10]:
            print(f"    - {name}")
        return

    tp = len(common & predicted_positive)
    fp = len(common - predicted_positive)
    jaccard = tp / len(common) if common else 0

    print(f"\n  On the {len(common)} overlapping pathways:")
    print(f"    Reactome + Current agree (relevant): {tp}")
    print(f"    Reactome says relevant, Current says other: {fp}")
    print(f"    Jaccard Overlap: {jaccard:.2%}")

    # =====================================================================
    # Step 6: Detailed comparison
    # =====================================================================
    print("\n[Step 6] Detailed comparison (Reactome gold standard pathways):")
    print(f"\n{'KEGG Pathway':<50} {'Current Framework':<20} {'Match?'}")
    print("-" * 80)

    for name in sorted(common):
        cat = current_classifications.get(name, "N/A")
        match = "YES" if cat in ("cross_talk", "signaling_backbone") else "NO"
        print(f"{name[:48]:<50} {cat:<20} {match}")

    missing_in_current = gold - set(current_classifications.keys())
    if missing_in_current:
        print(f"\n  Reactome pathways NOT in current framework ({len(missing_in_current)}):")
        for name in sorted(missing_in_current)[:20]:
            print(f"    - {name}")

    # =====================================================================
    # Step 7: Save results
    # =====================================================================
    output = {
        "search_terms": SEARCH_TERMS,
        "reactome_pathways_found": len(reactome_map),
        "reactome_with_kegg_mapping": mapped_count,
        "unique_kegg_ids": len(all_kegg_ids),
        "kegg_ids_mapped_to_name": len(gold_standard_names),
        "current_framework_kegg": len(current_classifications),
        "common_pathways": len(common),
        "metrics": {
            "jaccard": round(jaccard, 4),
            "true_positives": tp,
            "false_positives": fp,
        },
        "gold_standard_kegg_names": sorted(gold_standard_names),
        "common_pathway_details": [
            {"pathway": name, "current_category": current_classifications.get(name, "N/A")}
            for name in sorted(common)
        ],
        "reactome_kegg_mapping": {
            st_id: {"name": reactome_map[st_id], "kegg_ids": kids}
            for st_id, kids in reactome_to_kegg_ids.items() if kids
        },
    }

    out_path = (
        Path(__file__).resolve().parent.parent.parent
        / "src" / "outputs" / "sensitivity" / "reactome_validation.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    print(f"\n[OK] Results saved: {out_path}")


if __name__ == "__main__":
    main()
