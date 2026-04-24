"""
对比三种氧化应激判断方法在实际数据上的表现：
1. 旧方法：子串匹配（含 'ROS' 误匹配）
2. 改良方法：单词边界匹配
3. GO 层级方法：用 GO:0006979 后代节点
"""
import re
import json
import pandas as pd
from collections import Counter
from src.config import ANTIOXIDANT_KEYWORDS

# 加载 GO 氧化应激集合
go_ox_set = set()
with open('src/data/oxidative_stress_go_terms.txt', 'r', encoding='utf-8') as f:
    for line in f:
        go_id = line.strip().split('\t')[0]
        go_ox_set.add(go_id)

print(f"GO oxidative stress terms loaded: {len(go_ox_set)}")

def old_method(term: str) -> bool:
    """旧方法：子串匹配（有bug）"""
    term_lower = term.lower()
    return any(kw.lower() in term_lower for kw in ANTIOXIDANT_KEYWORDS)

def improved_method(term: str) -> bool:
    """改良方法：单词边界匹配 + 严格关键词"""
    term_lower = term.lower()
    # 严格关键词（去掉 ROS，避免 tyrosine 误匹配）
    strict_kws = [
        'oxidative stress', 'ferroptosis', 'glutathione',
        'peroxisome', 'superoxide', 'lipid peroxidation',
        'catalase', 'peroxidase', 'thioredoxin',
        'autophagy', 'mitophagy', 'NADPH oxidase',
    ]
    for kw in strict_kws:
        if re.search(r'\b' + re.escape(kw) + r'\b', term_lower):
            return True
    # KEGG 专用 curated 列表
    kegg_ox = {
        'AGE-RAGE signaling pathway in diabetic complications',
        'MAPK signaling pathway',
    }
    if term in kegg_ox:
        return True
    return False

def go_method(term: str) -> bool:
    """GO 层级方法：提取 GO ID，查是否是 GO:0006979 后代"""
    match = re.search(r'\(GO:(\d+)\)', term)
    if match:
        go_id = f"GO:{match.group(1)}"
        return go_id in go_ox_set
    return False

# 加载一个实际分析的通路结果
# 从 enrichment cache 中提取所有出现过的通路
all_terms = set()
cache_files = [
    'src/outputs/sensitivity/enrichment_cache/030adf2bc022.json',
    'src/outputs/sensitivity/enrichment_cache/0c98816d2046.json',
    'src/outputs/sensitivity/enrichment_cache/1434cb60bc66.json',
]
for cf in cache_files:
    with open(cf, 'r', encoding='utf-8') as f:
        data = json.load(f)
    for db_name, records in data.items():
        if isinstance(records, list):
            for rec in records:
                all_terms.add(rec.get('Term', ''))

all_terms = sorted(t for t in all_terms if t)
print(f"\n实际分析中出现的通路总数: {len(all_terms)}")

# 对比三种方法
old_results = []
imp_results = []
go_results = []

for term in all_terms:
    old_flag = old_method(term)
    imp_flag = improved_method(term)
    go_flag = go_method(term)
    
    if old_flag:
        old_results.append(term)
    if imp_flag:
        imp_results.append(term)
    if go_flag:
        go_results.append(term)

print(f"\n旧方法（子串匹配）标记为 OX 的数量: {len(old_results)}")
print(f"改良方法（单词边界）标记为 OX 的数量: {len(imp_results)}")
print(f"GO 层级方法标记为 OX 的数量: {len(go_results)}")

print("\n" + "="*80)
print("旧方法误标记的通路（被旧方法标记但 GO 方法不标记）:")
false_positives = [t for t in old_results if not go_method(t)]
for t in false_positives:
    matched = [kw for kw in ANTIOXIDANT_KEYWORDS if kw.lower() in t.lower()]
    print(f"  [FALSE+] {t}")
    print(f"     误匹配关键词: {matched}")

print("\n" + "="*80)
print("GO 层级方法正确标记的通路（真正与氧化应激响应相关）:")
for t in go_results:
    print(f"  [TRUE+] {t}")

print("\n" + "="*80)
print("三种方法同时标记的通路（交集，最可靠）:")
intersection = [t for t in all_terms if old_method(t) and improved_method(t) and go_method(t)]
for t in intersection:
    print(f"  [TRUE+] {t}")
