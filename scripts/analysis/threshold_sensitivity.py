#!/usr/bin/env python3
"""PubMed threshold sensitivity analysis -- offline computation from existing data."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

with open('src/outputs/sensitivity/pubmed_validation.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

details = data['pathway_details']

framework_positive = set()
for d in details:
    if d['framework_positive']:
        framework_positive.add(d['pathway'])

all_pathways = set(d['pathway'] for d in details)

def calc_metrics(threshold):
    pubmed_positive = set(d['pathway'] for d in details if d['pubmed_count'] >= threshold)
    
    tp = len(pubmed_positive & framework_positive)
    fp = len(framework_positive - pubmed_positive)
    fn = len(pubmed_positive - framework_positive)
    tn = len(all_pathways - pubmed_positive - framework_positive)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    jaccard = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0
    
    return {
        'threshold': threshold,
        'pubmed_pos': len(pubmed_positive),
        'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
        'precision': precision, 'recall': recall, 'f1': f1,
        'jaccard': jaccard,
    }

results = [calc_metrics(t) for t in [5, 10, 20]]

print("=" * 80)
print("PubMed Co-occurrence Threshold Sensitivity Analysis")
print("=" * 80)
print()
print("Threshold  PubMed+  TP  FP  FN  TN  Precision  Recall     F1        Jaccard")
print("-" * 80)
for r in results:
    print(f"  n>={r['threshold']:<5}   {r['pubmed_pos']:<7}  {r['tp']:<3} {r['fp']:<3} {r['fn']:<3} {r['tn']:<3} {r['precision']:<9.2%} {r['recall']:<9.2%} {r['f1']:<9.2%} {r['jaccard']:<9.2%}")

print()
print("=" * 80)
print("Per-pathway PubMed counts (sorted, with threshold brackets)")
print("=" * 80)
sorted_details = sorted(details, key=lambda x: x['pubmed_count'])
for d in sorted_details:
    count = d['pubmed_count']
    bracket = ""
    if 5 <= count < 10:
        bracket = "[5-9]   <-- changes at n=5"
    elif 10 <= count < 20:
        bracket = "[10-19] <-- changes at n=10"
    elif count >= 20:
        bracket = "[20+]"
    else:
        bracket = "[0-4]"
    fw = "Y" if d['framework_positive'] else "N"
    print(f"  n={count:>5} {bracket:<30} [FW={fw}] {d['pathway']}")

print()
print("=" * 80)
print("Detailed breakdown: pathways that change status across thresholds")
print("=" * 80)
# Show pathways with counts in the critical range 0-20
critical = [d for d in details if d['pubmed_count'] < 20]
critical.sort(key=lambda x: x['pubmed_count'])
for d in critical:
    count = d['pubmed_count']
    at5 = "POS" if count >= 5 else "neg"
    at10 = "POS" if count >= 10 else "neg"
    at20 = "POS" if count >= 20 else "neg"
    fw = "POS" if d['framework_positive'] else "neg"
    print(f"  n={count:>3} | n>=5:{at5:<4} n>=10:{at10:<4} n>=20:{at20:<4} | FW:{fw:<4} | {d['pathway']}")
