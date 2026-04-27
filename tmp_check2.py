import json

data = json.load(open(r'D:\科研\网络药理学参数敏感性分析\code\src\outputs\sensitivity\q1_enrichment_audit.json', 'r', encoding='utf-8'))

kegg = data['pathway_details']['KEGG']
print('KEGG total in audit:', len(kegg))

q3 = [p for p in kegg if 'Q3' in p['quadrant']]
print('KEGG Q3 (稳健但效应弱):', len(q3))
print()
for p in q3:
    print(f"  {p['pathway']}: freq={p['frequency']}, logP={p['avg_logp']:.2f}, cat={p.get('pathway_category','N/A')}")

print()
q1 = [p for p in kegg if 'Q1' in p['quadrant']]
print('KEGG Q1 (核心稳健):', len(q1))
for p in q1:
    print(f"  {p['pathway']}: freq={p['frequency']}, logP={p['avg_logp']:.2f}")
