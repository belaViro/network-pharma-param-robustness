import json

data = json.load(open(r'D:\科研\网络药理学参数敏感性分析\code\src\outputs\sensitivity\q1_enrichment_audit.json', 'r', encoding='utf-8'))

print('Total pathways:', len(data['pathways']))
print()

kegg = [p for p in data['pathways'] if p['source'] == 'KEGG']
print('KEGG total:', len(kegg))

q1_kegg = [p for p in kegg if p['quadrant'] == 'Q1: 核心稳健通路']
print('KEGG Q1:', len(q1_kegg))

print()
print('KEGG Q3 pathways (稳健但效应弱):')
q3_kegg = [p for p in kegg if p['quadrant'] == 'Q3: 稳健但效应弱']
for p in q3_kegg:
    print(f"  {p['pathway']}: freq={p['frequency']:.2f}, logP={p['avg_logp']:.2f}, cat={p.get('pathway_category','N/A')}")

print()
print('GO_BP Q1 pathways:')
go_bp = [p for p in data['pathways'] if p['source'] == 'GO_BP']
q1_bp = [p for p in go_bp if p['quadrant'] == 'Q1: 核心稳健通路']
for p in q1_bp:
    print(f"  {p['pathway']}: freq={p['frequency']:.2f}, logP={p['avg_logp']:.2f}, cat={p.get('pathway_category','N/A')}")

print()
print('GO_BP Q3 pathways:')
q3_bp = [p for p in go_bp if p['quadrant'] == 'Q3: 稳健但效应弱']
for p in q3_bp:
    print(f"  {p['pathway']}: freq={p['frequency']:.2f}, logP={p['avg_logp']:.2f}, cat={p.get('pathway_category','N/A')}")
