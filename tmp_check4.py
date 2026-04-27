import json, math

data = json.load(open(r'D:\科研\网络药理学参数敏感性分析\code\src\outputs\sensitivity\q1_enrichment_audit.json', 'r', encoding='utf-8'))

kegg = data['pathway_details']['KEGG']
print('KEGG total:', len(kegg))
print()

for p in kegg:
    logp = -math.log10(p['pvalue']) if p['pvalue'] > 0 else 999
    cat = p.get('category', 'N/A')
    print(f"  {p['pathway']}: freq={p['frequency']}, logP={logp:.2f}, quadrant={p['quadrant']}, cat={cat}")
