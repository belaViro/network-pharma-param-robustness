import json

data = json.load(open(r'D:\科研\网络药理学参数敏感性分析\code\src\outputs\sensitivity\q1_enrichment_audit.json', 'r', encoding='utf-8'))

kegg = data['pathway_details']['KEGG']
print('Sample KEGG entry:')
print(json.dumps(kegg[0], indent=2, ensure_ascii=False))
