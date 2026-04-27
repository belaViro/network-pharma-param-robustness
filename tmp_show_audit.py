import json

data = json.load(open(r'D:\科研\网络药理学参数敏感性分析\code\src\outputs\sensitivity\q1_enrichment_audit.json', 'r', encoding='utf-8'))

print('='*70)
print('Q1 核心稳健靶点富集分析结果')
print('='*70)
print('Q1 靶点数量:', data['q1_targets_count'], '个')
print('靶点:', ', '.join(data['q1_targets']))
print()
print('富集显著通路总数:', data['enrichment_summary']['total_significant'])
print('  - 直接OX:', data['enrichment_summary']['direct_ox'])
print('  - 交叉调控:', data['enrichment_summary']['cross_talk'])
print('  - 信号主干:', data['enrichment_summary']['backbone'])
print('  - 其他:', data['enrichment_summary']['other'])
print()

print('按来源分层统计:')
for src, stats in data['pathway_details'].items():
    q1 = sum(1 for p in stats if 'Q1' in p['quadrant'])
    q2 = sum(1 for p in stats if 'Q2' in p['quadrant'])
    q3 = sum(1 for p in stats if 'Q3' in p['quadrant'])
    q4 = sum(1 for p in stats if 'Q4' in p['quadrant'])
    new = sum(1 for p in stats if '未' in p['quadrant'] or 'NEW' in p['quadrant'])
    total = len(stats)
    print(f'  {src:8s}: 总数={total:2d} | Q1={q1:2d}({q1/total*100:5.1f}%) | Q2={q2} | Q3={q3:2d} | Q4={q4} | 未出现={new}')

print()
print('系统通路(KEGG+GO_BP) Q1 率:', data['key_metrics']['system_pathway_q1_rate'], '%')
print('全部通路 Q1 率:', data['key_metrics']['overall_q1_rate'], '%')

print()
print('='*70)
print('KEGG 通路明细')
print('='*70)
for p in data['pathway_details']['KEGG']:
    short_q = p['quadrant'].replace('Q1: 核心稳健通路', 'Q1').replace('Q3: 稳健但效应弱', 'Q3')
    print(f'  [{short_q}] {p["pathway"]:55s} cat={p["category"]:18s} freq={p["frequency"]:5.1%}')

print()
print('='*70)
print('GO_BP 通路明细')
print('='*70)
for p in data['pathway_details']['GO_BP']:
    short_q = p['quadrant'].replace('Q1: 核心稳健通路', 'Q1').replace('Q2: 低频高显著', 'Q2').replace('Q3: 稳健但效应弱', 'Q3')
    print(f'  [{short_q}] {p["pathway"]:55s} cat={p["category"]:18s} freq={p["frequency"]:5.1%}')
