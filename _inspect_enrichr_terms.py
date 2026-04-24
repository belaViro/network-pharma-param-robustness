import json
import glob

files = glob.glob('src/outputs/sensitivity/enrichment_cache/*_enrichment.json')
print(f"找到 {len(files)} 个缓存文件")

all_terms = set()
for f in sorted(files)[:5]:
    with open(f, 'r', encoding='utf-8') as fh:
        data = json.load(fh)
    for db_name, df_json in data.get("significant", {}).items():
        if isinstance(df_json, str):
            import pandas as pd
            df = pd.read_json(df_json)
            for term in df.get("Term", []):
                all_terms.add(term)

print(f"\n前5个缓存文件中提取了 {len(all_terms)} 条唯一通路")
print("\n通路名称样例（前30条）:")
for i, term in enumerate(sorted(all_terms)[:30]):
    print(f"  {i+1}. {term}")
