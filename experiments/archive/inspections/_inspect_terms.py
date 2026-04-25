import json
import pandas as pd

with open('src/outputs/sensitivity/enrichment_cache/030adf2bc022.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

print("Keys:", list(data.keys()))
for db_name in ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG']:
    if db_name in data and isinstance(data[db_name], list) and len(data[db_name]) > 0:
        df = pd.DataFrame(data[db_name])
        print(f"\n=== {db_name} ({len(df)} terms) ===")
        for term in df['Term'].head(10):
            print(f"  - {term}")
