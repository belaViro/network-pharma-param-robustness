import pandas as pd
from src.config import ANTIOXIDANT_KEYWORDS

df = pd.read_csv('src/outputs/sensitivity/mvp_v2_pathway_rank_stability.csv')

print("关键词列表:", ANTIOXIDANT_KEYWORDS)
print()

ox = df[df['is_oxidative'] == True].sort_values('avg_rank')
print(f"总共标记为氧化应激的通路数: {len(ox)}")
print()

for _, row in ox.iterrows():
    term = row['pathway'].lower()
    matched = [kw for kw in ANTIOXIDANT_KEYWORDS if kw.lower() in term]
    print(f"通路: {row['pathway']}")
    print(f"  匹配到的关键词: {matched}")
    print()
