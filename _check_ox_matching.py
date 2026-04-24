import pandas as pd
from src.config import ANTIOXIDANT_KEYWORDS

df = pd.read_csv('src/outputs/sensitivity/mvp_v2_pathway_rank_stability.csv')

print("关键词列表:")
print(ANTIOXIDANT_KEYWORDS)
print()

print("被标记为'氧化应激'的通路,及可能匹配的关键词:")
print("=" * 80)

ox = df[df['is_oxidative'] == True].sort_values('avg_rank')
for _, row in ox.iterrows():
    term = row['pathway'].lower()
    matched = [kw for kw in ANTIOXIDANT_KEYWORDS if kw.lower() in term]
    print(f"\n通路: {row['pathway']}")
    print(f"  匹配到的关键词: {matched}")
