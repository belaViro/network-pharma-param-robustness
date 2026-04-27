#!/usr/bin/env python3
"""检查靶点和通路的频率分布"""
import pandas as pd

# 靶点频率分布
tdf = pd.read_csv('src/outputs/sensitivity/target_robustness.csv')
print('=== 靶点频率分布 ===')
vc = tdf['frequency'].value_counts().sort_index(ascending=False)
for freq, count in vc.head(15).items():
    print(f'  freq={freq:.4f}: {count} 个靶点')
print(f'唯一频率值数量: {tdf["frequency"].nunique()} / {len(tdf)} 个靶点')

# 通路频率分布
pdf = pd.read_csv('src/outputs/sensitivity/pathway_rank_stability_augmented.csv')
print('\n=== 通路频率分布 ===')
vc2 = pdf['frequency'].value_counts().sort_index(ascending=False)
for freq, count in vc2.head(15).items():
    print(f'  freq={freq:.4f}: {count} 条通路')
print(f'唯一频率值数量: {pdf["frequency"].nunique()} / {len(pdf)} 条通路')
