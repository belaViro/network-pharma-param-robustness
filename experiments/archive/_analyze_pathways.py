import pandas as pd

df = pd.read_csv('src/outputs/sensitivity/mvp_v2_pathway_rank_stability.csv')

print('=' * 70)
print('氧化应激相关通路排名 (按平均排名升序)')
print('=' * 70)

ox = df[df['is_oxidative'] == True].sort_values('avg_rank')
print(f'共 {len(ox)} 条氧化应激相关通路\n')

print(f"{'通路名称':<55} {'排名':<8} {'频率':<8} {'判定'}")
print('-' * 80)

for _, row in ox.iterrows():
    name = row['pathway'][:52]
    stable = '核心稳健' if row['frequency'] >= 0.8 else ('重要' if row['frequency'] >= 0.5 else '潜在')
    print(f"{name:<55} {row['avg_rank']:<8.1f} {row['frequency']:<8.1%} {stable}")

print('\n' + '=' * 70)
print('前20名通路中氧化应激占比')
print('=' * 70)

top20 = df.head(20)
ox_top20 = top20[top20['is_oxidative'] == True]
print(f'前20名中氧化应激通路: {len(ox_top20)} 条 ({len(ox_top20)/20:.0%})')
for _, row in ox_top20.iterrows():
    print(f"  - {row['pathway']} (排名={row['avg_rank']:.1f}, 频率={row['frequency']:.1%})")

print('\n' + '=' * 70)
print('频率>=50%的通路统计')
print('=' * 70)

high_freq = df[df['frequency'] >= 0.5]
ox_high = high_freq[high_freq['is_oxidative'] == True]
non_ox_high = high_freq[high_freq['is_oxidative'] == False]

print(f"频率>=50%的通路总数: {len(high_freq)}")
print(f"  - 氧化应激相关: {len(ox_high)} 条")
print(f"  - 非氧化应激:   {len(non_ox_high)} 条")

print('\n频率>=50%且排名<=20的通路:')
stable = high_freq[high_freq['avg_rank'] <= 20].sort_values('avg_rank')
for _, row in stable.iterrows():
    ox_mark = '[OX]' if row['is_oxidative'] else '    '
    print(f"  {ox_mark} {row['pathway'][:50]} (排名={row['avg_rank']:.1f})")
