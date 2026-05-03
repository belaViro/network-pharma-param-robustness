import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

# 数据
receptors = ['STAT3', 'EGFR', 'SRC']
ligands = ['Catechin', 'Miquelianin', 'Procyanidin B1']
data = np.array([
    [-6.265, -7.234, -7.024],
    [-8.267, -9.499, -7.795],
    [-9.160, -9.764, -8.310],
])

fig, ax = plt.subplots(figsize=(8, 5.5))

# 热图
cmap = sns.color_palette("YlOrRd_r", as_cmap=True)
im = ax.imshow(data, cmap=cmap, aspect='auto', vmin=-10, vmax=-6)

# 设置刻度
ax.set_xticks(np.arange(len(ligands)))
ax.set_yticks(np.arange(len(receptors)))
ax.set_xticklabels(ligands, fontsize=11)
ax.set_yticklabels(receptors, fontsize=11)

# 旋转x标签
plt.setp(ax.get_xticklabels(), rotation=15, ha="right", rotation_mode="anchor")

# 在每个格子中标注数值
for i in range(len(receptors)):
    for j in range(len(ligands)):
        val = data[i, j]
        text_color = "white" if val < -8.0 else "black"
        ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                color=text_color, fontsize=13, fontweight='bold')

# 颜色条
cbar = ax.figure.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
cbar.ax.set_ylabel('Binding Affinity (kcal/mol)', rotation=-90, va="bottom", fontsize=11)

ax.set_title('Molecular Docking Results: Strawberry Leaf Compounds vs. Dual-Robust Targets',
             fontsize=12, fontweight='bold', pad=15)
ax.set_xlabel('Ligands', fontsize=11, labelpad=8)
ax.set_ylabel('Receptors', fontsize=11, labelpad=8)

# 添加阈值注释
ax.annotate('← Stronger binding', xy=(0.02, 0.02), xycoords='figure fraction',
            fontsize=9, style='italic', color='gray')

plt.tight_layout()
plt.savefig('results/docking_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('results/docking_heatmap.pdf', bbox_inches='tight')
print("热图已保存: docking_heatmap.png / .pdf")
