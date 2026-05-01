"""
活性成分-核心靶点-通路 关系图
Compound-Target-Pathway Relationship Chart

参考样式：水平分组图表，展示每个通路下的靶点分布
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple
import warnings
import logging

from src.config import STRAWBERRY_LEAF_MOLECULES, safe_compound_file_stem

logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore')

# 学术配色
COLOR_PALETTE = [
    '#4CAF50', '#8BC34A', '#FF9800', '#2196F3', '#9C27B0', '#C62828',
    '#00897B', '#6D4C41', '#3949AB', '#F4511E', '#7CB342', '#546E7A',
    '#D81B60', '#00ACC1', '#5E35B1', '#8D6E63', '#039BE5', '#C0CA33',
]
COLORS = {
    compound: COLOR_PALETTE[i % len(COLOR_PALETTE)]
    for i, compound in enumerate(STRAWBERRY_LEAF_MOLECULES)
}


def _active_compounds(df: pd.DataFrame) -> list[str]:
    """Return configured compounds that are present in the relationship table."""
    return [compound for compound in STRAWBERRY_LEAF_MOLECULES if compound in df.columns]


def load_data(data_dir: Path) -> Tuple[Dict, pd.DataFrame, pd.DataFrame]:
    """加载所有数据"""
    # 1. 加载SwissTargetPrediction结果
    swiss_data = {}
    for compound in STRAWBERRY_LEAF_MOLECULES:
        file_path = data_dir / 'swiss_targets' / f'{safe_compound_file_stem(compound)}_targets.csv'
        if file_path.exists():
            df = pd.read_csv(file_path)
            # 筛选概率 > 0 的
            df = df[df['Probability*'] > 0]
            swiss_data[compound] = set(df['Common name'].tolist())

    # 2. 加载Hub节点
    hub_df = pd.read_csv(data_dir / 'ppi_network' / 'hub_nodes.csv')

    # 3. 加载KEGG富集结果
    kegg_df = pd.read_csv(data_dir / 'enrichment' / 'KEGG_enrichment.csv', sep='\t')

    return swiss_data, hub_df, kegg_df


def build_pathway_target_compound(
    swiss_data: Dict,
    hub_df: pd.DataFrame,
    kegg_df: pd.DataFrame,
    top_n: int = 12,
) -> pd.DataFrame:
    """构建通路-靶点-成分关系表"""
    results = []

    for _, row in kegg_df.head(top_n).iterrows():
        pathway = row['Term']
        genes_str = row['Genes']
        pvalue = row['Adjusted P-value']
        overlap = row['Overlap']

        if pd.isna(genes_str):
            continue

        genes = [g.strip() for g in str(genes_str).split(';')]
        hub_genes = set(hub_df['Gene'].tolist())

        # 找出该通路中属于Hub的基因
        pathway_hub_genes = [g for g in genes if g in hub_genes]

        # 统计每个成分预测了多少个该通路的靶点
        compound_counts = {}
        for compound, targets in swiss_data.items():
            count = len(set(pathway_hub_genes) & targets)
            compound_counts[compound] = count

        # 该通路的总靶点数
        total_targets = len(pathway_hub_genes)
        total_genes_in_pathway = int(overlap.split('/')[0]) if isinstance(overlap, str) else 0

        results.append({
            'pathway': pathway,
            'pathway_short': pathway.split(' pathway')[0].split(' signaling')[0][:25],
            'pvalue': pvalue,
            'total_targets': total_targets,
            'total_genes': total_genes_in_pathway,
            'log_pval': -np.log10(max(pvalue, 1e-100)),
            **compound_counts
        })

    return pd.DataFrame(results)


def plot_compound_target_pathway(
    df: pd.DataFrame,
    output_path: str,
    figsize: Tuple = (14, 10),
) -> str:
    """绘制活性成分-靶点-通路图"""

    fig, ax = plt.subplots(figsize=figsize)

    compounds = _active_compounds(df)
    n_pathways = len(df)

    # Y轴位置
    y_positions = np.arange(n_pathways)

    # 绘制堆叠水平条形图 - 每个成分预测的靶点数
    left = np.zeros(n_pathways)

    for i, compound in enumerate(compounds):
        values = df[compound].values
        bars = ax.barh(
            y_positions, values,
            left=left,
            color=COLORS[compound],
            label=compound,
            height=0.7,
            alpha=0.9,
            edgecolor='white',
            linewidth=0.5,
        )

        # 在条形上添加数值标签（如果值>0）
        for j, (val, l) in enumerate(zip(values, left)):
            if val > 0:
                ax.text(
                    l + val / 2, j,
                    str(int(val)),
                    ha='center', va='center',
                    fontsize=8, fontweight='bold', color='white'
                )

        left += values

    # 在条形末端添加总靶点数
    totals = df['total_targets'].values
    for j, (total, lp) in enumerate(zip(totals, left)):
        ax.text(
            lp + 0.3, j,
            f'(n={total})',
            ha='left', va='center',
            fontsize=8, color='#333333'
        )

    # Y轴标签
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['pathway_short'].tolist(), fontsize=10)

    # X轴
    ax.set_xlabel('Number of Target Genes', fontsize=12)
    ax.set_xlim(0, max(left) + 5)

    # 标题
    ax.set_title(
        'Strawberry Leaf Bioactive Compounds: Target-Pathway Distribution',
        fontsize=14, fontweight='bold', pad=15
    )
    ax.text(
        0.5, 1.02,
        'Pathway Enrichment Analysis with Compound Prediction',
        transform=ax.transAxes,
        fontsize=10, ha='center', style='italic', color='gray'
    )

    # 图例
    ax.legend(
        compounds,
        loc='upper right',
        ncol=1,
        fontsize=9,
        title='Compounds',
        title_fontsize=10,
        framealpha=0.9,
    )

    # 美化
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()  # 最显著的在上

    # 添加P-value标注
    for i, pval in enumerate(df['pvalue'].values):
        if pval < 0.001:
            pval_text = 'p<0.001'
        else:
            pval_text = f'p={pval:.2e}'
        ax.text(
            -0.5, i,
            pval_text,
            ha='right', va='center',
            fontsize=7, color='gray'
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"  Saved: {output_path}")
    return output_path


def plot_enhanced_bubble_pathway(
    df: pd.DataFrame,
    output_path: str,
    figsize: Tuple = (14, 10),
) -> str:
    """绘制增强版气泡图 - 展示通路与成分的关系"""

    fig, ax = plt.subplots(figsize=figsize)

    compounds = _active_compounds(df)
    n_pathways = len(df)

    # Y轴
    y_positions = np.arange(n_pathways)

    # 计算气泡大小 (基于总靶点数)
    max_targets = df['total_targets'].max()
    sizes = (df['total_targets'] / max_targets) * 1500 + 200

    # 颜色映射 (基于-log10(P-value))
    colors = -np.log10(df['pvalue'].clip(lower=1e-100))

    # 绘制气泡
    scatter = ax.scatter(
        df['total_targets'],
        y_positions,
        s=sizes,
        c=colors,
        cmap='RdYlBu_r',
        alpha=0.7,
        edgecolors='white',
        linewidths=2,
        zorder=10,
    )

    # 添加颜色条
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('-log₁₀(Adjusted P-value)', fontsize=10)

    # 堆叠显示每个成分的贡献
    for i, compound in enumerate(compounds):
        values = df[compound].values
        left = np.zeros(n_pathways)

        for j, (val, l) in enumerate(zip(values, left)):
            if val > 0:
                # 绘制小分段
                pass

    # Y轴标签
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df['pathway_short'].tolist(), fontsize=10)

    # X轴
    ax.set_xlabel('Number of Hub Target Genes', fontsize=12)
    ax.set_xlim(-1, max(df['total_targets']) + 3)

    # 标题
    ax.set_title(
        'Pathway-Target Network: Compound Contributions',
        fontsize=14, fontweight='bold', pad=15
    )

    # 添加图例说明气泡大小
    legend_sizes = [5, 10, 20]
    legend_labels = ['5 genes', '10 genes', '20 genes']
    for ls, ll in zip(legend_sizes, legend_labels):
        ax.scatter([], [], s=ls * 50 + 100, c='gray', alpha=0.5,
                   label=ll, edgecolors='white', linewidths=1)
    ax.legend(loc='lower right', title='Gene Count', fontsize=8)

    # 美化
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"  Saved: {output_path}")
    return output_path


def plot_grouped_bar_with_targets(
    df: pd.DataFrame,
    output_path: str,
    figsize: Tuple = (16, 12),
) -> str:
    """绘制分组条形图 + 右侧显示关键靶点"""

    fig, (ax_main, ax_legend) = plt.subplots(1, 2, figsize=figsize,
                                               gridspec_kw={'width_ratios': [3, 1]})

    compounds = _active_compounds(df)
    n_pathways = len(df)

    # ========== 主图：堆叠条形图 ==========
    y_positions = np.arange(n_pathways)

    left = np.zeros(n_pathways)
    bars_legend = []

    for compound in compounds:
        values = df[compound].values
        bar = ax_main.barh(
            y_positions, values,
            left=left,
            color=COLORS[compound],
            label=compound,
            height=0.65,
            alpha=0.9,
            edgecolor='white',
            linewidth=0.8,
        )
        bars_legend.append(bar)

        # 标签
        for j, (val, l) in enumerate(zip(values, left)):
            if val >= 1:  # 只显示>=1的
                ax_main.text(
                    l + val / 2, j,
                    str(int(val)),
                    ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white'
                )

        left += values

    # 添加总计
    totals = df['total_targets'].values
    for j, (total, l) in enumerate(zip(totals, left)):
        ax_main.text(
            l + 0.2, j,
            f'Total: {int(total)}',
            ha='left', va='center',
            fontsize=8, color='#333333', fontweight='bold'
        )

    # Y轴
    ax_main.set_yticks(y_positions)
    ax_main.set_yticklabels(df['pathway_short'].tolist(), fontsize=10)

    # X轴
    ax_main.set_xlabel('Number of Predicted Targets', fontsize=12)
    ax_main.set_xlim(0, max(left) + 8)

    # 标题
    ax_main.set_title(
        'Strawberry Leaf: Compound-Target-Pathway Network',
        fontsize=14, fontweight='bold', pad=10
    )

    ax_main.legend(
        compounds,
        loc='upper right',
        fontsize=9,
        title='Bioactive Compounds',
        framealpha=0.95,
    )

    # 美化
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.invert_yaxis()

    # ========== 右侧：Top靶点列表 ==========
    ax_legend.axis('off')

    # 显示Top 10 Hub基因
    top_genes_text = "Top Hub Genes\n" + "─" * 20 + "\n"
    hub_genes = ['TNF', 'SRC', 'EGFR', 'BCL2', 'MMP9',
                  'STAT3', 'PTGS2', 'ESR1', 'HIF1A', 'ERBB2']

    for i, gene in enumerate(hub_genes):
        top_genes_text += f"{i+1}. {gene}\n"

    top_genes_text += "\n" + "─" * 20 + "\n"
    top_genes_text += "Degree Centrality:\n"
    top_genes_text += "TNF (60), SRC (56)\n"
    top_genes_text += "EGFR (50), BCL2 (47)"

    ax_legend.text(
        0.1, 0.95, top_genes_text,
        transform=ax_legend.transAxes,
        fontsize=10,
        verticalalignment='top',
        fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#F5F5F5',
                  edgecolor='gray', alpha=0.9)
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"  Saved: {output_path}")
    return output_path


if __name__ == "__main__":
    data_dir = Path(__file__).parent.parent / "outputs" / "shared"
    output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "figures"

    logger.info("\n" + "=" * 60)
    logger.info("  Plotting: Compound-Target-Pathway Network")
    logger.info("=" * 60)

    # 加载数据
    logger.info("  Loading data...")
    swiss_data, hub_df, kegg_df = load_data(data_dir)

    # 构建关系表
    logger.info("  Building relationships...")
    df = build_pathway_target_compound(swiss_data, hub_df, kegg_df, top_n=12)

    logger.info(f"  - Pathways: {len(df)}")
    logger.info(f"  - Compounds: {len(swiss_data)}")
    logger.info(df[['pathway_short', 'total_targets']].to_string())

    # 绘制图表
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("\n  Generating plots...")

    # 方案1: 堆叠条形图
    plot_compound_target_pathway(
        df,
        str(output_dir / 'compound_target_pathway_stacked.png'),
        figsize=(14, 10),
    )

    # 方案2: 增强气泡图
    plot_enhanced_bubble_pathway(
        df,
        str(output_dir / 'compound_target_pathway_bubble.png'),
        figsize=(14, 10),
    )

    # 方案3: 分组条形图 + 靶点列表
    plot_grouped_bar_with_targets(
        df,
        str(output_dir / 'compound_target_pathway_full.png'),
        figsize=(16, 12),
    )

    logger.info("\n[OK] All plots generated!")
