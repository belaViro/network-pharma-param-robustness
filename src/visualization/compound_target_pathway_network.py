"""
活性成分-抗氧化核心靶点-通路 网络图
Compound-Target-Pathway Network Visualization

三层网络布局:
Layer 1 (顶部): 草莓叶活性成分
Layer 2 (中部): 核心靶点 (Hub Genes from PPI)
Layer 3 (底部): 富集通路 (KEGG/GO)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings
import logging

from src.config import STRAWBERRY_LEAF_MOLECULES, safe_compound_file_stem

logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore')


# 配色方案 - 学术风格
COLORS = {
    'compound': '#2E7D32',       # 深绿 - 活性成分
    'compound_bg': '#E8F5E9',    # 浅绿背景
    'target_high': '#C62828',    # 深红 - 高重要性靶点
    'target_med': '#EF5350',     # 中红
    'target_low': '#FFCDD2',    # 浅红 - 低重要性靶点
    'pathway': '#1565C0',        # 深蓝 - 通路
    'pathway_bg': '#E3F2FD',     # 浅蓝背景
    'edge_compound': '#4CAF50',  # 成分到靶点边
    'edge_target': '#2196F3',    # 靶点到通路边
    'text_dark': '#212121',
    'text_light': '#FFFFFF',
    'grid': '#E0E0E0',
}


def load_swiss_targets(data_dir: Path) -> Dict[str, pd.DataFrame]:
    """加载各化合物的靶点预测结果"""
    targets_dict = {}

    for compound in STRAWBERRY_LEAF_MOLECULES:
        file_path = data_dir / 'swiss_targets' / f'{safe_compound_file_stem(compound)}_targets.csv'
        if file_path.exists():
            df = pd.read_csv(file_path)
            targets_dict[compound] = df

    return targets_dict


def load_hub_nodes(data_dir: Path, top_n: int = 20) -> pd.DataFrame:
    """加载Hub节点（按度中心性排序）"""
    hub_file = data_dir / 'ppi_network' / 'hub_nodes.csv'
    if hub_file.exists():
        df = pd.read_csv(hub_file)
        return df.head(top_n)
    return pd.DataFrame()


def load_enrichment(data_dir: Path, top_n: int = 10) -> pd.DataFrame:
    """加载KEGG富集结果"""
    kegg_file = data_dir / 'enrichment' / 'KEGG_enrichment.csv'
    if kegg_file.exists():
        df = pd.read_csv(kegg_file, sep='\t')
        return df.head(top_n)
    return pd.DataFrame()


def build_network_data(
    swiss_targets: Dict[str, pd.DataFrame],
    hub_nodes: pd.DataFrame,
    pathways: pd.DataFrame,
    prob_threshold: float = 0.1,
) -> Tuple[Dict, List, List]:
    """构建网络数据

    Returns:
        compound_target_edges: {(compound, gene): weight}
        target_pathway_edges: {(gene, pathway): weight}
        all_nodes: {'compounds': [], 'targets': [], 'pathways': []}
    """
    compound_target_edges = {}
    target_pathway_edges = {}
    all_nodes = {'compounds': [], 'targets': [], 'pathways': []}

    # 1. 化合物 -> 靶点 边
    hub_genes = set(hub_nodes['Gene'].tolist()) if len(hub_nodes) > 0 else set()

    for compound_display, df in swiss_targets.items():
        # 使用 'Common name' 作为基因名列
        if 'Common name' in df.columns and 'Probability*' in df.columns:
            for _, row in df.iterrows():
                gene = row['Common name']
                prob = float(row['Probability*'])
                if prob >= prob_threshold and gene in hub_genes:
                    compound_target_edges[(compound_display, gene)] = prob

    # 2. 靶点 -> 通路 边 (从富集结果解析)
    for _, row in pathways.iterrows():
        term = row['Term']
        genes_str = row.get('Genes', '')
        pval = float(row['Adjusted P-value'])

        # 简化通路名称
        pathway_short = term.split(' pathway')[0].split(' signaling')[0]
        pathway_short = pathway_short.split(' in cancer')[0].strip()

        if 'Genes' in row and pd.notna(genes_str):
            genes = str(genes_str).split(';')
            for gene in genes:
                gene = gene.strip()
                if gene in hub_genes:
                    # 使用-pval作为权重（越小越显著权重越大）
                    weight = -np.log10(max(pval, 1e-100))
                    if (gene, pathway_short) not in target_pathway_edges:
                        target_pathway_edges[(gene, pathway_short)] = weight

    # 节点列表
    all_nodes['compounds'] = list(set([c for c, _ in compound_target_edges.keys()]))
    all_nodes['targets'] = list(set([t for _, t in compound_target_edges.keys()]))
    all_nodes['pathways'] = list(set([p for _, p in target_pathway_edges.keys()]))

    return compound_target_edges, target_pathway_edges, all_nodes


def calculate_node_importance(
    hub_nodes: pd.DataFrame,
    compound_target_edges: Dict,
    target_pathway_edges: Dict,
) -> Dict[str, float]:
    """计算节点重要性得分"""
    importance = {}

    if len(hub_nodes) > 0:
        max_degree = hub_nodes['Degree'].max()
        min_degree = hub_nodes['Degree'].min()

        for _, row in hub_nodes.iterrows():
            gene = row['Gene']
            # 归一化度中心性
            norm_deg = (row['Degree'] - min_degree) / (max_degree - min_degree + 1)
            # 计算化合物连接数
            compound_count = sum(1 for (c, g), _ in compound_target_edges.items() if g == gene)
            # 计算通路连接数
            pathway_count = sum(1 for (g, p) in target_pathway_edges.keys() if g == gene)
            # 综合得分
            n_compounds = max(len({c for c, _ in compound_target_edges.keys()}), 1)
            importance[gene] = 0.5 * norm_deg + 0.3 * (compound_count / n_compounds) + 0.2 * (pathway_count / 10)

    return importance


def draw_layer_background(
    ax: plt.Axes,
    y: float,
    height: float,
    color: str,
    alpha: float = 0.15,
) -> None:
    """绘制层级背景区域"""
    ax.add_patch(FancyBboxPatch(
        (-0.5, y - height / 2),
        10, height,
        boxstyle="round,pad=0.02,rounding_size=0.1",
        facecolor=color,
        edgecolor='none',
        alpha=alpha,
        transform=ax.transData,
        zorder=0,
    ))


def draw_node(
    ax: plt.Axes,
    x: float,
    y: float,
    size: float,
    color: str,
    label: str,
    shape: str = 'circle',
    fontsize: int = 8,
) -> None:
    """绘制节点"""
    if shape == 'circle':
        circle = plt.Circle(
            (x, y), size,
            facecolor=color,
            edgecolor='white',
            linewidth=1.5,
            zorder=10,
        )
        ax.add_patch(circle)
        # 标签
        ax.text(x, y, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='white', zorder=11)
    elif shape == 'ellipse':
        ellipse = mpatches.Ellipse(
            (x, y), size * 1.5, size,
            facecolor=color,
            edgecolor='white',
            linewidth=1.5,
            zorder=10,
        )
        ax.add_patch(ellipse)
        ax.text(x, y, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='white', zorder=11)


def draw_edge(
    ax: plt.Axes,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    color: str,
    alpha: float = 0.5,
    linewidth: float = 1.0,
    style: str = 'solid',
) -> None:
    """绘制边"""
    linestyle = '-' if style == 'solid' else '--'
    ax.plot(
        [x1, x2], [y1, y2],
        color=color,
        alpha=alpha,
        linewidth=linewidth,
        linestyle=linestyle,
        zorder=5,
    )


def draw_three_layer_network(
    compound_target_edges: Dict,
    target_pathway_edges: Dict,
    all_nodes: Dict,
    importance: Dict[str, float],
    output_path: str = None,
    figsize: Tuple = (16, 12),
) -> str:
    """绘制三层网络图"""

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1, 13)
    ax.set_aspect('equal')
    ax.axis('off')

    # ========== 层级高度定义 ==========
    LAYER_COMPOUND = 11      # 化合物层
    LAYER_TARGET = 6.5       # 靶点层
    LAYER_PATHWAY = 2         # 通路层

    # ========== 绘制层级背景 ==========
    draw_layer_background(ax, LAYER_COMPOUND, 1.5, COLORS['compound_bg'])
    draw_layer_background(ax, LAYER_TARGET, 3.5, COLORS['target_low'], alpha=0.1)
    draw_layer_background(ax, LAYER_PATHWAY, 2.5, COLORS['pathway_bg'])

    # ========== 绘制层级标题 ==========
    ax.text(0.2, LAYER_COMPOUND, '🍓 Bioactive Compounds',
            fontsize=11, fontweight='bold', color=COLORS['compound'], va='center')
    ax.text(0.2, LAYER_TARGET, '🎯 Core Targets (Hub Genes)',
            fontsize=11, fontweight='bold', color='#C62828', va='center')
    ax.text(0.2, LAYER_PATHWAY, '📊 Enriched Pathways',
            fontsize=11, fontweight='bold', color=COLORS['pathway'], va='center')

    # ========== 计算节点位置 ==========

    # 化合物节点 (顶部横排)
    compounds = sorted(all_nodes['compounds'])
    n_compounds = len(compounds)
    compound_x = {c: 2 + i * (8 / max(n_compounds - 1, 1)) for i, c in enumerate(compounds)}

    # 靶点节点 (中间层，按重要性从左到右)
    targets = sorted(all_nodes['targets'], key=lambda x: importance.get(x, 0), reverse=True)
    n_targets = len(targets)
    target_x = {t: 1 + i * (8 / max(n_targets - 1, 1)) for i, t in enumerate(targets)}

    # 通路节点 (底部，按相似性分组)
    pathways = sorted(all_nodes['pathways'])
    n_pathways = len(pathways)
    pathway_x = {p: 1 + i * (8 / max(n_pathways - 1, 1)) for i, p in enumerate(pathways)}

    # ========== 绘制边 (先画边再画节点) ==========

    # 化合物 -> 靶点 边
    for (compound, target), weight in compound_target_edges.items():
        if compound in compound_x and target in target_x:
            alpha = 0.3 + 0.5 * weight
            linewidth = 0.5 + 1.5 * weight
            draw_edge(
                ax,
                compound_x[compound], LAYER_COMPOUND - 0.4,
                target_x[target], LAYER_TARGET + 0.4,
                COLORS['edge_compound'],
                alpha=alpha,
                linewidth=linewidth,
                style='solid',
            )

    # 靶点 -> 通路 边
    for (target, pathway), weight in target_pathway_edges.items():
        if target in target_x and pathway in pathway_x:
            alpha = min(0.8, 0.3 + 0.1 * weight)
            linewidth = 0.3 + 0.3 * min(weight / 5, 1)
            draw_edge(
                ax,
                target_x[target], LAYER_TARGET - 0.4,
                pathway_x[pathway], LAYER_PATHWAY + 0.3,
                COLORS['edge_target'],
                alpha=alpha,
                linewidth=linewidth,
                style='--',
            )

    # ========== 绘制节点 ==========

    # 化合物节点 (大椭圆)
    compound_sizes = {
        'Catechin': 0.45,
        'Procyanidin B1': 0.45,
        'Scutellarin': 0.4,
        'Sucrose': 0.35,
        'Trehalose': 0.35,
    }
    for compound in compounds:
        size = compound_sizes.get(compound, 0.4)
        draw_node(
            ax, compound_x[compound], LAYER_COMPOUND,
            size, COLORS['compound'],
            compound.replace(' ', '\n'), shape='ellipse', fontsize=7
        )

    # 靶点节点 (圆形，大小按重要性)
    min_size, max_size = 0.25, 0.45
    imp_values = [importance.get(t, 0) for t in targets]
    imp_min, imp_max = min(imp_values) if imp_values else 0, max(imp_values) if imp_values else 1

    for target in targets:
        imp = importance.get(target, 0)
        if imp_max > imp_min:
            size = min_size + (imp - imp_min) / (imp_max - imp_min) * (max_size - min_size)
        else:
            size = 0.3

        # 颜色按重要性
        if imp > 0.7:
            color = COLORS['target_high']
        elif imp > 0.4:
            color = COLORS['target_med']
        else:
            color = '#EF9A9A'

        draw_node(ax, target_x[target], LAYER_TARGET, size, color, target, fontsize=6)

    # 通路节点 (圆角矩形，用长条表示)
    for pathway in pathways:
        # 计算该通路有多少靶点连接
        connected_targets = sum(1 for (t, p) in target_pathway_edges.keys() if p == pathway)
        # 高度反映连接数
        box_height = 0.25 + 0.15 * min(connected_targets, 5)

        rect = FancyBboxPatch(
            (pathway_x[pathway] - 0.4, LAYER_PATHWAY - box_height / 2),
            0.8, box_height,
            boxstyle="round,pad=0.02,rounding_size=0.05",
            facecolor=COLORS['pathway'],
            edgecolor='white',
            linewidth=1,
            zorder=10,
        )
        ax.add_patch(rect)

        # 简短的标签
        label = pathway[:12] + '...' if len(pathway) > 12 else pathway
        ax.text(pathway_x[pathway], LAYER_PATHWAY, label,
                ha='center', va='center',
                fontsize=6, fontweight='bold', color='white', zorder=11)

    # ========== 添加图例 ==========
    legend_x = 9.5
    legend_y = 10

    # 节点类型图例
    ax.add_patch(plt.Circle((legend_x, legend_y), 0.15, facecolor=COLORS['compound'], zorder=10))
    ax.text(legend_x + 0.3, legend_y, 'Compound', fontsize=8, va='center')

    ax.add_patch(plt.Circle((legend_x, legend_y - 0.5), 0.15, facecolor=COLORS['target_med'], zorder=10))
    ax.text(legend_x + 0.3, legend_y - 0.5, 'Target (Hub)', fontsize=8, va='center')

    ax.add_patch(FancyBboxPatch((legend_x - 0.15, legend_y - 1.0 - 0.1), 0.3, 0.2,
                                boxstyle="round,pad=0.02,rounding_size=0.05",
                                facecolor=COLORS['pathway'], zorder=10))
    ax.text(legend_x + 0.3, legend_y - 1.0, 'Pathway', fontsize=8, va='center')

    # 边类型图例
    ax.plot([legend_x - 0.15, legend_x + 0.15], [legend_y - 2, legend_y - 2],
            color=COLORS['edge_compound'], linewidth=2, alpha=0.7)
    ax.text(legend_x + 0.3, legend_y - 2, 'Compound→Target', fontsize=8, va='center')

    ax.plot([legend_x - 0.15, legend_x + 0.15], [legend_y - 2.5, legend_y - 2.5],
            color=COLORS['edge_target'], linewidth=2, linestyle='--', alpha=0.7)
    ax.text(legend_x + 0.3, legend_y - 2.5, 'Target→Pathway', fontsize=8, va='center')

    # 节点大小图例
    ax.text(legend_x, legend_y - 3.5, 'Node size:', fontsize=7, fontweight='bold', va='center')
    ax.add_patch(plt.Circle((legend_x - 0.1, legend_y - 4), 0.1, facecolor='gray', alpha=0.5))
    ax.text(legend_x + 0.3, legend_y - 4, 'Low importance', fontsize=7, va='center')
    ax.add_patch(plt.Circle((legend_x - 0.1, legend_y - 4.5), 0.2, facecolor='gray', alpha=0.7))
    ax.text(legend_x + 0.3, legend_y - 4.5, 'High importance', fontsize=7, va='center')

    # ========== 添加统计信息 ==========
    stats_text = (
        f"Network Statistics:\n"
        f"• Compounds: {len(compounds)}\n"
        f"• Targets: {len(targets)}\n"
        f"• Pathways: {len(pathways)}\n"
        f"• C→T edges: {len(compound_target_edges)}\n"
        f"• T→P edges: {len(target_pathway_edges)}"
    )
    ax.text(0.5, 0.5, stats_text, fontsize=8, va='bottom', ha='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'),
            transform=ax.transData)

    # ========== 标题 ==========
    ax.text(5, 12.5, 'Strawberry Leaf: Compound-Target-Pathway Network',
            fontsize=14, fontweight='bold', ha='center', color='#1A237E')
    ax.text(5, 12.0, 'Network Pharmacology Analysis of Antioxidant Mechanisms',
            fontsize=10, ha='center', color='gray', style='italic')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        logger.info(f"  保存: {output_path}")

    plt.close()
    return output_path


def create_network_figure(data_dir: Path, output_dir: Path) -> str:
    """创建完整的网络图"""

    # 1. 加载数据
    logger.info("  加载数据...")
    swiss_targets = load_swiss_targets(data_dir)
    hub_nodes = load_hub_nodes(data_dir, top_n=20)
    pathways = load_enrichment(data_dir, top_n=12)

    if len(hub_nodes) == 0 or len(pathways) == 0:
        logger.info("  错误: 缺少Hub节点或通路数据")
        return None

    logger.info(f"  - 化合物: {len(swiss_targets)}种")
    logger.info(f"  - Hub靶点: {len(hub_nodes)}个")
    logger.info(f"  - 富集通路: {len(pathways)}条")

    # 2. 构建网络数据
    logger.info("  构建网络数据...")
    compound_target_edges, target_pathway_edges, all_nodes = build_network_data(
        swiss_targets, hub_nodes, pathways, prob_threshold=0.1
    )

    logger.info(f"  - 成分-靶点边: {len(compound_target_edges)}条")
    logger.info(f"  - 靶点-通路边: {len(target_pathway_edges)}条")

    # 3. 计算节点重要性
    importance = calculate_node_importance(hub_nodes, compound_target_edges, target_pathway_edges)

    # 4. 绘制网络图
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / 'compound_target_pathway_network.png'

    logger.info("  绘制网络图...")
    draw_three_layer_network(
        compound_target_edges,
        target_pathway_edges,
        all_nodes,
        importance,
        output_path=str(output_path),
        figsize=(18, 14),
    )

    return str(output_path)


if __name__ == "__main__":
    # 测试
    data_dir = Path(__file__).parent.parent / "outputs" / "shared"
    output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "figures"

    logger.info("\n" + "=" * 60)
    logger.info("  绘制: 活性成分-核心靶点-通路 网络图")
    logger.info("=" * 60)

    result = create_network_figure(data_dir, output_dir)

    if result:
        logger.info(f"\n[OK] 图表已生成: {result}")
    else:
        logger.info("\n[ERROR] 图表生成失败")
