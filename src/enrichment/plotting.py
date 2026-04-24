"""
Enrichment Plotting Module
绘制GO和KEGG富集分析图表
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict
import logging


logger = logging.getLogger(__name__)

def plot_go_enrichment(
    go_df,
    output_path: str = None,
    title: str = "GO Enrichment",
    figsize: tuple = (12, 8),
    top_n: int = 10,
) -> str:
    """绘制GO富集分析条形图

    Args:
        go_df: GO富集结果DataFrame
        output_path: 输出路径
        title: 图表标题
        figsize: 图像尺寸
        top_n: 显示前N条

    Returns:
        保存的图像路径
    """
    if output_path is None:
        output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "figures"
        output_path = output_dir / f"{title.lower().replace(' ', '_')}_barplot.png"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if go_df is None or len(go_df) == 0:
        logger.info("  无数据可绘图")
        return None

    df = go_df.head(top_n).copy()
    df = df.sort_values('Adjusted P-value')

    fig, ax = plt.subplots(figsize=figsize)

    y_pos = range(len(df))
    log_pval = -df['Adjusted P-value'].apply(lambda x: max(0.001, x)).apply(lambda x: -np.log10(x))

    colors = plt.cm.viridis(log_pval / log_pval.max())

    bars = ax.barh(y_pos, log_pval, color=colors, alpha=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([term[:60] + '...' if len(term) > 60 else term for term in df['Term']], fontsize=9)
    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11)
    ax.set_title(f'{title}\n(Top {top_n} by P-value)', fontsize=13, fontweight='bold')

    for i, (idx, row) in enumerate(df.iterrows()):
        ax.text(log_pval.iloc[i] + 0.1, i, f"p={row['Adjusted P-value']:.2e}",
                va='center', fontsize=8, color='#333333')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"  保存: {output_path}")
    return str(output_path)


def plot_kegg_enrichment(
    kegg_df,
    output_path: str = None,
    title: str = "KEGG Pathway Enrichment",
    figsize: tuple = (12, 10),
    top_n: int = 20,
) -> str:
    """绘制KEGG富集分析气泡图

    Args:
        kegg_df: KEGG富集结果DataFrame (需包含EF列)
        output_path: 输出路径
        title: 图表标题
        figsize: 图像尺寸
        top_n: 显示前N条

    Returns:
        保存的图像路径
    """
    if output_path is None:
        output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "figures"
        output_path = output_dir / "kegg_bubble_plot.png"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if kegg_df is None or len(kegg_df) == 0:
        logger.info("  无数据可绘图")
        return None

    df = kegg_df.nsmallest(top_n, 'Adjusted P-value').copy()

    df['Count'] = df['Overlap'].apply(
        lambda x: int(x.split('/')[0]) if isinstance(x, str) else int(x)
    )

    log_pval = -np.log10(df['Adjusted P-value'].clip(lower=1e-100))

    min_size, max_size = 200, 1500
    counts = df['Count'].values.astype(float)
    if counts.max() != counts.min():
        sizes = min_size + (counts - counts.min()) / (counts.max() - counts.min()) * (max_size - min_size)
    else:
        sizes = np.full(len(df), (min_size + max_size) / 2)

    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        df['EF'],
        range(len(df)),
        s=sizes,
        c=log_pval,
        cmap='RdYlBu_r',
        alpha=0.75,
        edgecolors='white',
        linewidths=1.0,
    )

    y_labels = df['Term'].apply(lambda x: x[:45] + '...' if len(x) > 45 else x).tolist()
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(y_labels, fontsize=9)

    ax.set_xlabel('Enrichment', fontsize=11)
    ax.set_ylabel('')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=15)
    ax.invert_yaxis()
    ax.set_xlim(left=0)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('-log₁₀(Adjusted P-value)', fontsize=10)

    for i, (idx, row) in enumerate(df.iterrows()):
        fontsize = 9 if int(row['Count']) >= 10 else 10
        ax.text(
            df['EF'].iloc[i],
            i,
            f"{int(row['Count'])}",
            va='center',
            ha='center',
            fontsize=fontsize,
            color='white',
            fontweight='bold',
        )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(left=False)
    ax.grid(axis='x', alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"  保存: {output_path}")
    return str(output_path)


def plot_go_grouped_barplot(
    go_results: Dict,
    output_path: str = None,
    top_n: int = 10,
    figsize: tuple = (18, 6),
) -> str:
    """绘制 GO 富集分析分组条形图（分面设计：BP、MF、CC 三列并排）

    Args:
        go_results: 包含 GO_BP、GO_MF、GO_CC 的字典
        output_path: 输出路径
        top_n: 每个类别取前 N 条
        figsize: 图像尺寸 (总宽度, 总高度)

    Returns:
        保存的图像路径
    """
    if output_path is None:
        output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "figures"
        output_path = output_dir / "go_grouped_barplot.png"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    categories = {
        'GO_BP': {'label': 'Biological Process', 'color': '#E74C3C'},
        'GO_MF': {'label': 'Molecular Function', 'color': '#3498DB'},
        'GO_CC': {'label': 'Cellular Component', 'color': '#27AE60'},
    }

    valid_cats = {k: v for k, v in categories.items()
                  if k in go_results and len(go_results[k]) > 0}

    if not valid_cats:
        logger.info("  无数据可绘图")
        return None

    n_cats = len(valid_cats)
    fig, axes = plt.subplots(1, n_cats, figsize=figsize, sharey=False)

    if n_cats == 1:
        axes = [axes]

    for ax, (cat_key, config) in zip(axes, valid_cats.items()):
        df = go_results[cat_key].head(top_n).copy()
        df = df.sort_values('Adjusted P-value')

        log_pval = -np.log10(df['Adjusted P-value'].clip(lower=1e-100))

        labels = df['Term'].apply(lambda x: x.split(' (GO:')[0]).tolist()
        labels = [t[:40] + '...' if len(t) > 40 else t for t in labels]

        y_pos = range(len(df))

        bars = ax.barh(y_pos, log_pval, color=config['color'], alpha=0.85,
                       edgecolor='white', linewidth=0.5, height=0.7)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()

        ax.set_xlabel('-log₁₀(Adjusted P-value)', fontsize=10)
        ax.set_title(f"{config['label']}", fontsize=12, fontweight='bold', color=config['color'])

        for bar, pval in zip(bars, df['Adjusted P-value']):
            ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                    f'p={pval:.1e}', va='center', ha='left', fontsize=7, color='#333333')

        ax.xaxis.grid(True, linestyle='--', alpha=0.4)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('GO Enrichment Analysis', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"  保存: {output_path}")
    return str(output_path)


def plot_all_enrichment(
    significant_results: Dict,
    output_dir: str = None,
) -> Dict:
    """绘制所有富集分析图

    Args:
        significant_results: run_enrichment_analysis返回的significant结果
        output_dir: 输出目录

    Returns:
        图像路径字典
    """
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "figures"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}

    if 'GO_BP' in significant_results and len(significant_results['GO_BP']) > 0:
        paths['GO_BP'] = plot_go_enrichment(
            significant_results['GO_BP'],
            output_dir / 'go_biological_process.png',
            title='GO Biological Process',
            top_n=10,
        )

    if 'GO_MF' in significant_results and len(significant_results['GO_MF']) > 0:
        paths['GO_MF'] = plot_go_enrichment(
            significant_results['GO_MF'],
            output_dir / 'go_molecular_function.png',
            title='GO Molecular Function',
            top_n=10,
        )

    if 'GO_CC' in significant_results and len(significant_results['GO_CC']) > 0:
        paths['GO_CC'] = plot_go_enrichment(
            significant_results['GO_CC'],
            output_dir / 'go_cellular_component.png',
            title='GO Cellular Component',
            top_n=10,
        )

    paths['GO_Grouped'] = plot_go_grouped_barplot(
        significant_results,
        output_dir / 'go_grouped_barplot.png',
        top_n=8,
    )

    if 'KEGG' in significant_results and len(significant_results['KEGG']) > 0:
        paths['KEGG'] = plot_kegg_enrichment(
            significant_results['KEGG'],
            output_dir / 'kegg_bubble_plot.png',
            title='KEGG Pathway Enrichment',
            top_n=20,
        )

    return paths


__all__ = [
    "plot_go_enrichment",
    "plot_kegg_enrichment",
    "plot_go_grouped_barplot",
    "plot_all_enrichment",
]