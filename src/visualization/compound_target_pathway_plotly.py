"""
活性成分-核心靶点-通路 图 (Plotly版本)
更高质量的学术图表
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

from src.config import STRAWBERRY_LEAF_MOLECULES, safe_compound_file_stem


# 化合物配色
COLOR_PALETTE = [
    '#2E7D32', '#66BB6A', '#FF9800', '#2196F3', '#9C27B0', '#C62828',
    '#00897B', '#6D4C41', '#3949AB', '#F4511E', '#7CB342', '#546E7A',
    '#D81B60', '#00ACC1', '#5E35B1', '#8D6E63', '#039BE5', '#C0CA33',
]

COMPOUNDS = STRAWBERRY_LEAF_MOLECULES
COLORS = {compound: COLOR_PALETTE[i % len(COLOR_PALETTE)] for i, compound in enumerate(COMPOUNDS)}


def _active_compounds(df: pd.DataFrame) -> list[str]:
    """Return compounds that are present as columns in the relationship table."""
    return [compound for compound in COMPOUNDS if compound in df.columns]


def load_data(data_dir: Path):
    """加载数据"""
    swiss_data = {}

    for compound_display in COMPOUNDS:
        compound_file = safe_compound_file_stem(compound_display)
        file_path = data_dir / 'swiss_targets' / f'{compound_file}_targets.csv'
        if file_path.exists():
            df = pd.read_csv(file_path)
            df = df[df['Probability*'] > 0]
            swiss_data[compound_display] = set(df['Common name'].tolist())

    hub_df = pd.read_csv(data_dir / 'ppi_network' / 'hub_nodes.csv')
    kegg_df = pd.read_csv(data_dir / 'enrichment' / 'KEGG_enrichment.csv', sep='\t')

    return swiss_data, hub_df, kegg_df


def build_relationship_df(swiss_data, hub_df, kegg_df, top_n=12):
    """构建关系数据框"""
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
        pathway_hub_genes = [g for g in genes if g in hub_genes]

        compound_counts = {}
        for compound, targets in swiss_data.items():
            count = len(set(pathway_hub_genes) & targets)
            compound_counts[compound] = count

        results.append({
            'pathway': pathway,
            'pathway_short': pathway.split(' pathway')[0].split(' signaling')[0],
            'pvalue': pvalue,
            'log_pval': -np.log10(max(pvalue, 1e-100)),
            'total_targets': len(pathway_hub_genes),
            **compound_counts
        })

    return pd.DataFrame(results)


def plot_stacked_bar(df, output_path):
    """堆叠条形图"""
    fig = go.Figure()
    compounds = _active_compounds(df)

    # 添加堆叠条形
    for compound in compounds:
        fig.add_trace(go.Bar(
            name=compound,
            y=df['pathway_short'],
            x=df[compound],
            orientation='h',
            marker_color=COLORS[compound],
            text=df[compound],
            textposition='inside',
            textfont=dict(color='white', size=11, family='Arial Black'),
            hovertemplate=f'<b>{compound}</b><br>' +
                          'Targets: %{x}<br>' +
                          '<extra></extra>',
        ))

    # 添加总计
    fig.add_trace(go.Bar(
        name='Total',
        y=df['pathway_short'],
        x=[0] * len(df),  # 空条占位
        orientation='h',
        marker_color='rgba(0,0,0,0)',
        showlegend=False,
        text=[f"<b>Total: {t}</b>" for t in df['total_targets']],
        textposition='outside',
        textfont=dict(color='#333333', size=11, family='Arial'),
    ))

    # 布局
    fig.update_layout(
        title=dict(
            text='<b>Strawberry Leaf: Compound-Target-Pathway Network</b><br>' +
                 '<span style="font-size:12px; color:gray">Pathway Enrichment Analysis</span>',
            x=0.5,
            font=dict(size=16)
        ),
        barmode='stack',
        yaxis=dict(
            title='',
            autorange='reversed',
            tickfont=dict(size=11),
        ),
        xaxis=dict(
            title='Number of Target Genes',
            tickfont=dict(size=10),
        ),
        legend=dict(
            title='<b>Bioactive Compounds</b>',
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='center',
            x=0.5,
            font=dict(size=11),
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=600,
        width=900,
        margin=dict(l=150, r=80, t=100, b=60),
        font=dict(family='Arial'),
        showlegend=True,
    )

    # 添加P-value标注
    for i, (pval, pathway) in enumerate(zip(df['pvalue'], df['pathway_short'])):
        if pval < 0.001:
            pval_text = 'p<0.001'
        else:
            pval_text = f'p={pval:.2e}'
        fig.add_annotation(
            x=-2,
            y=i,
            text=pval_text,
            showarrow=False,
            font=dict(size=9, color='gray'),
            xanchor='right',
        )

    fig.update_xaxes(range=[-3, df[compounds].sum(axis=1).max() + 5])
    fig.update_yaxes(range=[-1, len(df)])

    fig.write_image(output_path, scale=2)
    logger.info(f"  Saved: {output_path}")
    return output_path


def plot_grouped_dot(df, output_path):
    """点图/气泡图 - 展示成分-靶点贡献"""
    fig = go.Figure()
    compounds = _active_compounds(df)

    for compound in compounds:
        fig.add_trace(go.Scatter(
            name=compound,
            x=df[compound],
            y=df['pathway_short'],
            mode='markers+text',
            marker=dict(
                size=df[compound] * 5 + 10,
                color=COLORS[compound],
                opacity=0.8,
                line=dict(color='white', width=1.5),
            ),
            text=[str(v) if v > 0 else '' for v in df[compound]],
            textposition='middle center',
            textfont=dict(color='white', size=10, family='Arial Black'),
            hovertemplate=f'<b>{compound}</b><br>' +
                          'Targets: %{text}<br>' +
                          '<extra></extra>',
            showlegend=True,
        ))

    # 添加连接线
    for i, pathway in enumerate(df['pathway_short']):
        values = [df[c].iloc[i] for c in compounds]
        total = sum(values)
        if total > 0:
            cumulative = 0
            for compound in compounds:
                val = df[compound].iloc[i]
                if val > 0:
                    cumulative += val
            # 添加总计数
            fig.add_annotation(
                x=total + 0.5,
                y=i,
                text=f'<b>{total}</b>',
                showarrow=False,
                font=dict(size=12, color='#333333'),
                xanchor='left',
            )

    fig.update_layout(
        title=dict(
            text='<b>Compound Contributions to Each Pathway</b><br>' +
                 '<span style="font-size:12px; color:gray">Bubble size indicates target count per compound</span>',
            x=0.5,
            font=dict(size=16)
        ),
        xaxis=dict(
            title='Number of Predicted Targets',
            tickfont=dict(size=10),
            range=[-1, df[compounds].sum(axis=1).max() + 3],
        ),
        yaxis=dict(
            title='',
            autorange='reversed',
            tickfont=dict(size=11),
        ),
        legend=dict(
            title='<b>Compounds</b>',
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='center',
            x=0.5,
        ),
        plot_bgcolor='white',
        height=600,
        width=900,
        margin=dict(l=150, r=100, t=100, b=60),
        font=dict(family='Arial'),
        showlegend=True,
    )

    fig.write_image(output_path, scale=2)
    logger.info(f"  Saved: {output_path}")
    return output_path


def plot_heatmap_style(df, output_path):
    """热力图风格 - 展示成分-通路矩阵"""
    # 构建矩阵
    compounds = _active_compounds(df)
    matrix = df[compounds].values

    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        y=df['pathway_short'],
        x=compounds,
        colorscale=[
            [0, '#FFFFFF'],
            [0.1, '#E8F5E9'],
            [0.5, '#66BB6A'],
            [1, '#1B5E20']
        ],
        showscale=True,
        colorbar=dict(
            title='Targets',
            tickfont=dict(size=10),
        ),
        text=matrix,
        texttemplate='<b>%{text}</b>',
        textfont=dict(color='white', size=12, family='Arial Black'),
        hovertemplate='<b>%{y}</b><br>' +
                      'Compound: %{x}<br>' +
                      'Targets: %{z}<br>' +
                      '<extra></extra>',
    ))

    # 添加总计行
    totals = df[compounds].sum(axis=1).values

    fig.update_layout(
        title=dict(
            text='<b>Compound-Pathway Target Matrix</b><br>' +
                 '<span style="font-size:12px; color:gray">Number of shared target genes</span>',
            x=0.5,
            font=dict(size=16)
        ),
        yaxis=dict(
            title='',
            autorange='reversed',
            tickfont=dict(size=10),
        ),
        xaxis=dict(
            title='',
            tickfont=dict(size=11),
        ),
        plot_bgcolor='white',
        height=550,
        width=800,
        margin=dict(l=180, r=80, t=100, b=60),
        font=dict(family='Arial'),
    )

    fig.write_image(output_path, scale=2)
    logger.info(f"  Saved: {output_path}")
    return output_path


import logging

logger = logging.getLogger(__name__)

def plot_combined_panel(df, output_path, hub_csv_path):
    """组合面板图 - 完整展示"""
    from plotly.subplots import make_subplots
    compounds = _active_compounds(df)
    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.7, 0.3],
        horizontal_spacing=0.08,
        subplot_titles=['<b>Pathway Enrichment</b>', '<b>Top Hub Genes</b>']
    )

    # 左侧: 堆叠条形图
    for compound in compounds:
        fig.add_trace(go.Bar(
            name=compound,
            y=df['pathway_short'],
            x=df[compound],
            orientation='h',
            marker_color=COLORS[compound],
            text=df[compound],
            textposition='inside',
            textfont=dict(color='white', size=10, family='Arial Black'),
            showlegend=True,
            legendgroup='compounds',
        ), row=1, col=1)

    # 右侧: Hub基因列表
    hub_df = pd.read_csv(hub_csv_path)
    top_genes = hub_df.head(15)

    # 创建表格数据
    fig.add_trace(go.Table(
        header=dict(
            values=['<b>Rank</b>', '<b>Gene</b>', '<b>Degree</b>'],
            fill_color='#1A237E',
            font=dict(color='white', size=12, family='Arial Black'),
            align='center',
            height=28,
        ),
        cells=dict(
            values=[
                top_genes['Rank'].tolist(),
                top_genes['Gene'].tolist(),
                top_genes['Degree'].tolist(),
            ],
            fill_color=[
                ['#F5F5F5'] * len(top_genes),
                ['#E3F2FD'] * len(top_genes),
                ['#BBDEFB'] * len(top_genes),
            ],
            font=dict(size=11, family='Arial'),
            align='center',
            height=25,
        ),
        name='Hub Genes',
    ), row=1, col=2)

    # 更新布局
    fig.update_layout(
        title=dict(
            text='<b>Strawberry Leaf: Compound-Target-Pathway Network Analysis</b><br>' +
                 '<span style="font-size:12px; color:gray">Network Pharmacology of Antioxidant Mechanisms</span>',
            x=0.5,
            font=dict(size=16)
        ),
        barmode='stack',
        yaxis1=dict(
            autorange='reversed',
            tickfont=dict(size=10),
            title='',
        ),
        xaxis1=dict(
            title='Number of Target Genes',
            tickfont=dict(size=10),
        ),
        legend=dict(
            title='<b>Compounds</b>',
            orientation='h',
            yanchor='bottom',
            y=1.06,
            xanchor='center',
            x=0.35,
            font=dict(size=10),
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=650,
        width=1100,
        margin=dict(l=80, r=40, t=100, b=50),
        font=dict(family='Arial'),
        showlegend=True,
    )

    fig.write_image(output_path, scale=2)
    logger.info(f"  Saved: {output_path}")
    return output_path


if __name__ == "__main__":
    # Script is in src/visualization/, outputs is in src/
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'outputs' / 'shared'
    output_dir = base_dir / 'outputs' / 'shared' / 'figures'

    logger.info("\n" + "=" * 60)
    logger.info("  Plotting: Compound-Target-Pathway (Plotly)")
    logger.info("=" * 60)

    logger.info("  Loading data...")
    swiss_data, hub_df, kegg_df = load_data(data_dir)

    logger.info("  Building relationships...")
    df = build_relationship_df(swiss_data, hub_df, kegg_df, top_n=12)

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("\n  Generating plots...")
    plot_stacked_bar(df, str(output_dir / 'compound_target_pathway_plotly_stacked.png'))
    plot_grouped_dot(df, str(output_dir / 'compound_target_pathway_plotly_dot.png'))
    plot_heatmap_style(df, str(output_dir / 'compound_target_pathway_plotly_heatmap.png'))
    plot_combined_panel(df, str(output_dir / 'compound_target_pathway_plotly_panel.png'),
                        data_dir / 'ppi_network' / 'hub_nodes.csv')

    logger.info("\n[OK] All Plotly plots generated!")
