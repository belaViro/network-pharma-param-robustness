"""
Enrichment Analysis Module
GO and KEGG enrichment analysis for PPI network hub genes
"""

from src.enrichment.enrichment_analysis import run_enrichment_analysis
from src.enrichment.plotting import (
    plot_go_enrichment,
    plot_kegg_enrichment,
    plot_go_grouped_barplot,
    plot_all_enrichment,
)

__all__ = [
    "run_enrichment_analysis",
    "plot_go_enrichment",
    "plot_kegg_enrichment",
    "plot_go_grouped_barplot",
    "plot_all_enrichment",
]
