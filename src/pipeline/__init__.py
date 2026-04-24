"""
pipeline - 流水线编排
"""

from src.pipeline.stage1_fetch import stage1_fetch_all
from src.pipeline.stage2_predict import stage2_predict_all
from src.pipeline.stage3_merge_targets import merge_and_deduplicate_targets
from src.pipeline.stage4_merge_disease_genes import merge_disease_genes
from src.pipeline.stage5_intersection_analysis import intersect_targets_with_disease_genes
from src.pipeline.stage6_ppi_network import build_ppi_network

__all__ = [
    "stage1_fetch_all",
    "stage2_predict_all",
    "merge_and_deduplicate_targets",
    "merge_disease_genes",
    "intersect_targets_with_disease_genes",
    "build_ppi_network",
]
