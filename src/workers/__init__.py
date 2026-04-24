"""
workers - Stage worker 函数（用于 ThreadPool/ProcessPool）
"""

from src.workers.pubchem_worker import fetch_smiles
from src.workers.swiss_worker import predict_target

__all__ = [
    "fetch_smiles",
    "predict_target",
]
