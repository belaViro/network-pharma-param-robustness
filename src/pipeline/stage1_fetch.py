"""
Stage 1 - PubChem 查询流水线
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict
from src.workers import fetch_smiles
import logging


logger = logging.getLogger(__name__)

def stage1_fetch_all(molecule_names: List[str]) -> List[Dict]:
    """多线程获取所有分子的 SMILES

    Args:
        molecule_names: 分子名称列表

    Returns:
        每个分子的查询结果列表
    """
    results = []

    with ThreadPoolExecutor(max_workers=len(molecule_names)) as executor:
        futures = {
            executor.submit(fetch_smiles, name): name
            for name in molecule_names
        }

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            status_icon = "[+]" if result["status"] == "成功" else "[-]"
            logger.info(f"  [Stage1] {status_icon} {result['name']}: {result['status']}")

    return results