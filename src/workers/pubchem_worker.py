"""
Stage 1 Worker - PubChem 查询
用于 ThreadPoolExecutor，查询分子 SMILES
"""

import logging
from typing import Dict, Optional
from src.connectors.pubchem import PubChemClient

logger = logging.getLogger(__name__)


def fetch_smiles(molecule_name: str) -> Dict:
    """根据分子名称从 PubChem 获取 SMILES

    Args:
        molecule_name: 分子名称（如 "Aspirin"）

    Returns:
        包含 name、smiles、status 字段的字典
    """
    client = PubChemClient()

    try:
        results = client.get_compound_by_name(molecule_name)
        if results and results[0].get("smiles"):
            smiles = results[0]["smiles"]
            cid = results[0].get("cid")
            logger.info(f"成功获取 {molecule_name} 的 SMILES: {smiles[:30]}... (CID: {cid})")
            return {"name": molecule_name, "smiles": smiles, "cid": cid, "status": "成功"}
        else:
            logger.warning(f"无法获取 {molecule_name} 的 SMILES（无结果）")
            return {"name": molecule_name, "smiles": None, "cid": None, "status": "未找到"}
    except Exception as e:
        logger.warning(f"查询 {molecule_name} 失败: {e}")
        return {"name": molecule_name, "smiles": None, "cid": None, "status": f"错误: {e}"}
