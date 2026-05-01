"""
Stage 1 Worker - PubChem 查询
用于 ThreadPoolExecutor，查询分子 SMILES
"""

import logging
from typing import Dict, Optional
from src.connectors.pubchem import PubChemClient
from src.config import (
    STRAWBERRY_LEAF_COMPOUND_ALIASES,
    STRAWBERRY_LEAF_COMPOUND_MANUAL_SMILES,
)

logger = logging.getLogger(__name__)


def fetch_smiles(molecule_name: str) -> Dict:
    """根据分子名称从 PubChem 获取 SMILES

    Args:
        molecule_name: 分子名称（如 "Aspirin"）

    Returns:
        包含 name、smiles、status 字段的字典
    """
    manual = STRAWBERRY_LEAF_COMPOUND_MANUAL_SMILES.get(molecule_name)
    if manual:
        smiles = str(manual["smiles"])
        cid = manual.get("cid")
        source = str(manual.get("source", "manual_curated"))
        logger.info(
            f"使用手工校准结构获取 {molecule_name} 的 SMILES: {smiles[:30]}... "
            f"(CID: {cid}, source: {source})"
        )
        return {
            "name": molecule_name,
            "query_name": source,
            "smiles": smiles,
            "cid": cid,
            "status": "成功",
        }

    client = PubChemClient()
    query_names = [molecule_name]
    query_names.extend(STRAWBERRY_LEAF_COMPOUND_ALIASES.get(molecule_name, []))

    try:
        for query_name in query_names:
            results = client.get_compound_by_name(query_name)
            if results and results[0].get("smiles"):
                smiles = results[0]["smiles"]
                cid = results[0].get("cid")
                logger.info(
                    f"成功获取 {molecule_name} 的 SMILES: {smiles[:30]}... "
                    f"(CID: {cid}, query: {query_name})"
                )
                return {
                    "name": molecule_name,
                    "query_name": query_name,
                    "smiles": smiles,
                    "cid": cid,
                    "status": "成功",
                }

        logger.warning(f"无法获取 {molecule_name} 的 SMILES（无结果）")
        return {"name": molecule_name, "query_name": "", "smiles": None, "cid": None, "status": "未找到"}
    except Exception as e:
        logger.warning(f"查询 {molecule_name} 失败: {e}")
        return {"name": molecule_name, "query_name": "", "smiles": None, "cid": None, "status": f"错误: {e}"}
