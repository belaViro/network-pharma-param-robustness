"""
Stage 2 Worker - Swiss 靶点预测
用于 ProcessPoolExecutor，每个进程独立运行
"""

from typing import Dict
from src.connectors.swiss_target import SwissTargetPredictor


def predict_target(smiles: str, download_path: str, organism: str = "Homo_sapiens", molecule_name: str = None, proxy: str = None) -> Dict:
    """单进程 Swiss 靶点预测 worker

    Args:
        smiles: 分子的 SMILES 字符串
        download_path: CSV 文件保存目录
        organism: 物种（默认 Homo_sapiens）
        molecule_name: 分子名称，用于命名CSV文件
        proxy: 代理服务器地址，如 "http://127.0.0.1:7890"

    Returns:
        包含 smiles、name、status、csv 字段的字典
    """
    predictor = SwissTargetPredictor(download_path=download_path)

    try:
        csv_path = predictor.run(smiles=smiles, organism=organism, molecule_name=molecule_name, proxy=proxy)
        if csv_path:
            return {"smiles": smiles, "name": molecule_name or smiles, "status": "成功", "csv": csv_path}
        else:
            return {"smiles": smiles, "name": molecule_name or smiles, "status": "失败", "csv": ""}
    except Exception as e:
        return {"smiles": smiles, "name": molecule_name or smiles, "status": f"错误: {e}", "csv": ""}
    finally:
        predictor.close()
