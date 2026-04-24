"""
Stage 2 - Swiss 靶点预测流水线
"""

import multiprocessing as mp
import time
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict
from src.workers import predict_target
import logging


logger = logging.getLogger(__name__)

def stage2_predict_all(
    smiles_results: List[Dict],
    download_path: str,
    max_workers: int = None,
    proxy: str = None,
) -> List[Dict]:
    """多进程 Swiss 靶点预测

    ⚠️ SwissTargetPrediction 是单表单提交网站，强烈建议 max_workers=1
    并行提交会被 WAF 识别为机器人行为，导致 IP 被封。

    Args:
        smiles_results: Stage 1 的结果列表（只处理 status=="成功" 的）
                      每个元素包含 name 和 smiles 字段
        download_path: CSV 保存目录
        max_workers: 最大进程数（默认 1，强烈不建议超过 2）
        proxy: 代理服务器地址，如 "http://127.0.0.1:7890"

    Returns:
        每个分子的预测结果列表
    """
    valid_results = [r for r in smiles_results if r["status"] == "成功"]
    if not valid_results:
        return []

    if max_workers is None:
        max_workers = 1  # 安全默认值

    smiles_name_list = [(r["smiles"], r.get("name")) for r in valid_results]
    results = []

    # 单进程串行执行，避免触发反爬，复用同一个浏览器实例
    if max_workers == 1:
        from src.connectors.swiss_target import SwissTargetPredictor
        predictor = SwissTargetPredictor(download_path=download_path)
        batch_results = predictor.run_batch(
            smiles_list=smiles_name_list,
            organism="Homo_sapiens",
            proxy=proxy,
            max_retries=3,
        )
        for result in batch_results:
            status_icon = "[+]" if result["status"] == "成功" else "[-]"
            name = result.get("name", result.get("smiles", ""))
            logger.info(f"  [Stage2] {status_icon} {name}: {result['status']}")
        return batch_results

    # 多进程模式（不推荐，保留用于兼容旧配置）
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(predict_target, s, download_path, None, name, proxy): (s, name)
            for s, name in smiles_name_list
        }

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            status_icon = "[+]" if result["status"] == "成功" else "[-]"
            name = result.get("name", result["smiles"])
            logger.info(f"  [Stage2] {status_icon} {name}: {result['status']}")

    return results
