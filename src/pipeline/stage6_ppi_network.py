"""
Stage 6 - PPI网络构建
使用STRING Database API构建蛋白质相互作用网络
"""

import pandas as pd
from pathlib import Path
from typing import Dict

from ..connectors.string_connector import StringConnector


def build_ppi_network(
    targets_file: str = None,
    output_dir: str = None,
    min_score: float = None,
) -> Dict:
    """构建PPI网络

    Args:
        targets_file: 靶点文件路径，默认 intersected_targets.csv
        output_dir: 输出目录，默认 src/outputs/shared/ppi_network
        min_score: STRING最低置信度分数 (0-1，默认从配置读取)

    Returns:
        包含网络统计信息的字典
    """
    import logging
    from src.config import STRING_MIN_SCORE, INTERSECTED_TARGETS_CSV, PPI_NETWORK_DIR
    logger = logging.getLogger(__name__)

    if targets_file is None:
        targets_file = INTERSECTED_TARGETS_CSV
    else:
        targets_file = Path(targets_file)

    if output_dir is None:
        output_dir = PPI_NETWORK_DIR
    else:
        output_dir = Path(output_dir)

    if min_score is None:
        min_score = STRING_MIN_SCORE

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Stage 6: PPI网络构建 (min_score={min_score})")

        # 读取靶点
        targets_df = pd.read_csv(targets_file)
        gene_names = targets_df["Target_Gene"].tolist()
        logger.info(f"读取 {len(gene_names)} 个靶点")

        # 构建PPI网络
        connector = StringConnector(output_dir=str(output_dir))
        result = connector.build_ppi_network(
            gene_names=gene_names,
            min_score=min_score,
            output_prefix="string_interactions",
        )

        if result["status"] == "成功":
            interactions_df = result["interactions_df"]

            # 保存节点列表
            nodes = sorted(set(interactions_df["protein1"]) | set(interactions_df["protein2"]))
            nodes_df = pd.DataFrame({"Gene": nodes})
            nodes_df.to_csv(output_dir / "network_nodes.csv", index=False)

            # 统计节点度数
            from collections import Counter
            all_genes = list(interactions_df["protein1"]) + list(interactions_df["protein2"])
            degree_counts = Counter(all_genes)
            degree_df = pd.DataFrame([
                {"Gene": g, "Degree": d} for g, d in degree_counts.items()
            ]).sort_values("Degree", ascending=False)
            degree_df.to_csv(output_dir / "node_degrees.csv", index=False)

            logger.info(f"Stage 6 完成: {len(nodes)} 节点, {len(interactions_df)} 边")

        return result

    except FileNotFoundError as e:
        logger.error(f"[Stage6] 靶点文件未找到: {e}")
        return {"status": "错误", "error": f"靶点文件未找到: {e}"}
    except Exception as e:
        logger.exception(f"[Stage6] PPI网络构建失败: {e}")
        return {"status": "错误", "error": str(e)}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Stage 6: PPI网络构建")
    parser.add_argument(
        "--targets",
        type=str,
        default=None,
        help="靶点文件路径"
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.4,
        help="STRING最低置信度分数 (0-1)"
    )
    args = parser.parse_args()

    result = build_ppi_network(
        targets_file=args.targets,
        min_score=args.min_score,
    )
