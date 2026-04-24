"""
STRING Database API Connector
用于获取蛋白质相互作用网络数据

STRING API 文档: https://string-db.org/help/api/
"""

import requests
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
import time
from urllib.parse import quote

# STRING API 端点
STRING_API_URL = "https://string-db.org/api/json/network"


class StringConnector:
    """STRING Database API 连接器"""

    def __init__(
        self,
        species: int = 9606,  # Homo sapiens
        output_dir: Optional[str] = None,
    ):
        """
        Args:
            species: NCBI Taxonomy ID (默认 9606 = Human)
            output_dir: 输出目录
        """
        self.species = species
        if output_dir is None:
            output_dir = Path(__file__).parent.parent / "outputs" / "shared" / "ppi_network"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def get_ppi_network(
        self,
        gene_names: List[str],
        min_confidence: float = 0.4,
    ) -> Dict:
        """获取 PPI 网络

        Args:
            gene_names: 基因名列表
            min_confidence: 最低置信度 (0-1)

        Returns:
            包含网络数据的字典
        """
        logger.info(f"\n{'='*50}")
        logger.info(f"  获取 PPI 网络")
        logger.info(f"{'='*50}")
        logger.info(f"  基因数: {len(gene_names)}")
        logger.info(f"  最低置信度: {min_confidence}")

        # 构建基因列表，用 %0A 分隔
        identifiers = quote("\n".join(gene_names))
        params = {
            "species": self.species,
            "identifiers": identifiers,
            "format": "json",
        }

        try:
            response = requests.get(
                STRING_API_URL,
                params=params,
                timeout=300,  # STRING API 较慢
            )
            response.raise_for_status()
            data = response.json()

            if not isinstance(data, list):
                logger.info("  API 返回数据格式错误")
                return {"status": "错误", "interactions": [], "nodes": []}

            logger.info(f"  获取到 {len(data)} 条相互作用")

            # 构建相互作用 DataFrame
            interactions = []
            nodes_set = set()

            for item in data:
                score = float(item.get("score", 0))
                if score < min_confidence:
                    continue

                p1 = item.get("preferredName_A", "")
                p2 = item.get("preferredName_B", "")

                if p1 and p2:
                    interactions.append({
                        "protein1": p1,
                        "protein2": p2,
                        "combined_score": score,
                    })
                    nodes_set.add(p1)
                    nodes_set.add(p2)

            df = pd.DataFrame(interactions)
            logger.info(f"  置信度过滤后: {len(df)} 条相互作用")
            logger.info(f"  网络节点数: {len(nodes_set)}")

            return {
                "status": "成功",
                "data": data,
                "interactions": df,
                "nodes": list(nodes_set),
            }

        except Exception as e:
            logger.info(f"  获取 PPI 网络错误: {e}")
            return {"status": "错误", "error": str(e), "interactions": [], "nodes": []}

    def build_ppi_network(
        self,
        gene_names: List[str],
        min_score: float = 0.4,
        output_prefix: str = "ppi_network",
    ) -> Dict:
        """构建并保存 PPI 网络

        Args:
            gene_names: 基因名列表
            min_score: 最低置信度分数
            output_prefix: 输出文件前缀

        Returns:
            包含网络统计信息的字典
        """
        result = self.get_ppi_network(gene_names, min_confidence=min_score)

        if result["status"] != "成功":
            return result

        interactions_df = result["interactions"]
        nodes_list = result["nodes"]

        if len(interactions_df) == 0:
            logger.info("  未获取到相互作用")
            return {
                "status": "无数据",
                "genes": len(gene_names),
                "interactions": 0,
            }

        # 保存相互作用
        output_path = self.output_dir / f"{output_prefix}_interactions.csv"
        interactions_df.to_csv(output_path, index=False, encoding="utf-8-sig")
        logger.info(f"  已保存相互作用: {output_path}")

        # 保存节点列表
        nodes_df = pd.DataFrame({"Gene": sorted(nodes_list)})
        nodes_df.to_csv(self.output_dir / f"{output_prefix}_nodes.csv", index=False)

        # 计算节点度数
        from collections import Counter
        all_genes = list(interactions_df["protein1"]) + list(interactions_df["protein2"])
        degree_counts = Counter(all_genes)
        degree_df = pd.DataFrame([
            {"Gene": g, "Degree": d} for g, d in degree_counts.items()
        ]).sort_values("Degree", ascending=False)
        degree_df.to_csv(self.output_dir / f"{output_prefix}_degrees.csv", index=False)

        # 统计
        avg_degree = sum(degree_counts.values()) / max(len(nodes_list), 1)

        logger.info(f"\n{'='*50}")
        logger.info(f"  PPI 网络统计")
        logger.info(f"{'='*50}")
        logger.info(f"  输入基因数: {len(gene_names)}")
        logger.info(f"  网络节点数: {len(nodes_list)}")
        logger.info(f"  相互作用数: {len(interactions_df)}")
        logger.info(f"  平均度数: {avg_degree:.2f}")

        return {
            "status": "成功",
            "genes": len(gene_names),
            "network_nodes": len(nodes_list),
            "interactions": len(interactions_df),
            "avg_degree": avg_degree,
            "output_path": str(output_path),
            "interactions_df": interactions_df,
            "degree_df": degree_df,
        }


def build_ppi_from_targets(
    targets_file: str,
    output_dir: str = None,
    min_score: float = 0.4,
) -> Dict:
    """从靶点文件构建 PPI 网络

    Args:
        targets_file: 靶点文件路径
        output_dir: 输出目录
        min_score: 最低置信度分数

    Returns:
        网络统计信息
    """
    df = pd.read_csv(targets_file)
    gene_names = df["Target_Gene"].tolist()

    connector = StringConnector(output_dir=output_dir)
    return connector.build_ppi_network(gene_names, min_score=min_score)


import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="STRING PPI 网络构建")
    parser.add_argument("--targets", type=str, default=None, help="靶点文件路径")
    parser.add_argument("--min-score", type=float, default=0.4, help="最低置信度分数 (默认 0.4)")
    args = parser.parse_args()

    if args.targets is None:
        targets_file = Path(__file__).parent.parent / "outputs" / "shared" / "intersected_targets.csv"
    else:
        targets_file = Path(args.targets)

    logger.info(f"\n从 {targets_file} 构建 PPI 网络...\n")
    result = build_ppi_from_targets(targets_file, min_score=args.min_score)

    if result["status"] == "成功":
        logger.info(f"\n完成! 相互作用数: {result['interactions']}")