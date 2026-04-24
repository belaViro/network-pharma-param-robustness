"""
connectors - 外部服务连接器（API 客户端 + 浏览器自动化）
"""

from src.connectors.pubchem import PubChemClient, PubChemConfig, PubChemError
from src.connectors.swiss_target import SwissTargetPredictor, SwissTargetConfig
from src.connectors.string_connector import StringConnector

__all__ = [
    "PubChemClient",
    "PubChemConfig",
    "PubChemError",
    "SwissTargetPredictor",
    "SwissTargetConfig",
    "StringConnector",
]
