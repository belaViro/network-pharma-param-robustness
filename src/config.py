"""
全局配置中心
集中管理路径、API参数、分析阈值等常量
"""

from pathlib import Path
from typing import Dict, Any

# ============================================================================
# 项目路径
# ============================================================================

# 项目根目录 (code/)
PROJECT_ROOT = Path(__file__).parent.parent

# 源代码目录
SRC_DIR = PROJECT_ROOT / "src"

# 输出目录
OUTPUTS_DIR = SRC_DIR / "outputs"
SHARED_DIR = OUTPUTS_DIR / "shared"
TRADITIONAL_DIR = OUTPUTS_DIR / "traditional"
REPORTS_DIR = OUTPUTS_DIR / "reports"
SENSITIVITY_DIR = OUTPUTS_DIR / "sensitivity"

# Stage 输出子目录
SWISS_TARGETS_DIR = SHARED_DIR / "swiss_targets"
DISEASE_GENES_DIR = SHARED_DIR / "disease_genes"
PPI_NETWORK_DIR = SHARED_DIR / "ppi_network"
SHARED_FIGURES_DIR = SHARED_DIR / "figures"
TRADITIONAL_ENRICHMENT_DIR = TRADITIONAL_DIR / "enrichment"
TRADITIONAL_FIGURES_DIR = TRADITIONAL_DIR / "figures"
TRADITIONAL_EVAL_DIR = TRADITIONAL_DIR / "evaluation"

# 关键文件路径
MERGED_TARGETS_CSV = SHARED_DIR / "merged_targets.csv"
MERGED_DISEASE_GENES_CSV = DISEASE_GENES_DIR / "merged_disease_genes.csv"
INTERSECTED_TARGETS_CSV = SHARED_DIR / "intersected_targets.csv"
PPI_EDGES_CSV = PPI_NETWORK_DIR / "string_interactions_interactions.csv"
PPI_NODES_CSV = PPI_NETWORK_DIR / "network_nodes.csv"
PPI_DEGREES_CSV = PPI_NETWORK_DIR / "node_degrees.csv"
PPI_HUB_NODES_CSV = PPI_NETWORK_DIR / "hub_nodes.csv"

# ============================================================================
# 外部 API 配置
# ============================================================================

# PubChem
PUBCHEM_API_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/JSON"

# SwissTargetPrediction
SWISS_URL = "https://www.swisstargetprediction.ch/"
SWISS_ORGANISM = "Homo_sapiens"
SWISS_TIMEOUT = 15000  # ms
SWISS_MAX_RETRIES = 3

# STRING
STRING_API_URL = "https://string-db.org/api/json/network"
STRING_SPECIES = 9606  # Homo sapiens
STRING_MIN_SCORE = 0.4
STRING_API_TIMEOUT = 300  # seconds

# ============================================================================
# 分析参数
# ============================================================================

# 靶点预测
PROBABILITY_THRESHOLD = 0.0  # SwissTarget Probability > 0

# 富集分析
PVALUE_CUTOFF = 0.01
PVALUE_CUTOFF_RELAXED = 0.05
ENRICHMENT_TOP_N = 20
BACKGROUND_GENES = 20500  # 人类蛋白编码基因总数

# Hub 基因筛选
TOP_K_HUB = 30

# KEGG 抗氧化通路关键词
ANTIOXIDANT_KEYWORDS: list[str] = [
    'oxidative stress',
    'ROS',
    'ferroptosis',
    'glutathione',
    'catalase',
    'peroxidase',
    'superoxide',
    'NF-kappa',
    'MAPK',
    'p53',
    'apoptosis',
    'NADPH oxidase',
    'thioredoxin',
    'lipid peroxidation',
    'autophagy',
    'mitophagy',
    'AGE-RAGE',
    'IL-17',
    'TNF',
    'NOD-like',
]

# 分析参数目录
SENSITIVITY_FIGURES_DIR = SENSITIVITY_DIR / "figures"
SENSITIVITY_EVAL_DIR = SENSITIVITY_DIR / "evaluation"

# ============================================================================
# 日志配置
# ============================================================================

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
LOG_LEVEL = "INFO"

# ============================================================================
# 化合物列表
# ============================================================================

STRAWBERRY_LEAF_MOLECULES: list[str] = [
    'Sucrose',
    'Trehalose',
    'Catechin',
    'Procyanidin B1',
    'Scutellarin',
]


# ============================================================================
# 辅助函数
# ============================================================================

def ensure_output_dirs() -> None:
    """确保所有输出目录存在"""
    dirs = [
        SHARED_DIR,
        SWISS_TARGETS_DIR,
        DISEASE_GENES_DIR,
        PPI_NETWORK_DIR,
        SHARED_FIGURES_DIR,
        TRADITIONAL_DIR,
        TRADITIONAL_ENRICHMENT_DIR,
        TRADITIONAL_FIGURES_DIR,
        TRADITIONAL_EVAL_DIR,
        REPORTS_DIR,
        SENSITIVITY_DIR,
        SENSITIVITY_FIGURES_DIR,
        SENSITIVITY_EVAL_DIR,
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def get_report_filename(prefix: str, suffix: str = ".md") -> Path:
    """生成带时间戳的报告文件名"""
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return REPORTS_DIR / f"{prefix}_{timestamp}{suffix}"
