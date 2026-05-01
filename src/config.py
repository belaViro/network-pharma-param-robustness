"""
全局配置中心
集中管理路径、API参数、分析阈值等常量
"""

from pathlib import Path

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
PRIMARY_PROBABILITY_THRESHOLD = 0.05
PROBABILITY_THRESHOLD = PRIMARY_PROBABILITY_THRESHOLD
PROBABILITY_THRESHOLDS: list[float] = [0.0, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]

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
SENSITIVITY_FULL_TARGETS_CSV = SENSITIVITY_DIR / "intersected_targets_probability_gt0.csv"
SENSITIVITY_PPI_NETWORK_DIR = SENSITIVITY_DIR / "ppi_network_full"
SENSITIVITY_PPI_EDGES_CSV = SENSITIVITY_PPI_NETWORK_DIR / "string_interactions_interactions.csv"

# ============================================================================
# 日志配置
# ============================================================================

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
LOG_LEVEL = "INFO"

# ============================================================================
# 化合物列表
# ============================================================================

# LC-MS/MS 原始鉴定成分为 19 条记录；Trehalose 与 alpha,alpha-Trehalose
# 为同一化合物，正式分析使用 18 个去重后的 canonical compound。
STRAWBERRY_LEAF_RAW_COMPOUNDS: list[str] = [
    "Miquelianin",
    "Sucrose",
    "Trehalose",
    "Oleamide",
    "Catechin",
    "Indole-3-acrylic acid",
    "DL-Tryptophan",
    "(3beta,9xi)-3-(beta-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide",
    "Denin",
    "2-Methylglutaric acid",
    "Procyanidin B1",
    "Scutellarin",
    "Choline",
    "alpha,alpha-Trehalose",
    "DL-Malic acid",
    "Citric acid",
    "A-Ketoglutaric Acid",
    "L-Allo-Isoleucine",
    "Hexadecanamide",
]

STRAWBERRY_LEAF_COMPOUND_CANONICAL_MAP: dict[str, str] = {
    "alpha,alpha-Trehalose": "Trehalose",
}

STRAWBERRY_LEAF_COMPOUND_ALIASES: dict[str, list[str]] = {
    "Miquelianin": [
        "Quercetin 3-O-glucuronide",
        "Quercetin-3-O-glucuronide",
    ],
    "Trehalose": [
        "alpha,alpha-Trehalose",
        "D-Trehalose",
    ],
    "Indole-3-acrylic acid": [
        "3-Indoleacrylic acid",
        "Indoleacrylic acid",
    ],
    "DL-Tryptophan": [
        "Tryptophan",
    ],
    "(3beta,9xi)-3-(beta-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide": [
        "3-(beta-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide",
    ],
    "A-Ketoglutaric Acid": [
        "alpha-Ketoglutaric acid",
        "2-Oxoglutaric acid",
    ],
    "L-Allo-Isoleucine": [
        "L-Alloisoleucine",
    ],
    "Hexadecanamide": [
        "Palmitamide",
    ],
}

STRAWBERRY_LEAF_COMPOUND_MANUAL_SMILES: dict[str, dict[str, str | int]] = {
    # PubChem name search does not reliably resolve "Denin"; use the curated
    # structure recorded by KNApSAcK/CropMetabolome for C33H56O14 Denin.
    "Denin": {
        "smiles": "CC(O)[C@H]1CC[C@]2(O)[C@@H]3CC[C@H]4C[C@@H](O[C@@H]5O[C@@H](C)[C@@H](O[C@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)[C@@H]5O)CC[C@]4(C)[C@H]3[C@H](O)[C@@H](O)[C@]12C",
        "cid": 101268475,
        "source": "manual_curated",
    },
}


def canonicalize_compound_name(name: str) -> str:
    """Return the deduplicated canonical compound name used in analysis."""
    return STRAWBERRY_LEAF_COMPOUND_CANONICAL_MAP.get(name, name)


def get_strawberry_leaf_molecules(include_duplicates: bool = False) -> list[str]:
    """Return strawberry leaf compound names for pipeline execution.

    Parameters
    ----------
    include_duplicates:
        If True, return all 19 LC-MS/MS raw records. If False, collapse
        duplicate chemical identities and return the 18 canonical compounds
        used for target prediction and robustness analysis.
    """
    if include_duplicates:
        return list(STRAWBERRY_LEAF_RAW_COMPOUNDS)

    molecules: list[str] = []
    seen: set[str] = set()
    for raw_name in STRAWBERRY_LEAF_RAW_COMPOUNDS:
        canonical = canonicalize_compound_name(raw_name)
        if canonical not in seen:
            molecules.append(canonical)
            seen.add(canonical)
    return molecules


def safe_compound_file_stem(name: str) -> str:
    """Return the filename stem used for per-compound SwissTarget CSV files."""
    safe_name = "".join(c for c in name if c.isalnum() or c in (" ", "-", "_")).strip()
    return safe_name.replace(" ", "_")


STRAWBERRY_LEAF_MOLECULES: list[str] = get_strawberry_leaf_molecules()


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
