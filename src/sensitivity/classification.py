"""
GO BP/KEGG 功能富集条目分类体系（论文方法学核心）— GO 层级追溯 + KEGG 官方分类 + 关键词版

兼容说明
--------
历史函数名和字段名仍保留 pathway，以兼容既有 CSV/JSON 与绘图脚本。
实际分类对象是 GO BP/KEGG functional enrichment terms；其中 GO BP 是
biological process term，不是 pathway。只有 KEGG 条目可严格称为 pathway。

分类逻辑（优先级从高到低）：
1. direct_ox: GO 层级追溯（种子节点 + 祖先链匹配）
2. cross_talk / signaling_backbone: KEGG BRITE 官方分类自动映射
3. cross_talk: 人工精选的已知与氧化应激有功能交叉的 KEGG 通路/机制条目
4. signaling_backbone: 核心信号转导 + 氧化应激相关关键词（正则单词边界匹配）
5. other: 以上都不是

升级记录
--------
2026-04-26: direct_ox 层从"扁平 GO ID 查表"升级为"GO DAG 向上追溯祖先"。
2026-04-26: 新增 KEGG BRITE 官方分类层。基于 KEGG 官方 CLASS（Human Diseases /
           Signal transduction / Immune system 等）自动判定 KEGG 条目类型，解决
           "Pathways in cancer" 等大量 KEGG 条目被误归为 other 的问题。
2026-04-26: 大幅扩展 BACKBONE_KEYWORDS（从 12 个扩展到 25+ 个），覆盖 cancer、
           apoptosis、HIF-1、inflammation 等氧化应激相关语义域。
2026-04-26 (P0): 修复 " Hedgehog" 前导空格 bug；补充 "MAP kinase", "phosphorylation",
           "kinase cascade" 以覆盖高频 GO_BP 信号调控条目；精简 CROSS_TALK_PATHWAYS，
           删除已由 KEGG CLASS 映射覆盖的具体癌种（Bladder/Breast/Gastric/Prostate/
           Melanoma），避免人工列表与官方分类重复覆盖。
"""

import os
import sys
from io import StringIO
from pathlib import Path
import re
from typing import Dict, List, Optional

from goatools.obo_parser import GODag

# ---------------------------------------------------------------------------
# GO 氧化应激种子集合（来自 oxidative_stress_go_terms.txt 的扁平白名单）
# ---------------------------------------------------------------------------
_GO_OX_FILE = Path(__file__).resolve().parent.parent / "data" / "oxidative_stress_go_terms.txt"
GO_OX_IDS: set = set()
if _GO_OX_FILE.exists():
    with open(_GO_OX_FILE, "r", encoding="utf-8") as f:
        for line in f:
            go_id = line.strip().split("\t")[0]
            GO_OX_IDS.add(go_id)

# ---------------------------------------------------------------------------
# GO DAG 层级追溯
# ---------------------------------------------------------------------------
_OBO_FILE = Path(__file__).resolve().parent.parent / "data" / "go-basic.obo"
_godag: Optional[GODag] = None


def _load_godag() -> Optional[GODag]:
    """延迟加载 GO DAG（第一次调用时加载，后续复用）。"""
    global _godag
    if _godag is None and _OBO_FILE.exists():
        # 捕获 goatools 的加载输出，避免污染 stdout
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            _godag = GODag(str(_OBO_FILE))
        finally:
            sys.stdout = old_stdout
    return _godag


def _is_ox_related_by_hierarchy(go_id: str) -> bool:
    """
    通过 GO DAG 向上追溯祖先，判断该 GO term 是否与氧化应激相关。

    逻辑
    ----
    1. 快速路径：GO ID 本身在种子集合中（精确匹配）。
    2. 层级路径：GO ID 的任一祖先（is_a / part_of）在种子集合中。

    示例
    ----
    GO:0034599 (cellular response to oxidative stress)
        └── 祖先包含 GO:0006979 (response to oxidative stress) → ✅ 命中

    参数
    ----
    go_id : str
        形如 "GO:0006979" 的 GO 标识符。

    返回
    ----
    bool
    """
    # 快速路径：精确匹配
    if go_id in GO_OX_IDS:
        return True

    godag = _load_godag()
    if godag is None or go_id not in godag:
        return False

    # 向上追溯所有祖先（含间接祖先）
    ancestors = godag[go_id].get_all_parents()
    return bool(ancestors & GO_OX_IDS)


# ---------------------------------------------------------------------------
# KEGG 官方分类映射（自动加载）
# ---------------------------------------------------------------------------
_KEGG_CLASS_FILE = Path(__file__).resolve().parent.parent / "data" / "kegg_pathway_classes.json"
KEGG_PATHWAY_CLASSES: Dict[str, str] = {}
if _KEGG_CLASS_FILE.exists():
    import json
    with open(_KEGG_CLASS_FILE, "r", encoding="utf-8") as f:
        KEGG_PATHWAY_CLASSES = json.load(f)

# KEGG CLASS → 分类标签映射规则
# 基于 KEGG BRITE 官方层级分类，将 KEGG 条目自动映射为 cross_talk / signaling_backbone / other
KEGG_CLASS_CATEGORY_MAP = {
    # 疾病大类（含传染病）→ cross_talk
    "Human Diseases; Cancer:": "cross_talk",
    "Human Diseases; Cardiovascular disease": "cross_talk",
    "Human Diseases; Neurodegenerative disease": "cross_talk",
    "Human Diseases; Endocrine and metabolic disease": "cross_talk",
    "Human Diseases; Infectious disease:": "cross_talk",
    # 信号转导 → signaling_backbone
    "Environmental Information Processing; Signal transduction": "signaling_backbone",
    # 免疫系统 → cross_talk
    "Organismal Systems; Immune system": "cross_talk",
    # 内分泌系统 → cross_talk（激素信号与氧化应激密切相关）
    "Organismal Systems; Endocrine system": "cross_talk",
    # 神经系统 → cross_talk（神经退行与氧化应激相关）
    "Organismal Systems; Nervous system": "cross_talk",
    # 排泄系统 → other
    "Organismal Systems; Excretory system": "other",
    # 代谢 → other（除非有明确关联）
    "Metabolism;": "other",
    # 细胞过程 → cross_talk（细胞黏附、通讯与氧化应激相关）
    "Cellular Processes; Cellular community": "cross_talk",
}


def _classify_by_kegg_class(term: str, source: Optional[str] = None) -> Optional[str]:
    """
    基于 KEGG 官方 BRITE CLASS 自动分类。

    仅对 KEGG 来源条目生效，GO 来源返回 None（由 GO 层级和关键词处理）。
    """
    if source != "KEGG" or not KEGG_PATHWAY_CLASSES:
        return None

    class_str = KEGG_PATHWAY_CLASSES.get(term)
    if not class_str:
        return None

    # 前缀匹配 CLASS 映射规则
    for prefix, category in KEGG_CLASS_CATEGORY_MAP.items():
        if class_str.startswith(prefix):
            return category
    return None


# ---------------------------------------------------------------------------
# 交叉调控条目（已知与氧化应激存在功能交叉）
# ---------------------------------------------------------------------------
CROSS_TALK_PATHWAYS: set = {
    # 精确已知与氧化应激存在功能交叉的 KEGG 通路/机制条目
    # （KEGG BRITE 官方分类未覆盖或需要明确确认的边缘 case）
    "AGE-RAGE signaling pathway in diabetic complications",
    "IL-17 signaling pathway",
    "Ferroptosis",
    "Glutathione metabolism",
    "Peroxisome",
    "Fluid shear stress and atherosclerosis",
    "NOD-like receptor signaling pathway",
    "TNF signaling pathway",
    "Toll-like receptor signaling pathway",
    "NF-kappa B signaling pathway",
    # 泛癌通路（/overview / specific mechanisms，而非具体癌种）
    # 具体癌种（Bladder/Breast/Gastric/Prostate/Melanoma）已由 KEGG CLASS 映射覆盖
    "Pathways in cancer",
    "Proteoglycans in cancer",
    "Chemical carcinogenesis",
    "MicroRNAs in cancer",
    "Central carbon metabolism in cancer",
    # 神经退行 / 代谢 / 心血管（KEGG CLASS 已覆盖，此处作为冗余兜底）
    "Alzheimer disease",
    "Hepatitis B",
    "Lipid and atherosclerosis",
}

# ---------------------------------------------------------------------------
# 信号网络主干关键词（核心信号转导和相关 GO BP 过程）
# 依据：KEGG BRITE "Signal Transduction" 类别 (hsa04000–hsa04999)
# 匹配方式：正则单词边界 \bkeyword\b
# ---------------------------------------------------------------------------
BACKBONE_KEYWORDS: List[str] = [
    # 核心信号转导（原有）
    "MAPK",
    "PI3K-Akt",
    "Rap1",
    "Ras signaling",
    "receptor tyrosine kinase",
    "protein tyrosine kinase",
    "cytokine-mediated",
    "calcium signaling",
    "cAMP signaling",
    "chemokine signaling",
    "JAK-STAT",
    "EGFR tyrosine kinase inhibitor resistance",
    # 疾病相关（新增）
    "cancer",
    "carcinogenesis",
    "melanoma",
    "atherosclerosis",
    "Alzheimer",
    "Hepatitis",
    "diabetic complications",
    "cardiovascular",
    # 氧化应激核心生物学过程（新增）
    "apoptosis",
    "apoptotic",
    "autophagy",
    "senescence",
    "inflammation",
    "inflammatory",
    "immune",
    "cytokine",
    "hypoxia",
    # 关键信号通路（新增）
    "HIF-1",
    "Notch",
    "focal adhesion",
    "gap junction",
    "Wnt",
    "TGF-beta",
    "TGF-β",
    "Hedgehog",
    # 补充遗漏的上位信号术语（P0 fix）
    "MAP kinase",        # 匹配 GO_BP 中 "MAP kinase" 写法（与 MAPK 等价）
    "phosphorylation",   # 匹配磷酸化调控过程
    "kinase cascade",    # 匹配信号级联通路
]


def classify_pathway(term: str, source: Optional[str] = None) -> str:
    """GO BP/KEGG 功能富集条目分类（GO 层级追溯版）。

    函数名保留 pathway 是为了兼容历史调用；GO BP 条目不是通路。

    参数
    ----------
    term : str
        富集条目名称（如 "MAPK signaling pathway" 或
        "positive regulation of phosphorylation (GO:0042327)"）。
    source : str, optional
        来源数据库（"KEGG", "GO_BP", "GO_MF" 等）。
        如果为 "GO_MF"，不会进入 signaling_backbone 分类。

    返回
    -------
    str
        "direct_ox" | "cross_talk" | "signaling_backbone" | "other"
    """
    term_lower = term.lower()

    # 1. 直接氧化应激响应（GO 层级追溯）
    match = re.search(r"\(GO:(\d+)\)", term)
    if match:
        go_id = f"GO:{match.group(1)}"
        if _is_ox_related_by_hierarchy(go_id):
            return "direct_ox"

    # 2. KEGG 官方分类自动映射（基于 BRITE CLASS）
    kegg_cat = _classify_by_kegg_class(term, source)
    if kegg_cat is not None:
        return kegg_cat

    # 3. 交叉调控条目（固定列表精确匹配）
    if term in CROSS_TALK_PATHWAYS:
        return "cross_talk"

    # 4. 信号网络主干条目（单词边界匹配）
    #    注意：GO_MF（分子功能）条目不进入此类
    if source == "GO_MF":
        return "other"

    for kw in BACKBONE_KEYWORDS:
        if re.search(r"\b" + re.escape(kw.lower()) + r"\b", term_lower):
            return "signaling_backbone"

    return "other"


def classify_pathways_batch(
    terms: List[str],
    sources: Optional[List[str]] = None,
) -> List[str]:
    """批量分类。"""
    if sources is None:
        sources = [None] * len(terms)
    return [classify_pathway(t, s) for t, s in zip(terms, sources)]
