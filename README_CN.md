# 网络药理学参数与成分双扰动鲁棒性框架

[English](./README.md) | **中文**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> **网络药理学靶点优先级识别的双扰动框架：** 参数扰动用于改良固定阈值筛选流程，成分留一扰动用于解析多成分体系中的机制依赖性，并以草莓叶抗氧化活性成分为应用案例。（DOI: 10.7501/j.issn.0253-2670.2026.22.017）

## 概述

本仓库包含一项网络药理学双扰动鲁棒性分析的可复现包。项目同时关注两类不确定性：**参数敏感性**，即人为阈值变化可能改变靶点和 GO BP/KEGG 功能富集结果；**成分输入依赖性**，即候选机制可能由某一个或少数输入成分驱动，而不一定代表整体多成分体系。

### 核心贡献

- **参数扰动作为方法改良**：将单一固定阈值下的 TopK 靶点列表，改造为参数空间中的复现频率和网络位置优先级排序。
- **成分扰动作为机制分析**：通过逐化合物留一运行，识别不依赖单一输入成分的靶点和 GO BP/KEGG 功能富集条目。
- **112 组参数组合**：穷举 SwissTargetPrediction probability 阈值（7）× STRING 置信度阈值（4）× Hub 基因截断（4）的全因子网格。
- **2,016 组成分扰动运行**：18 个结构规范化后的草莓叶化合物 × 112 组参数组合。
- **可审计证据链**：PR-QI 靶点/功能条目优先级、Top30 消融、PubMed 共现审计、degree-matched permutation 和分子对接 seed 稳定性输出。

### 关键发现

| 指标 | 数值 |
|------|------|
| 参数组合总数 | **112**（7 × 4 × 4） |
| 化合物 LOO 组合总数 | **2,016**（18 × 112） |
| SwissTarget 并集靶点 | 684 |
| 正概率疾病靶点 | 286 |
| 主报告阈值靶点（Probability > 0.05） | 219 |
| PR-QI 参数稳健候选靶点 | 50 |
| PR-QI GO BP/KEGG 功能条目 | 70 |
| **双稳健靶点** | **29**（Q1 保留率 58%） |
| **Tier 1 靶点** | **STAT3**（参数频率 100%，LOO 频率 100%） |
| 双稳健 GO BP/KEGG 功能条目 | 58 |
| 最高化合物影响度 | Cardiac glycoside（得分 0.306） |

---

## 快速开始

### 1. 克隆仓库

```bash
git clone https://github.com/belaViro/network-pharma-param-robustness.git
cd network-pharma-param-robustness
```

### 2. 安装依赖

```bash
# 方式 A：pip
pip install -r requirements.txt

# 方式 B：可编辑安装（推荐）
pip install -e .

# 可选：重新获取 PubChem/SwissTargetPrediction 原始数据
pip install -e ".[acquisition]"
```

需要 Python ≥ 3.10。

### 3. 手动运行主分析

```bash
# 传统固定阈值基线
python scripts/run_traditional_pipeline.py

# 参数网格敏感性分析（112 组组合）
python scripts/analysis/run_sensitivity_analysis.py

# 化合物留一法敏感性分析（2,016 组组合）
python scripts/analysis/run_compound_sensitivity_analysis.py

# 外部证据审计与 Top30 消融分析
python scripts/analysis/robustness_validation_bundle.py

# 论文图表重生成
python scripts/analysis/generate_paper_figures.py

# 论文补充表导出
python scripts/analysis/export_supplementary_tables.py
```

### 4. 一键复现

```bash
# 从仓库内置输入文件复现完整论文工作流（默认）
python reproduce.py

# 通过 Bash 包装脚本运行同一流程
bash reproduce.sh

# Windows PowerShell 包装脚本
.\reproduce.ps1

# 只快速复现参数鲁棒性核心流程
python reproduce.py --core
```

默认复现入口会从仓库内置输入文件开始，依次运行论文级结果工作流：传统固定阈值基线 → 112 组参数网格 → 稳健 vs. 参数化对比 → GO BP/KEGG 功能条目象限增强 → 化合物 LOO 扰动 → 外部证据与 degree-matched validation bundle → 最终论文图件 → 补充表。

常用模式：

| 命令 | 用途 |
|------|------|
| `python reproduce.py` | 从内置输入和缓存复现完整论文工作流 |
| `python reproduce.py --core` | 只快速复现参数鲁棒性核心流程 |
| `python reproduce.py --skip-install` | 复用当前环境，不重新安装本包 |
| `python reproduce.py --pubmed-mode query` | 通过 NCBI E-utilities 重新获取 PubMed 共现计数 |
| `python reproduce.py --full --include-acquisition` | 同时重新运行 PubChem/SwissTargetPrediction 获取阶段 |

PubMed 验证默认使用缓存计数，不主动联网查询。重新查询 PubMed 或重新获取原始预测数据都需要网络，并且耗时会明显增加。

---

## 仓库结构

```
.
├── scripts/
│   ├── analysis/
│   │   ├── run_sensitivity_analysis.py           # 主程序：112 组参数网格
│   │   ├── run_compound_sensitivity_analysis.py  # 主程序：18 化合物 LOO
│   │   ├── compare_robust_vs_parametric.py       # 稳健 vs. 参数化对比
│   │   ├── augment_pathway_robustness.py         # GO BP/KEGG 鲁棒性象限增强
│   │   ├── robustness_validation_bundle.py        # Top30 消融与外部证据审计
│   │   ├── generate_paper_figures.py             # 论文图表生成
│   │   ├── export_supplementary_tables.py        # 论文补充表导出
│   │   ├── validate_with_pubmed.py               # 独立 PubMed 诊断脚本
│   │   └── ...                                   # 其他诊断脚本
│   ├── run_pipeline.py                           # 可选 PubChem/SwissTarget 获取
│   ├── run_full_pipeline.py                      # 旧版分阶段获取与基线流水线
│   └── run_traditional_pipeline.py               # 基线方法对比
│
├── src/
│   ├── config.py                                 # 全局配置中心
│   ├── connectors/                               # API 封装
│   │   ├── swiss_target.py                       # SwissTargetPrediction
│   │   ├── string_connector.py                   # STRING PPI
│   │   ├── pubchem.py                            # PubChem
│   │   └── pubmed_connector.py                   # PubMed 共现审计
│   ├── enrichment/
│   │   ├── enrichment_analysis.py                # GO/KEGG 富集（gseapy）
│   │   └── plotting.py                           # 富集可视化
│   ├── visualization/                            # 绘图工具
│   ├── data/                                     # 参考数据
│   │   ├── oxidative_stress_go_terms.txt         # 氧化应激相关 GO 条目
│   │   └── kegg_pathway_classes.json             # KEGG 通路分类
│   └── outputs/                                  # 全部分析输出
│       ├── paper_figures/                        # 论文成图
│       ├── sensitivity/                          # 参数敏感性结果
│       │   ├── target_robustness.csv
│       │   ├── pathway_robustness.csv
│       │   ├── param_grid_results.csv
│       │   └── compound_sensitivity/             # LOO 化合物结果
│       ├── shared/                               # 跨阶段共享数据
│       └── traditional/                          # 基线方法输出
│
├── docs/                                         # 项目文档
├── 分子对接/                                      # 分子对接拼图与结果图
├── requirements.txt
├── pyproject.toml
├── reproduce.py                                  # 跨平台复现总控脚本
├── reproduce.ps1                                 # Windows PowerShell 包装脚本
└── reproduce.sh                                  # Bash 包装脚本
```

---

## 输出文件

### 核心结果（参数扰动）

| 文件 | 说明 |
|------|------|
| `src/outputs/sensitivity/target_robustness.csv` | 112 组组合下各靶点的频率与度统计 |
| `src/outputs/sensitivity/pathway_robustness.csv` | GO BP/KEGG 功能富集条目鲁棒性指标 |
| `src/outputs/sensitivity/param_grid_results.csv` | 每组组合的富集汇总 |
| `src/outputs/sensitivity/pathway_rank_stability.csv` | 功能富集条目排名在参数空间中的一致性 |

### 核心结果（成分扰动）

| 文件 | 说明 |
|------|------|
| `src/outputs/sensitivity/compound_sensitivity/compound_influence_scores.csv` | 各化合物的影响度评分（0–1） |
| `src/outputs/sensitivity/compound_sensitivity/dual_robust_targets.csv` | 29 个双稳健靶点及其 LOO 保留率 |
| `src/outputs/sensitivity/compound_sensitivity/dual_robust_pathways.csv` | 58 个双稳健 GO BP/KEGG 功能条目及其 LOO 保留率 |
| `src/outputs/sensitivity/compound_sensitivity/leave_one_out_pathway_robustness.csv` | 各化合物剔除后的功能富集条目稳健性 |

### 外部证据与论文图表

兼容说明：历史文件名和部分 CSV 字段仍保留 `pathway_*`，但本仓库中这些字段实际表示 **KEGG 通路 + GO BP 生物过程条目** 的混合功能富集条目；GO BP 条目不是通路。

| 文件 | 说明 |
|------|------|
| `src/outputs/sensitivity/robustness_validation/` | Top30 消融、PubMed 共现审计和 degree-matched permutation 输出 |
| `src/outputs/paper_figures/` | 论文成图，包括工作流、参数热图、成分 LOO 热图、PubMed 阈值审计和分子对接图 |
| `src/outputs/paper_figures/figure-catalog.md` | 图件用途、数据来源与解释提示 |
| `src/outputs/reproduction_manifest.json` | `reproduce.py` 生成的运行清单；已被 Git 忽略 |

### 检查点

- `src/outputs/sensitivity/compound_sensitivity/checkpoints/combo_results/` — 2,016 个独立 JSON 文件（每个化合物 × 组合一份），支持全量审计与断点续跑。

### 论文补充表

| 文件 | 说明 |
|------|------|
| `src/outputs/supplementary/supplementary_tables_v1.xlsx` | Table S1-S6 汇总工作簿 |
| `src/outputs/supplementary/csv/` | 每张补充表的 CSV 副本 |

生成命令：

```bash
python scripts/analysis/export_supplementary_tables.py
```

---

## 参数网格

| 参数 | 取值 | 说明 |
|-----------|--------|-------------|
| `Prob`（SwissTarget） | 0.0, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20 | 靶点纳入概率阈值 |
| `STRING`（PPI） | 0.4, 0.5, 0.7, 0.9 | STRING 综合得分阈值 |
| `TopK`（Hub） | 10, 20, 30, 50 | 按度排名选取的 Hub 基因数 |
| `Pvalue`（固定） | 0.01 | 富集显著性截断 |

**总组合数：7 × 4 × 4 = 112**

---

## 草莓叶活性成分（n = 18）

LC-MS/MS 共鉴定 19 条原始记录；经结构规范化去重后（Trehalose ≡ alpha,alpha-Trehalose）得到 18 个化合物：

1. Miquelianin（槲皮素-3-O-葡萄糖苷酸）
2. Sucrose（蔗糖）
3. Trehalose（海藻糖）
4. Oleamide（油酰胺）
5. Catechin（儿茶素）
6. Indole-3-acrylic acid（吲哚-3-丙烯酸）
7. DL-Tryptophan（DL-色氨酸）
8. (3β,9ξ)-3-(β-D-Glucopyranosyloxy)-14-hydroxycard-20(22)-enolide *(强心苷类)*
9. Denin
10. 2-Methylglutaric acid（2-甲基戊二酸）
11. Procyanidin B1（原花青素 B1）
12. Scutellarin（野黄芩苷）
13. Choline（胆碱）
14. DL-Malic acid（DL-苹果酸）
15. Citric acid（柠檬酸）
16. A-Ketoglutaric Acid（α-酮戊二酸）
17. L-Allo-Isoleucine（L-别异亮氨酸）
18. Hexadecanamide（十六碳酰胺）

---

## 靶点分类体系

| 层级 | 判定标准 | 数量 | 代表靶点 |
|------|----------|------|----------|
| **Tier 1** | 参数频率 = 100% 且 LOO 频率 = 100% | 1 | **STAT3** |
| **Tier 2** | 参数频率 ≥ 50% 且双稳健 | 13 | EGFR, HSP90AA1, TNF, BCL2, SRC, ... |
| **Tier 3** | 参数频率 < 50% 但双稳健 | 15 | MAPK14, MET, IL2, KIT, CA1, ... |
| **化合物敏感型** | 参数频率较高但剔除关键化合物后消失 | 3 | CASP3, MAPK1, PIK3CA |
| **外围型** | 参数频率低且非双稳健 | 18 | APP, MAPK8, BCL2L1, MTOR, ... |

---

## 依赖

核心分析：
- pandas ≥ 2.0
- numpy ≥ 1.24
- networkx ≥ 3.0
- scipy ≥ 1.10
- gseapy ≥ 1.0（富集分析）
- goatools ≥ 1.3（GO 分析）

可视化：
- matplotlib ≥ 3.6
- seaborn ≥ 0.12

数据采集（可选）：
- playwright ≥ 1.40
- selenium ≥ 4.0

详见 `requirements.txt` 或 `pyproject.toml`。

---

## 数据可用性

本仓库已包含全部源代码、原始预测数据、中间检查点及最终结果。

- **仓库地址**：https://github.com/belaViro/network-pharma-param-robustness
- **可复现包**：仓库克隆本身已包含 `reproduce.py` 所需的内置输入、缓存中间文件、检查点和论文级输出。

使用的外部数据库：
- SwissTargetPrediction（http://www.swisstargetprediction.ch/）
- STRING v12（https://string-db.org/）
- Enrichr（https://maayanlab.cloud/Enrichr/）
- Gene Ontology（http://purl.obolibrary.org/obo/go-basic.obo）

---

## 引用

若使用本代码或框架，请引用：

> *Parameter and compound dual-perturbation robustness framework for network pharmacology target prioritization: strawberry leaf antioxidant compounds as a case study.*

---

## 许可证

[MIT License](./LICENSE)
