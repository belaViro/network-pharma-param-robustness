# 网络药理学参数与成分双扰动鲁棒性框架

(./README.md) | **中文**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> **网络药理学靶点优先级识别的双扰动框架：** 参数扰动用于改良固定阈值筛选流程，成分留一扰动用于解析多成分体系中的机制依赖性，并以草莓叶抗氧化活性成分为应用案例。

## 概述

本仓库包含一项网络药理学双扰动鲁棒性分析的可复现包。项目同时关注两类不确定性：**参数敏感性**，即人为阈值变化可能改变靶点和 GO BP/KEGG 功能富集结果；**成分输入依赖性**，即候选机制可能由某一个或少数输入成分驱动，而不一定代表整体多成分体系。

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

## 方法学创新思路

### 1. 从"点估计"到"分布估计"的范式转变

传统网络药理学依赖单一固定阈值组合（如 Probability > 0.1、STRING > 0.4、TopK = 20），本质上是靶点优先级的一次"点估计"——分析结论完全取决于这些人为选择的阈值。本工作将这一流程重构为参数空间中的"分布估计"：在 112 组参数组合构成的全因子网格中系统采样，以靶点在参数空间中的**复现频率**和**网络拓扑一致性**替代单一的"在/不在"二元判断。这使得研究者能够区分"在大多数合理参数下都稳定的靶点"和"仅在特定阈值下才出现的靶点"。

### 2. 双维度正交扰动框架

本工作提出并实现了一个双维度扰动框架，两个维度相互正交且互补：

| 扰动维度 | 核心问题 | 方法 | 输出指标 |
|----------|----------|------|----------|
| **参数扰动** | 方法选择是否影响结论？ | 全因子网格（7×4×4 = 112 组） | 参数复现频率、度中心性排序 |
| **成分扰动** | 输入选择是否影响结论？ | 逐化合物留一法（18×112 = 2,016 组） | LOO 保留率、化合物影响度评分 |

- **参数维度**系统性地改变 SwissTarget 概率阈值、STRING 置信度阈值和 Hub 基因截断数，考察分析结论对这些人为选择参数的敏感性。
- **成分维度**通过逐化合物留一法，考察多成分体系的分析结论是否由某一个或少数几个输入成分驱动。

该框架的核心价值在于：将传统上被视为"技术细节"的阈值选择和成分纳入问题，提升为方法学的核心议题，使研究者能够明确量化并报告其分析结论的稳健性边界。

### 3. PR-QI 优先级排序体系

提出**参数鲁棒性定量指标**（Parameter Robustness - Quantitative Index, PR-QI），综合靶点在参数空间中的出现频率和网络拓扑中心性，为靶点优先级提供**连续、可比较**的排序依据。相比传统方法中"高于阈值则入选、低于阈值则排除"的离散分类，PR-QI 使研究者能够：

- 根据稳健性排序合理设定靶点筛选截断
- 识别"边界"靶点（borderline targets）并单独分析
- 为后续实验验证提供基于稳健性的优先级指导

### 4. 成分影响度量化

定义并计算**化合物影响度评分**（Compound Influence Score, 0–1），量化每个输入成分对整体分析结论的边际贡献。该指标可识别：

- **关键驱动化合物**：剔除后导致大量靶点或功能条目失稳的化合物（如本案例中的 Cardiac glycoside, score = 0.306）
- **可替代化合物**：剔除后对整体结论影响甚微的化合物

这一量化方法为多成分体系的机制解析提供了新的分析维度——从"所有成分共同作用"的均质化假设，走向"量化每个成分的差异化贡献"。

### 5. 多层次靶点分类体系

基于双维度鲁棒性指标，建立靶点的四类分级体系（Tier 1–3 + 化合物敏感型 + 外围型），超越了传统"显著/不显著"的二分法：

| 层级 | 判定标准 | 生物学意义 |
|------|----------|------------|
| **Tier 1** | 参数频率 = 100% 且 LOO 频率 = 100% | 方法学上最可信的核心靶点 |
| **Tier 2** | 参数频率 ≥ 50% 且双稳健 | 高置信度稳健靶点 |
| **Tier 3** | 参数频率 < 50% 但双稳健 | 参数敏感但仍具成分鲁棒性的靶点 |
| **化合物敏感型** | 参数频率较高但剔除关键化合物后消失 | 依赖特定成分的靶点，需谨慎解读 |
| **外围型** | 参数频率低且非双稳健 | 仅在少数参数组合下出现的靶点 |

该分类体系帮助研究者在报告结论时明确标注每个靶点的证据强度，避免将"仅在特定条件下出现"的靶点与"在所有合理条件下均稳定"的靶点混为一谈。

### 6. 外部证据多模态审计

超越计算层面的鲁棒性验证，引入三种独立的外部证据线进行交叉验证：

- **PubMed 共现审计**：通过 PubMed 文献数据库验证预测靶点与疾病在已有文献中的共现情况，区分"计算预测新颖靶点"和"已有文献支持的已知靶点"。
- **Degree-matched permutation 检验**：通过度匹配随机置换，排除网络拓扑中心性偏差对靶点优先级排序的潜在混淆效应。
- **分子对接 seed 稳定性**：对关键靶点-化合物对进行多 seed 分子对接（每对 3 个 seed），从结构生物学角度验证结合模式的 seed 稳定性，区分"稳定结合模式"和"seed 敏感的偶然对接结果"。

### 7. 敏感性分析的角色升级：从"验证工具"到"发现工具"

传统上，敏感性分析被视为分析流程末尾的"稳健性检查"——确认结论是否稳定。本工作将敏感性分析**前置并整合为分析框架的核心**，使其成为靶点和功能条目发现的有机组成部分：

- 参数空间中的高频靶点**就是**高置信度候选靶点
- LOO 扰动下保留率高的靶点**就是**不依赖单一成分的稳健靶点
- 参数空间中功能条目排名的稳定性**就是**通路层面结论的可信度指标

这种角色转换使敏感性分析从"事后检查"升级为"事前发现"，为网络药理学方法论提供了新的范式参考。

---

## 复现代码对网络药理学研究方法的贡献

### 1. 提供传统网络药理学的一键全流程分析脚本

本仓库将传统网络药理学分析全流程——化合物结构规范化 → PubChem/SwissTarget 靶点预测 → STRING PPI 网络构建 → Hub 基因筛选 → GO BP/KEGG 富集分析 → 可视化——封装为 `run_traditional_pipeline.py`。用户仅需准备化合物列表，即可一键完成从原始数据到可发表图表的全流程分析。

- 降低了网络药理学方法的使用门槛，尤其对非生物信息学背景的药学研究者友好
- 提供了标准化的分析流程，避免了不同研究者因参数选择差异导致的结论不可比问题
- 脚本可作为教学工具，帮助初学者理解网络药理学的完整分析逻辑链

### 2. 将网络药理学分析转化为参数鲁棒性和成分鲁棒性问题

本工作首次系统性地将**敏感性分析**方法引入网络药理学流程，将传统上被视为"技术细节"的阈值选择问题提升为方法学的核心议题。通过参数网格全因子设计和成分 LOO 设计，为网络药理学研究提供了**可量化的不确定性评估框架**。

- 使研究者能够回答"如果选择不同的阈值，结论会改变吗？"这一关键但常被忽略的问题
- 将网络药理学从"单一快照"升级为"参数空间中的稳健快照集合"
- 为审稿人和读者提供了评估研究结论可信度的量化依据

### 3. 提供成分扰动下的鲁棒性分析新思路

传统网络药理学将多成分体系的靶点预测结果视为"整体协同效应"，隐含假设所有成分均对最终结论有贡献，未考虑单个成分的差异化驱动作用。本工作提出的**成分 LOO 扰动分析**填补了这一方法论空白。

- 使研究者能够区分"体系的整体效应"和"个别成分的驱动效应"
- 识别体系中哪些化合物是"关键不可替代"的，哪些是"可被其他成分补偿"的
- 为多成分体系（如中药复方、植物提取物）的机制研究提供了新的分析维度
- 指导后续实验验证中的化合物优先级选择

### 4. 完整的可复现研究包

本仓库从方法论层面回应了计算生物学领域日益关注的"可复现性危机"，提供从原始数据到最终论文图表的全链条可复现代码：

| 可复现要素 | 具体实现 |
|------------|----------|
| **跨平台复现** | `reproduce.py`（Python）/ `reproduce.sh`（Bash）/ `reproduce.ps1`（PowerShell） |
| **断点续跑** | 2,016 个独立 JSON 检查点，支持中断后从断点继续 |
| **内置数据** | 仓库包含全部中间输入、缓存文件和检查点 |
| **分层复现** | 支持 `--core`（核心流程）和 `--full`（完整流程）两种复现深度 |
| **输出透明** | 详细的输出文件说明、图表目录和结果锁定文件 |

- 为网络药理学研究的透明性和可复现性设立了实践标准
- 使其他研究者能够精确复现本工作的每一个分析步骤和每一个数值
- 提供了一种"可复现论文"的参考模板

### 5. 模块化架构设计便于社区复用

代码采用分层模块化架构，各组件职责清晰、可独立复用：

```
connectors/   → API 封装层（SwissTarget, STRING, PubChem, PubMed）
enrichment/   → 功能富集层（GO/KEGG 富集与可视化）
visualization/ → 可视化层（论文级图表生成）
```

- 其他研究者可直接复用单个模块（如 SwissTarget 连接器、STRING PPI 构建器），无需理解整个框架
- 各模块的接口清晰、文档齐全，降低了二次开发的成本
- 模块可被替换或扩展（如将 STRING 替换为其他 PPI 数据库），框架保持灵活性

### 6. 方法通用性与领域迁移能力

虽然本工作以草莓叶抗氧化活性成分为应用案例，但整个双扰动框架是**完全通用**的——只需替换输入化合物列表和疾病靶点集，即可应用于：

- 中药复方的网络药理学机制解析
- 多酚类、黄酮类等天然产物混合物的靶点预测
- 任何多成分体系的靶点优先级排序和稳健性评估

### 7. 补充性贡献：实验验证指导、阴性结果价值与质量控制维度

除上述核心贡献外，本工作还提供了几项补充性方法论贡献：

- **实验验证指导**：成分影响度排序和靶点 Tier 分级可直接指导后续体外/体内实验的优先级——优先验证 Tier 1 靶点和关键驱动化合物，提高实验资源的利用效率。
- **阴性结果的价值化**：传统分析流程中"因阈值变化而消失的靶点"通常被静默丢弃。本框架使这些信息变得可见且可解释：化合物敏感型靶点提示研究者该靶点可能仅在特定成分存在时才具有生物学意义，而非"假阳性"。
- **质量控制维度**：本框架为网络药理学研究增加了"置信度"维度——研究者不仅可以报告"发现了哪些靶点"，更重要的是可以报告"每个靶点的发现有多可靠"（参数频率 + LOO 保留率）。
- **跨学科桥梁**：本工作将统计学中敏感性分析、不确定性量化的思想，以及计算机科学中鲁棒性测试的方法论，系统性地引入网络药理学领域，搭建了跨学科方法论的桥梁。

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
