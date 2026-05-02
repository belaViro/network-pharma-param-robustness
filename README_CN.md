# 网络药理学参数敏感性分析

[English](./README.md) | **中文**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> **网络药理学鲁棒靶点识别的双扰动框架：** 参数网格敏感性分析与逐化合物留一法验证，应用于草莓叶活性成分抗氧化的靶点识别。

## 概述

本仓库包含一项系统性网络药理学管道鲁棒性分析的完整可复现包。我们关注一个关键但常被忽视的问题：**参数敏感性**——同一组学数据在不同阈值设定下可能产生截然不同的靶点和通路结果。

### 核心贡献

- **双扰动框架（Dual-perturbation）**：同时扰动 (1) 分析参数 与 (2) 输入化合物，识别对两类变异均稳健的靶点。
- **112 组参数组合**：穷举概率阈值（7）× STRING 置信度（4）× Hub 基因截断（4）的全因子网格。
- **2,016 组化合物敏感性运行**：18 种草莓叶活性成分的逐化合物留一法（Leave-One-Out, LOO）分析。
- **三层靶点分类体系**：核心双稳健靶点（29 个）vs. 化合物特异性靶点 vs. 参数条件性靶点。

### 关键发现

| 指标 | 数值 |
|------|------|
| 参数组合总数 | **112**（7 × 4 × 4） |
| 化合物 LOO 组合总数 | **2,016**（18 × 112） |
| SwissTarget 并集靶点 | 684 |
| 疾病交集靶点 | 639 |
| Q1 高置信候选靶点 | 50 |
| **双稳健靶点** | **29**（Q1 保留率 58%） |
| **Tier 1 靶点** | **STAT3**（参数频率 100%，LOO 频率 100%） |
| 双稳健通路 | 58 |
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
```

需要 Python ≥ 3.10。

### 3. 运行主分析

```bash
# 参数网格敏感性分析（112 组组合）
python scripts/analysis/run_sensitivity_analysis.py

# 化合物留一法敏感性分析（2,016 组组合）
python scripts/analysis/run_compound_sensitivity_analysis.py
```

### 4. 一键复现

```bash
bash reproduce.sh
```

该脚本依次执行：参数网格分析 → 稳健 vs. 参数化对比 → 通路象限增强 → 论文图表重生成。

---

## 仓库结构

```
.
├── scripts/
│   ├── analysis/
│   │   ├── run_sensitivity_analysis.py           # 主程序：112 组参数网格
│   │   ├── run_compound_sensitivity_analysis.py  # 主程序：18 化合物 LOO
│   │   ├── compare_robust_vs_parametric.py       # 稳健 vs. 参数化对比
│   │   ├── augment_pathway_robustness.py         # 通路稳健性象限增强
│   │   ├── plot_sensitivity_figures.py           # 论文图表生成
│   │   ├── validate_with_pubmed.py               # 外部验证（PubMed）
│   │   └── ...                                   # 其他诊断脚本
│   ├── run_full_pipeline.py                      # 端到端流水线
│   └── run_traditional_pipeline.py               # 基线方法对比
│
├── src/
│   ├── config.py                                 # 全局配置中心
│   ├── connectors/                               # API 封装
│   │   ├── swiss_target.py                       # SwissTargetPrediction
│   │   ├── string_connector.py                   # STRING PPI
│   │   ├── pubchem.py                            # PubChem
│   │   └── pubmed_connector.py                   # PubMed 验证
│   ├── enrichment/
│   │   ├── enrichment_analysis.py                # GO/KEGG 富集（gseapy）
│   │   └── plotting.py                           # 富集可视化
│   ├── visualization/                            # 绘图工具
│   ├── data/                                     # 参考数据
│   │   ├── go-basic.obo                          # Gene Ontology 本体
│   │   └── kegg_pathway_classes.json             # KEGG 通路分类
│   └── outputs/                                  # 全部分析输出
│       ├── sensitivity/                          # 参数敏感性结果
│       │   ├── target_robustness.csv
│       │   ├── pathway_robustness.csv
│       │   ├── param_grid_results.csv
│       │   └── compound_sensitivity/             # LOO 化合物结果
│       ├── shared/                               # 跨阶段共享数据
│       └── traditional/                          # 基线方法输出
│
├── tests/                                        # 单元与集成测试
├── docs/                                         # 项目文档
├── requirements.txt
├── pyproject.toml
└── reproduce.sh                                  # 一键复现脚本
```

---

## 输出文件

### 核心结果（参数敏感性）

| 文件 | 说明 |
|------|------|
| `src/outputs/sensitivity/target_robustness.csv` | 112 组组合下各靶点的频率与度统计 |
| `src/outputs/sensitivity/pathway_robustness.csv` | 通路稳健性指标（信号主干、直接氧化、交叉调控） |
| `src/outputs/sensitivity/param_grid_results.csv` | 每组组合的富集汇总 |
| `src/outputs/sensitivity/pathway_rank_stability.csv` | 通路排名在参数空间中的一致性 |

### 核心结果（化合物敏感性）

| 文件 | 说明 |
|------|------|
| `src/outputs/sensitivity/compound_sensitivity/compound_influence_scores.csv` | 各化合物的影响度评分（0–1） |
| `src/outputs/sensitivity/compound_sensitivity/dual_robust_targets.csv` | 29 个双稳健靶点及其 LOO 保留率 |
| `src/outputs/sensitivity/compound_sensitivity/loo_pathway_robustness.csv` | 化合物剔除后的通路保留情况 |

### 检查点

- `src/outputs/sensitivity/compound_sensitivity/checkpoints/combo_results/` — 2,016 个独立 JSON 文件（每个化合物 × 组合一份），支持全量审计与断点续跑。

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

LC-MS/MS 共鉴定 19 条原始记录；去重后（Trehalose ≡ alpha,alpha-Trehalose）得到 18 个 canonical 化合物：

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
15. alpha-Linolenic acid（α-亚麻酸）
16. Palmitic acid（棕榈酸）
17. Linoleic acid（亚油酸）
18. Stearic acid（硬脂酸）

---

## 靶点分类体系

| 层级 | 判定标准 | 数量 | 代表靶点 |
|------|----------|------|----------|
| **Tier 1** | 参数频率 = 100% 且 LOO 频率 = 100% | 1 | **STAT3** |
| **Tier 2** | 双稳健（参数稳健 & 化合物稳健） | 28 | AKT1, MAPK3, EGFR, ... |
| **化合物特异性** | 剔除强心苷后 LOO 频率降至 0 | — | CASP3, MAPK1, PIK3CA |
| **参数条件性** | 低频率双稳健 | — | MET, IL2, KIT |

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
- **完整可复现包**（约 33 MB）：可通过 GitHub Release [v1.0.0](https://github.com/belaViro/network-pharma-param-robustness/releases) 下载

使用的外部数据库：
- SwissTargetPrediction（http://www.swisstargetprediction.ch/）
- STRING v12（https://string-db.org/）
- Enrichr（https://maayanlab.cloud/Enrichr/）
- Gene Ontology（http://purl.obolibrary.org/obo/go-basic.obo）

---

## 引用

若使用本代码或框架，请引用：

> *Parameter robustness analysis for network pharmacology target identification: a dual-perturbation framework applied to strawberry leaf compounds against oxidative stress.*

---

## 许可证

[MIT License](./LICENSE)
