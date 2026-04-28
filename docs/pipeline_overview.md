# 网络药理学PPI筛选改进 - Pipeline 文档

## Pipeline 概述

```
Stage 1: PubChem获取SMILES
Stage 2: SwissTargetPrediction靶点预测
Stage 3: 靶点合并去重 + Probability>0过滤
Stage 4: 疾病基因数据库合并 (GeneCards + OMIM + TTD)
Stage 5: 靶点与疾病基因交集分析
Stage 6: PPI网络构建 (String DB)
Stage 7: 参数稳健性分析（Prob × STRING × TopK）
Stage 8: 富集分析与通路稳健性评估 (GO/KEGG)
```

---

## 文件结构

```
code/
├── docs/                           # 文档
│   ├── pipeline_overview.md        # 本文档
│   ├── stage1_pubchem.md           # Stage 1 记录
│   ├── stage2_swiss_prediction.md   # Stage 2 记录
│   ├── stage3_merge_targets.md      # Stage 3 记录
│   ├── stage4_merge_disease.md      # Stage 4 记录
│   ├── stage5_intersection.md       # Stage 5 记录
│   └── stage6_ppi_network.md       # Stage 6 记录
│
├── src/
│   ├── connectors/                  # 数据库连接器
│   │   ├── pubchem.py               # PubChem API
│   │   ├── swiss_target.py          # SwissTarget Playwright 自动化
│   │   ├── string_connector.py      # STRING API ✅
│   │   └── kegg_antioxidant.py      # KEGG 抗氧化通路筛选
│   │
│   ├── pipeline/                    # Pipeline脚本
│   │   ├── stage1_fetch.py          # Stage 1: PubChem 查询
│   │   ├── stage2_predict.py        # Stage 2: Swiss 靶点预测
│   │   ├── stage3_merge_targets.py
│   │   ├── stage4_merge_disease_genes.py
│   │   ├── stage5_intersection_analysis.py
│   │   ├── stage6_ppi_network.py    # ✅ 已实现
│   │   └── enrichment/              # Stage 8: 富集分析
│   │       └── enrichment_analysis.py
│   │
│   ├── visualization/               # 可视化模块
│   │   ├── ppi_network.py
│   │   ├── compound_target_pathway_network.py
│   │   ├── venn_diagram.py
│   │   └── enrichment_plotting.py
│   │
│   ├── workers/                     # 并发 Worker
│   │   ├── pubchem_worker.py
│   │   └── swiss_worker.py
│   │
│   └── outputs/                     # 输出结果
│       ├── shared/                  # 共享中间结果
│       │   ├── swiss_targets/       # 各成分靶点原始数据
│       │   │   ├── Catechin_targets.csv
│       │   │   ├── Procyanidin_B1_targets.csv
│       │   │   ├── Scutellarin_targets.csv
│       │   │   ├── Sucrose_targets.csv
│       │   │   └── Trehalose_targets.csv
│       │   ├── disease_genes/       # 疾病基因原始数据
│       │   │   ├── GeneCards-SearchResults.csv
│       │   │   ├── OMIM-Entry-Search.csv
│       │   │   └── TTD_Search_Results.xlsx
│       │   ├── merged_targets.csv           # Stage 3 输出
│       │   ├── merged_disease_genes.csv     # Stage 4 输出
│       │   ├── intersected_targets.csv      # Stage 5 输出
│       │   └── ppi_network/                 # Stage 6 输出 ✅
│       │       ├── string_interactions_interactions.csv
│       │       ├── string_interactions_nodes.csv
│       │       ├── string_interactions_degrees.csv
│       │       ├── network_nodes.csv
│       │       ├── node_degrees.csv
│       │       └── hub_nodes.csv
│       ├── traditional/             # 传统方法输出
│       │   ├── enrichment/
│       │   ├── figures/
│       │   └── evaluation/
│       └── sensitivity/             # 参数稳健性分析输出
│           ├── param_grid_results.csv
│           ├── target_robustness.csv
│           ├── pathway_rank_stability.csv
│           └── figures/
```

---

## 运行命令

```bash
# 运行单个Stage
python -m src.pipeline.stage1_fetch
python -m src.pipeline.stage2_predict
python -m src.pipeline.stage3_merge_targets
python -m src.pipeline.stage4_merge_disease_genes
python -m src.pipeline.stage5_intersection_analysis
python -m src.pipeline.stage6_ppi_network

# 运行完整流水线 (Stage 1-2)
python scripts/run_pipeline.py

# 运行传统方法 (Stage 3-6 + 分支A)
python scripts/run_traditional_pipeline.py

# 运行完整流水线 (Stage 1-6 + 分支A)
python scripts/run_full_pipeline.py

# 运行参数稳健性分析
python scripts/analysis/run_sensitivity_analysis.py
```

---

## 各Stage状态

| Stage | 状态 | 文档 |
|-------|------|------|
| Stage 1 | ✅ 完成 | `docs/stage1_pubchem.md` |
| Stage 2 | ✅ 完成 | `docs/stage2_swiss_prediction.md` |
| Stage 3 | ✅ 完成 | `docs/stage3_merge_targets.md` |
| Stage 4 | ✅ 完成 | `docs/stage4_merge_disease.md` |
| Stage 5 | ✅ 完成 | `docs/stage5_intersection.md` |
| Stage 6 | ✅ 完成 | `docs/stage6_ppi_network.md` |
| Stage 7 | ✅ 当前主线 (参数稳健性分析) | `scripts/analysis/run_sensitivity_analysis.py` |
| Stage 8 | ✅ 当前主线 (富集与通路稳健性) | `src/enrichment/enrichment_analysis.py` / `scripts/analysis/augment_pathway_robustness.py` |

> 注：早期 GNN/Hybrid 方案已转为探索记录，不作为当前投稿主工作流。
