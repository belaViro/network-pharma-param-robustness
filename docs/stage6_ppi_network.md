# Stage 6: PPI网络构建 (STRING Database)

## 概述

使用STRING Database API构建蛋白质相互作用(PPI)网络，获取112个交集靶点的蛋白互作关系。

## STRING API配置

- **API端点**: `https://string-db.org/api/json/network`
- **物种**: Homo sapiens (NCBI Taxonomy ID: 9606)
- **最低置信度**: 0.4 (默认阈值)

## API调用参数

```python
params = {
    "species": 9606,
    "identifiers": gene_names_separated_by_newline,
    "format": "json",
}
```

## 网络构建结果

| 指标 | 数值 |
|------|------|
| 输入靶点数 | 112 |
| 网络节点数 | 111 |
| 相互作用数 | 830 |
| 平均度数 | 14.95 |

**说明**: 1个靶点基因在STRING数据库中无相互作用记录

## Top 20 Hub基因 (按度数排序)

| 排名 | 基因 | 度数 | 功能相关 |
|------|------|------|----------|
| 1 | TNF | 60 | 肿瘤坏死因子，炎症反应核心 |
| 2 | SRC | 56 | 酪氨酸激酶，信号转导 |
| 3 | EGFR | 50 | 表皮生长因子受体 |
| 4 | BCL2 | 47 | 抗凋亡蛋白 |
| 5 | IL6 | 42 | 白细胞介素6，炎症因子 |
| 6 | MMP9 | 41 | 基质金属蛋白酶9 |
| 7 | STAT3 | 38 | 转录因子，氧化应激响应 |
| 8 | HIF1A | 37 | 低氧诱导因子 |
| 9 | MMP2 | 36 | 基质金属蛋白酶2 |
| 10 | VEGFA | 35 | 血管内皮生长因子A |
| 11 | PTGS2 | 34 | 环氧化酶2 (COX-2) |
| 12 | ESR1 | 33 | 雌激素受体 |
| 13 | KDR | 32 | 血管内皮生长因子受体2 |
| 14 | ERBB2 | 31 | HER2受体 |
| 15 | HSP90AA1 | 30 | 热休克蛋白 |
| 16 | FGF2 | 28 | 成纤维细胞生长因子 |
| 17 | APP | 27 | 淀粉样前体蛋白 |
| 18 | IL2 | 25 | 白细胞介素2 |
| 19 | KIT | 24 | 干细胞因子受体 |
| 20 | MAPK14 | 23 | MAP激酶 |

## 输出文件

### 核心文件 (位于 `src/outputs/shared/ppi_network/`)
- `string_interactions_interactions.csv` - 相互作用列表
  - 列: protein1, protein2, combined_score

- `string_interactions_nodes.csv` - 网络节点列表
  - 列: Gene

- `string_interactions_degrees.csv` - 节点度数统计
  - 列: Gene, Degree

- `network_nodes.csv` - 网络节点列表（兼容副本）
- `node_degrees.csv` - 节点度数统计（兼容副本）
- `hub_nodes.csv` - Top Hub 基因列表

## 网络特征分析

1. **高度数节点集中于**:
   - 炎症相关: TNF, PTGS2 (COX-2), MMP9, MMP2
   - 信号转导: SRC, EGFR, STAT3, MAPK14
   - 血管生成: HIF1A, KDR, VEGFA, FGF2
   - 细胞凋亡: BCL2

2. **氧化应激关联**: 多个hub基因与氧化应激直接相关(TNF, HIF1A, NRF2等)

## 运行命令

```bash
# 运行Stage 6
python -m src.pipeline.stage6_ppi_network

# 指定最低置信度
python -m src.pipeline.stage6_ppi_network --min-score 0.7

# 指定靶点文件
python -m src.pipeline.stage6_ppi_network --targets path/to/targets.csv
```

## 后续处理

Stage 7将使用以下方法进行节点筛选对比:
- **GNN图神经网络筛选**: 基于图结构学习节点重要性
- **Degree Centrality基准**: 基于度数中心性筛选

消融实验设计:
- 实验组: GNN筛选出的核心节点
- 对照组: Degree Centrality筛选的节点
- 评估: GO/KEGG富集分析比较

## 状态

✅ 已完成
