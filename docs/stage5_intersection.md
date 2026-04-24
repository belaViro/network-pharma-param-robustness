# Stage 5: 靶点与疾病基因交集分析

## 概述

将 SwissTargetPrediction 预测的靶点与氧化应激疾病基因取交集，筛选出与氧化应激相关的靶点。

## 输入

1. **靶点数据**: `src/outputs/shared/merged_targets.csv`
   - 来源: SwissTargetPrediction
   - 数量: 116 个靶点 (Probability > 0，已修正多基因拆分)

2. **疾病基因数据**: `src/outputs/shared/disease_genes/merged_disease_genes.csv`
   - 来源: GeneCards + OMIM + TTD
   - 数量: 14,766 个基因

## 交集结果

| 项目 | 数量 |
|------|------|
| SwissTargetPrediction 靶点 | 116 |
| 疾病基因 | 14,766 |
| **交集靶点** | **112** |

**覆盖率**: 112/116 = 96.6% 的预测靶点与氧化应激相关

## 成分来源分布

| 成分 | 交集靶点数 |
|------|-----------|
| Catechin (儿茶素) | 56 |
| Sucrose (蔗糖) | 23 |
| Scutellarin (野黄芩苷) | 12 |
| Procyanidin B1 (原花青素B1) | 9 |
| Trehalose (海藻糖) | 4 |

## 高概率交集靶点 (Top 20)

| 靶点 | 来源成分 | 概率 | 疾病来源 |
|------|----------|------|----------|
| CA7 | Catechin | 1.0 | GeneCards |
| CA6 | Catechin | 1.0 | GeneCards |
| CA2 | Catechin | 1.0 | GeneCards |
| CA12 | Catechin | 1.0 | GeneCards |
| CA1 | Catechin | 1.0 | GeneCards |
| CA5A | Catechin | 1.0 | GeneCards |
| CA4 | Catechin | 1.0 | GeneCards |
| CA5B | Catechin | 1.0 | GeneCards |
| CA9 | Catechin | 1.0 | GeneCards |
| CA3 | Catechin | 1.0 | GeneCards |
| AKR1B1 | Scutellarin | 0.712 | GeneCards |
| MMP2 | Procyanidin B1 | 0.577 | GeneCards |
| BACE1 | Catechin | 0.567 | GeneCards |
| PTGS1 | Catechin | 0.518 | GeneCards |
| PGF | Procyanidin B1 | 0.352 | GeneCards |
| VEGFA | Procyanidin B1 | 0.352 | GeneCards |
| CDK1 | Trehalose | 0.243 | GeneCards |
| MMP9 | Procyanidin B1 | 0.208 | GeneCards |
| ABCB1 | Catechin | 0.156 | GeneCards |
| BCL2 | Catechin | 0.156 | GeneCards |

## 分析

1. **儿茶素主导**: CA家族(碳酸酐酶)全部由儿茶素贡献，概率均为1.0
2. **疾病相关性高**: 96.6%的预测靶点与氧化应激相关
3. **血管生成相关**: VEGFA, PGF, MMP2, MMP9 等参与血管生成

## 输出文件

- `src/outputs/shared/intersected_targets.csv`

## 运行命令

```bash
python -m src.pipeline.stage5_intersection_analysis
```

## 状态

✅ 已完成
