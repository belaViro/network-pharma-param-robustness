# 参数网格设计文档

## 当前主分析参数网格

本文件记录当前锁定的参数敏感性分析网格。所有论文、README、补充材料和结果解释均应以本页数字为准。

| 参数 | 含义 | 当前取值 | 组合数 |
|:---|:---|:---|---:|
| `Prob` | SwissTargetPrediction Probability 阈值 | `[0, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]` | 7 |
| `STRING` | STRING combined_score 阈值 | `[0.40, 0.50, 0.70, 0.90]` | 4 |
| `TopK` | PPI Degree 排名前 K 个 hub genes | `[10, 20, 30, 50]` | 4 |
| `Pvalue` | 富集分析显著性阈值 | `0.01` 固定 | 1 |
| **总计** |  |  | **112** |

## 当前锁定口径

- 参数组合总数：112。
- 成分扰动分析也使用同一 112 组合网格。
- 18 个 canonical 化合物 leave-one-compound-out 后总运行数：18 × 112 = 2016。
- Q1 candidate robust targets：50。
- dual-robust targets：29。
- Tier 1 dual-robust target：1，STAT3。
- Q1 candidate robust functional terms：70。
- dual-robust functional terms：58。

## 设计说明

1. `Prob=0` 使用 `Probability > threshold` 逻辑，因此表示保留所有正概率预测靶点。
2. `Prob=0.15/0.20` 保留在主网格中作为压力测试点，用于观察高阈值导致的靶点池收缩和富集边界。
3. `STRING` 使用 0.40、0.50、0.70、0.90，覆盖 STRING 中等到高置信度区间。
4. `TopK` 使用 10、20、30、50，覆盖紧缩核心、常用 Top30 和较宽松 hub 集合。
5. `Pvalue` 不参与多取值扫描，因为它影响富集显著性判断，不影响上游靶点筛选；如需审稿回应，可单独做 `0.05` 敏感性说明。

## 历史口径处理

早期 36 组合、72 组合、33 个 Q1 靶点、9/13 个 Tier 1 等均为探索阶段结果，不进入当前主文、补充表或结果锁定口径。若旧报告中出现这些数字，应标记为历史探索，不得与 112 组合主结果混用。
