# Stage 3 靶点分析记录

> **注意**: 本文档记录了 2026-04-03 的运行结果（109 个靶点）。后续修正多基因拆分后最新结果为 **116 个靶点**。详见 [`数据修正记录.md`](数据修正记录.md)。

## 日期: 2026-04-03

## 靶点统计结果

### 原始数据
- SwissTargetPrediction 5种成分靶点预测
- Probability > 0 过滤阈值（与原论文一致）

### 合并去重后靶点数量

| 来源 | 靶点数 |
|------|--------|
| **完整5成分** | **109** |
| 排除Catechin (4成分) | 89 |
| Catechin 新增贡献 | 20 |

### 各成分贡献明细

#### 完整5成分 (109个靶点)

| 成分 | 靶点数 |
|------|--------|
| Catechin (儿茶素) | 58 |
| Sucrose (蔗糖) | 23 |
| Scutellarin (野黄芩苷) | 12 |
| Procyanidin B1 (原花青素B1) | 9 |
| Trehalose (海藻糖) | 7 |

#### 排除Catechin后 (89个靶点)

| 成分 | 靶点数 |
|------|--------|
| Procyanidin B1 | 45 |
| Sucrose | 23 |
| Scutellarin | 14 |
| Trehalose | 7 |

### 重要发现

1. **儿茶素贡献20个独有靶点**
   - 排除儿茶素：89个靶点
   - 包含儿茶素：109个靶点
   - 新增：20个（主要为CA家族碳酸酐酶靶点）

2. **儿茶素概率覆盖效应**
   - 由于儿茶素预测概率高，很多靶点被其"抢走"
   - 排除儿茶素后，Procyanidin B1 从9个增至45个
   - 这是正常的去重逻辑（保留最高概率来源）

3. **与原论文对比**
   - 原论文（4成分，无儿茶素）：87个靶点
   - 本测试（4成分，无儿茶素）：89个靶点
   - 差异来源：数据库版本更新

### 去重逻辑说明

合并去重逻辑：
1. 同一靶点可能出现在多个成分中
2. 保留**概率最高**的成分作为来源
3. 过滤 **Probability > 0** 的靶点

这意味着：
- 如果靶点A在Catechin中概率=1.0，在Procyanidin B1中概率=0.5
- 最终记录：A, Catechin, 1.0
- 即使Procyanidin B1对该靶点也有贡献，也会被标记为Catechin来源

### 文件位置

- 原始CSV: `src/outputs/shared/swiss_targets/`
- 合并结果: `src/outputs/shared/merged_targets.csv`
- 疾病基因: `src/outputs/shared/disease_genes/merged_disease_genes.csv`

### 验证命令

```bash
python -m src.pipeline.stage3_merge_targets
```
