# Stage 4: 疾病基因数据库合并

## 概述

合并 GeneCards、OMIM、TTD 三个氧化应激相关基因数据库。

## 数据来源

### GeneCards
- 网址: https://www.genecards.org/
- 关键词: "oxidative stress"
- 搜索结果: 14,756 个基因
- 文件: `GeneCards-SearchResults.csv`

### OMIM
- 网址: https://www.omim.org/
- 关键词: "oxidative stress"
- 筛选: 只保留 `*` 开头的基因（代表基因，`#` 开头代表疾病表型）
- 文件: `OMIM-Entry-Search.csv`

### TTD (Therapeutic Target Database)
- 网址: http://db.idtbru.com/tdt/
- 关键词: "oxidative stress"
- 搜索结果: 2 条记录
- 文件: `TTD_Search_Results.xlsx`

## 数据处理

### GeneCards
- 直接提取 `Gene Symbol` 列

### OMIM
- Title 格式: `"描述文字; 基因符号"`
- 提取方法: 按分号分割，取后半部分
- 示例: `OXIDATIVE STRESS-RESPONSIVE PROTEIN, SERINE-RICH, 1; OSER1` → `OSER1`

### TTD
- Gene Names 格式: `"OXSR1,KIAA1101,OSR1"` (逗号分隔多个基因)
- 提取方法: 按逗号分割，逐个提取

## 合并结果

| 数据库 | 原始记录 | 去重后 |
|--------|----------|--------|
| GeneCards | 14,756 | 14,756 |
| OMIM | 171 | 9 (其余在GeneCards中已存在) |
| TTD | 3 | 0 (全部在GeneCards中已存在) |
| **总计** | 14,930 | **14,766** |

### OMIM 独有基因 (9个)

这些基因在 GeneCards 中不存在:
- ACTR5, C10ORF88, CDADC1, DJ1, ERO1L, LOXL1AS1, MTATP6, RECQL2, SMYD1

### TTD 分析
- 第一条记录: `OXSR1,KIAA1101,OSR1` → 3个基因，全部在GeneCards中
- 第二条记录: `N.A.` → 无效数据

## 去重逻辑

1. 合并三个数据库的基因
2. 按基因符号去重
3. 优先级: GeneCards > OMIM > TTD
4. 保留优先级最高的来源信息

## 输出文件

- `src/outputs/shared/disease_genes/merged_disease_genes.csv`

## 运行命令

```bash
python -m src.pipeline.stage4_merge_disease_genes
```

## 状态

✅ 已完成
