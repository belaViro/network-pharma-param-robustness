# 参数网格敏感性分析报告

**生成时间**: 2026-05-07 14:56:35

## 1. 参数网格设置

| 参数 | 取值 | 组合数 |
|------|------|--------|
| Prob (SwissTarget) | [0.0, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2] | 7 |
| STRING (PPI) | [0.4, 0.5, 0.7, 0.9] | 4 |
| TopK (Hub基因) | [10, 20, 30, 50] | 4 |
| Pvalue 固定 | 0.01 | 1 |
| **总计** | | **112** |

## 2. 数据加载概况

- SwissTarget 并集: 684 靶点
- 疾病基因: 14765 个
- 疾病交集后: 639 靶点
- STRING PPI 边数: 3180
- 成功组合: 112 / 112

## 3. 各概率阈值靶点数

| Prob 阈值 | 疾病相关靶点数 |
|-----------|---------------|
| > 0.0 | 286 |
| > 0.03 | 272 |
| > 0.05 | 219 |
| > 0.08 | 194 |
| > 0.1 | 156 |
| > 0.15 | 42 |
| > 0.2 | 33 |

## 4. 每组参数富集分析详情

| 组合 | 靶点数 | PPI节点 | PPI边 | 显著通路 | 主干道 | 交叉调控 | 直接OX | 其他 | 平均EF | 最小P |
|------|--------|--------|-------|---------|-------|---------|--------|-------|--------|-------|
| P0.0_S0.4_K10 | 286 | 284 | 3157 | 40 | 11 | 14 | 0 | 15 | 100.7 | 6.71e-13 |
| P0.0_S0.4_K20 | 286 | 284 | 3157 | 40 | 13 | 15 | 0 | 12 | 47.7 | 1.17e-19 |
| P0.0_S0.4_K30 | 286 | 284 | 3157 | 40 | 15 | 13 | 0 | 12 | 43.5 | 1.58e-24 |
| P0.0_S0.4_K50 | 286 | 284 | 3157 | 40 | 15 | 13 | 0 | 12 | 34.2 | 3.90e-29 |
| P0.0_S0.5_K10 | 286 | 278 | 2095 | 40 | 12 | 13 | 0 | 15 | 90.2 | 6.57e-13 |
| P0.0_S0.5_K20 | 286 | 278 | 2095 | 40 | 12 | 14 | 0 | 14 | 52.9 | 5.72e-22 |
| P0.0_S0.5_K30 | 286 | 278 | 2095 | 40 | 14 | 14 | 0 | 12 | 45.4 | 1.57e-24 |
| P0.0_S0.5_K50 | 286 | 278 | 2095 | 40 | 15 | 14 | 0 | 11 | 34.5 | 2.27e-34 |
| P0.0_S0.7_K10 | 286 | 260 | 919 | 40 | 8 | 10 | 0 | 22 | 96.0 | 6.43e-13 |
| P0.0_S0.7_K20 | 286 | 260 | 919 | 40 | 10 | 16 | 1 | 13 | 61.5 | 5.72e-22 |
| P0.0_S0.7_K30 | 286 | 260 | 919 | 40 | 14 | 15 | 1 | 10 | 39.3 | 8.93e-23 |
| P0.0_S0.7_K50 | 286 | 260 | 919 | 40 | 14 | 15 | 0 | 11 | 34.0 | 4.67e-38 |
| P0.0_S0.9_K10 | 286 | 205 | 370 | 40 | 9 | 12 | 0 | 19 | 102.3 | 1.08e-12 |
| P0.0_S0.9_K20 | 286 | 205 | 370 | 40 | 13 | 11 | 1 | 15 | 55.2 | 4.40e-18 |
| P0.0_S0.9_K30 | 286 | 205 | 370 | 40 | 13 | 12 | 3 | 12 | 41.1 | 1.14e-24 |
| P0.0_S0.9_K50 | 286 | 205 | 370 | 40 | 15 | 12 | 1 | 12 | 27.7 | 1.30e-32 |
| P0.03_S0.4_K10 | 272 | 270 | 2884 | 40 | 11 | 14 | 0 | 15 | 100.7 | 6.71e-13 |
| P0.03_S0.4_K20 | 272 | 270 | 2884 | 40 | 11 | 17 | 0 | 12 | 45.2 | 5.69e-22 |
| P0.03_S0.4_K30 | 272 | 270 | 2884 | 40 | 15 | 12 | 0 | 13 | 44.3 | 7.88e-25 |
| P0.03_S0.4_K50 | 272 | 270 | 2884 | 40 | 15 | 13 | 0 | 12 | 32.8 | 1.36e-32 |
| P0.03_S0.5_K10 | 272 | 264 | 1911 | 40 | 12 | 13 | 0 | 15 | 90.2 | 6.57e-13 |
| P0.03_S0.5_K20 | 272 | 264 | 1911 | 40 | 12 | 14 | 0 | 14 | 52.9 | 5.72e-22 |
| P0.03_S0.5_K30 | 272 | 264 | 1911 | 40 | 14 | 13 | 0 | 13 | 43.5 | 7.79e-25 |
| P0.03_S0.5_K50 | 272 | 264 | 1911 | 40 | 14 | 14 | 0 | 12 | 33.0 | 2.26e-34 |
| P0.03_S0.7_K10 | 272 | 249 | 858 | 40 | 8 | 10 | 0 | 22 | 96.0 | 6.43e-13 |
| P0.03_S0.7_K20 | 272 | 249 | 858 | 40 | 10 | 16 | 1 | 13 | 61.5 | 5.72e-22 |
| P0.03_S0.7_K30 | 272 | 249 | 858 | 40 | 15 | 13 | 2 | 10 | 38.0 | 1.11e-24 |
| P0.03_S0.7_K50 | 272 | 249 | 858 | 40 | 14 | 14 | 0 | 12 | 37.0 | 5.97e-40 |
| P0.03_S0.9_K10 | 272 | 197 | 355 | 40 | 9 | 12 | 0 | 19 | 102.3 | 1.08e-12 |
| P0.03_S0.9_K20 | 272 | 197 | 355 | 40 | 13 | 11 | 1 | 15 | 55.2 | 4.40e-18 |
| P0.03_S0.9_K30 | 272 | 197 | 355 | 40 | 13 | 12 | 3 | 12 | 41.1 | 1.14e-24 |
| P0.03_S0.9_K50 | 272 | 197 | 355 | 40 | 15 | 12 | 1 | 12 | 44.1 | 1.30e-32 |
| P0.05_S0.4_K10 | 219 | 218 | 2159 | 40 | 11 | 14 | 0 | 15 | 100.7 | 6.71e-13 |
| P0.05_S0.4_K20 | 219 | 218 | 2159 | 40 | 17 | 13 | 0 | 10 | 43.8 | 4.46e-18 |
| P0.05_S0.4_K30 | 219 | 218 | 2159 | 40 | 13 | 13 | 0 | 14 | 44.7 | 1.25e-26 |
| P0.05_S0.4_K50 | 219 | 218 | 2159 | 40 | 17 | 13 | 0 | 10 | 31.9 | 2.24e-34 |
| P0.05_S0.5_K10 | 219 | 213 | 1442 | 40 | 12 | 13 | 0 | 15 | 90.2 | 6.57e-13 |
| P0.05_S0.5_K20 | 219 | 213 | 1442 | 40 | 10 | 14 | 0 | 16 | 55.4 | 2.21e-24 |
| P0.05_S0.5_K30 | 219 | 213 | 1442 | 40 | 14 | 15 | 0 | 11 | 39.8 | 1.19e-26 |
| P0.05_S0.5_K50 | 219 | 213 | 1442 | 40 | 14 | 15 | 0 | 11 | 29.6 | 3.43e-36 |
| P0.05_S0.7_K10 | 219 | 203 | 665 | 40 | 12 | 11 | 0 | 17 | 98.6 | 7.42e-13 |
| P0.05_S0.7_K20 | 219 | 203 | 665 | 40 | 10 | 16 | 1 | 13 | 61.5 | 5.72e-22 |
| P0.05_S0.7_K30 | 219 | 203 | 665 | 40 | 14 | 13 | 1 | 12 | 38.7 | 7.38e-25 |
| P0.05_S0.7_K50 | 219 | 203 | 665 | 40 | 13 | 14 | 0 | 13 | 30.1 | 5.91e-40 |
| P0.05_S0.9_K10 | 219 | 160 | 279 | 40 | 9 | 12 | 0 | 19 | 102.3 | 1.08e-12 |
| P0.05_S0.9_K20 | 219 | 160 | 279 | 40 | 10 | 13 | 1 | 16 | 57.2 | 4.40e-18 |
| P0.05_S0.9_K30 | 219 | 160 | 279 | 40 | 10 | 11 | 1 | 18 | 41.0 | 9.08e-23 |
| P0.05_S0.9_K50 | 219 | 160 | 279 | 40 | 13 | 12 | 0 | 15 | 33.1 | 3.18e-36 |
| P0.08_S0.4_K10 | 194 | 193 | 1845 | 40 | 11 | 14 | 0 | 15 | 100.7 | 6.71e-13 |
| P0.08_S0.4_K20 | 194 | 193 | 1845 | 40 | 14 | 15 | 0 | 11 | 40.7 | 1.07e-19 |
| P0.08_S0.4_K30 | 194 | 193 | 1845 | 40 | 13 | 14 | 0 | 13 | 45.1 | 1.14e-28 |
| P0.08_S0.4_K50 | 194 | 193 | 1845 | 40 | 17 | 12 | 0 | 11 | 31.8 | 2.25e-34 |
| P0.08_S0.5_K10 | 194 | 189 | 1242 | 40 | 12 | 13 | 0 | 15 | 90.2 | 6.57e-13 |
| P0.08_S0.5_K20 | 194 | 189 | 1242 | 40 | 10 | 14 | 0 | 16 | 56.4 | 2.22e-24 |
| P0.08_S0.5_K30 | 194 | 189 | 1242 | 40 | 14 | 15 | 0 | 11 | 39.8 | 1.19e-26 |
| P0.08_S0.5_K50 | 194 | 189 | 1242 | 40 | 14 | 13 | 0 | 13 | 32.3 | 6.10e-40 |
| P0.08_S0.7_K10 | 194 | 180 | 591 | 40 | 8 | 10 | 0 | 22 | 96.0 | 6.43e-13 |
| P0.08_S0.7_K20 | 194 | 180 | 591 | 40 | 10 | 16 | 1 | 13 | 61.5 | 5.72e-22 |
| P0.08_S0.7_K30 | 194 | 180 | 591 | 40 | 13 | 15 | 2 | 10 | 39.6 | 7.38e-25 |
| P0.08_S0.7_K50 | 194 | 180 | 591 | 40 | 14 | 13 | 0 | 13 | 31.8 | 5.85e-40 |
| P0.08_S0.9_K10 | 194 | 142 | 251 | 40 | 9 | 12 | 0 | 19 | 102.3 | 1.08e-12 |
| P0.08_S0.9_K20 | 194 | 142 | 251 | 40 | 9 | 12 | 1 | 18 | 55.6 | 1.23e-19 |
| P0.08_S0.9_K30 | 194 | 142 | 251 | 40 | 10 | 12 | 1 | 17 | 38.3 | 9.08e-23 |
| P0.08_S0.9_K50 | 194 | 142 | 251 | 40 | 15 | 13 | 0 | 12 | 45.9 | 1.28e-32 |
| P0.1_S0.4_K10 | 156 | 155 | 1191 | 40 | 11 | 15 | 0 | 14 | 87.7 | 1.25e-10 |
| P0.1_S0.4_K20 | 156 | 155 | 1191 | 40 | 10 | 14 | 0 | 16 | 52.0 | 1.20e-19 |
| P0.1_S0.4_K30 | 156 | 155 | 1191 | 40 | 15 | 13 | 0 | 12 | 41.6 | 1.14e-28 |
| P0.1_S0.4_K50 | 156 | 155 | 1191 | 40 | 15 | 13 | 0 | 12 | 27.8 | 2.24e-34 |
| P0.1_S0.5_K10 | 156 | 151 | 818 | 40 | 11 | 15 | 0 | 14 | 87.7 | 1.25e-10 |
| P0.1_S0.5_K20 | 156 | 151 | 818 | 40 | 13 | 13 | 0 | 14 | 49.9 | 5.75e-22 |
| P0.1_S0.5_K30 | 156 | 151 | 818 | 40 | 15 | 12 | 0 | 13 | 37.6 | 1.26e-26 |
| P0.1_S0.5_K50 | 156 | 151 | 818 | 40 | 14 | 12 | 0 | 14 | 28.8 | 1.33e-32 |
| P0.1_S0.7_K10 | 156 | 143 | 395 | 40 | 10 | 14 | 0 | 16 | 77.0 | 1.25e-10 |
| P0.1_S0.7_K20 | 156 | 143 | 395 | 40 | 12 | 15 | 0 | 13 | 52.4 | 1.13e-19 |
| P0.1_S0.7_K30 | 156 | 143 | 395 | 40 | 14 | 14 | 0 | 12 | 37.2 | 1.15e-28 |
| P0.1_S0.7_K50 | 156 | 143 | 395 | 40 | 14 | 12 | 0 | 14 | 34.5 | 1.30e-32 |
| P0.1_S0.9_K10 | 156 | 110 | 170 | 40 | 11 | 11 | 0 | 18 | 99.3 | 1.08e-12 |
| P0.1_S0.9_K20 | 156 | 110 | 170 | 40 | 14 | 9 | 0 | 17 | 55.4 | 1.13e-18 |
| P0.1_S0.9_K30 | 156 | 110 | 170 | 40 | 16 | 10 | 0 | 14 | 42.8 | 2.35e-26 |
| P0.1_S0.9_K50 | 156 | 110 | 170 | 40 | 15 | 12 | 0 | 13 | 45.8 | 3.63e-29 |
| P0.15_S0.4_K10 | 42 | 41 | 98 | 40 | 6 | 15 | 0 | 19 | 89.4 | 5.30e-06 |
| P0.15_S0.4_K20 | 42 | 41 | 98 | 25 | 2 | 10 | 0 | 13 | 91.4 | 5.61e-16 |
| P0.15_S0.4_K30 | 42 | 41 | 98 | 30 | 3 | 5 | 0 | 22 | 96.1 | 2.64e-20 |
| P0.15_S0.4_K50 | 42 | 41 | 98 | 31 | 3 | 6 | 0 | 22 | 75.4 | 9.80e-22 |
| P0.15_S0.5_K10 | 42 | 40 | 65 | 40 | 6 | 15 | 0 | 19 | 98.0 | 6.43e-09 |
| P0.15_S0.5_K20 | 42 | 40 | 65 | 32 | 3 | 10 | 0 | 19 | 102.4 | 4.99e-13 |
| P0.15_S0.5_K30 | 42 | 40 | 65 | 33 | 4 | 7 | 0 | 22 | 85.9 | 2.62e-20 |
| P0.15_S0.5_K50 | 42 | 40 | 65 | 31 | 2 | 6 | 0 | 23 | 78.7 | 7.36e-22 |
| P0.15_S0.7_K10 | 42 | 37 | 29 | 40 | 6 | 16 | 0 | 18 | 132.9 | 3.75e-06 |
| P0.15_S0.7_K20 | 42 | 37 | 29 | 39 | 4 | 12 | 0 | 23 | 105.7 | 5.61e-16 |
| P0.15_S0.7_K30 | 42 | 37 | 29 | 32 | 3 | 7 | 0 | 22 | 94.9 | 2.50e-20 |
| P0.15_S0.7_K50 | 42 | 37 | 29 | 31 | 3 | 5 | 0 | 23 | 84.8 | 2.65e-19 |
| P0.15_S0.9_K10 | 42 | 25 | 12 | 40 | 7 | 15 | 0 | 18 | 143.7 | 2.36e-09 |
| P0.15_S0.9_K20 | 42 | 25 | 12 | 40 | 4 | 12 | 0 | 24 | 103.8 | 8.09e-08 |
| P0.15_S0.9_K30 | 42 | 25 | 12 | 33 | 4 | 8 | 0 | 21 | 91.7 | 8.28e-09 |
| P0.15_S0.9_K50 | 42 | 25 | 12 | 33 | 4 | 8 | 0 | 21 | 91.7 | 8.28e-09 |
| P0.2_S0.4_K10 | 33 | 32 | 58 | 31 | 2 | 10 | 0 | 19 | 117.5 | 3.62e-12 |
| P0.2_S0.4_K20 | 33 | 32 | 58 | 20 | 3 | 2 | 0 | 15 | 116.1 | 2.31e-22 |
| P0.2_S0.4_K30 | 33 | 32 | 58 | 12 | 0 | 3 | 0 | 9 | 97.6 | 2.20e-23 |
| P0.2_S0.4_K50 | 33 | 32 | 58 | 7 | 0 | 2 | 0 | 5 | 112.1 | 4.82e-23 |
| P0.2_S0.5_K10 | 33 | 31 | 36 | 33 | 2 | 10 | 0 | 21 | 123.2 | 1.81e-15 |
| P0.2_S0.5_K20 | 33 | 31 | 36 | 13 | 1 | 2 | 0 | 10 | 129.5 | 2.39e-22 |
| P0.2_S0.5_K30 | 33 | 31 | 36 | 12 | 0 | 3 | 0 | 9 | 97.6 | 2.22e-23 |
| P0.2_S0.5_K50 | 33 | 31 | 36 | 11 | 0 | 2 | 0 | 9 | 101.8 | 3.28e-23 |
| P0.2_S0.7_K10 | 33 | 28 | 15 | 24 | 3 | 4 | 0 | 17 | 125.8 | 5.16e-09 |
| P0.2_S0.7_K20 | 33 | 28 | 15 | 26 | 3 | 3 | 0 | 20 | 100.3 | 2.08e-22 |
| P0.2_S0.7_K30 | 33 | 28 | 15 | 15 | 1 | 2 | 0 | 12 | 99.0 | 1.24e-20 |
| P0.2_S0.7_K50 | 33 | 28 | 15 | 15 | 1 | 2 | 0 | 12 | 99.0 | 1.24e-20 |
| P0.2_S0.9_K10 | 33 | 17 | 5 | 36 | 4 | 11 | 0 | 21 | 102.0 | 4.64e-06 |
| P0.2_S0.9_K20 | 33 | 17 | 5 | 17 | 2 | 1 | 0 | 14 | 103.6 | 1.44e-05 |
| P0.2_S0.9_K30 | 33 | 17 | 5 | 17 | 2 | 1 | 0 | 14 | 103.6 | 1.44e-05 |
| P0.2_S0.9_K50 | 33 | 17 | 5 | 17 | 2 | 1 | 0 | 14 | 103.6 | 1.44e-05 |

### 4.1 各组合富集通路明细

<details>
<summary><b>P0.03_S0.4_K10</b> — 显著40条, 主干11, 交叉14, 直接OX0, EF=100.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.71e-13 | 80.0 |
| Pathways in cancer | 交叉调控 | 7.20e-10 | 30.9 |
| Hepatitis B | 交叉调控 | 2.51e-09 | 75.9 |
| Lipid and atherosclerosis | 交叉调控 | 1.04e-08 | 57.2 |
| Human cytomegalovirus infection | 其他 | 1.10e-08 | 54.7 |
| Chemical carcinogenesis | 交叉调控 | 1.32e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.55e-08 | 26.4 |
| MicroRNAs in cancer | 交叉调控 | 5.38e-08 | 39.7 |
| Estrogen signaling pathway | 交叉调控 | 6.14e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.85e-07 | 29.8 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.06e-07 | 53.1 |
| IL-17 signaling pathway | 交叉调控 | 1.34e-06 | 87.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.57e-06 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.95e-06 | 75.2 |
| TNF signaling pathway | 交叉调控 | 1.95e-06 | 73.2 |
| Toxoplasmosis | 其他 | 1.95e-06 | 73.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.33e-06 | 59.0 |
| negative regulation of production of miRNAs involved in | 其他 | 4.77e-06 | 512.5 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 6.39e-06 | 25.4 |
| negative regulation of apoptotic signaling pathway (GO: | 信号主干通路 | 6.39e-06 | 105.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.4_K20</b> — 显著40条, 主干11, 交叉17, 直接OX0, EF=45.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.69e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 2.23e-18 | 29.0 |
| Chemical carcinogenesis | 交叉调控 | 4.69e-13 | 42.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.51e-12 | 21.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.91e-12 | 25.5 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.96e-12 | 47.8 |
| Lipid and atherosclerosis | 交叉调控 | 8.45e-12 | 42.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.02e-11 | 27.9 |
| Estrogen signaling pathway | 交叉调控 | 1.28e-11 | 59.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.28e-11 | 59.0 |
| Bladder cancer | 交叉调控 | 4.13e-11 | 150.0 |
| Prostate cancer | 交叉调控 | 7.12e-11 | 74.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.98e-11 | 71.8 |
| MicroRNAs in cancer | 交叉调控 | 1.04e-10 | 29.8 |
| HIF-1 signaling pathway | 信号主干通路 | 1.23e-10 | 65.8 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 1.81e-10 | 27.6 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.31e-10 | 36.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.88e-10 | 26.1 |
| Human cytomegalovirus infection | 其他 | 3.12e-10 | 36.4 |
| positive regulation of cellular biosynthetic process (G | 其他 | 7.61e-10 | 45.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.4_K30</b> — 显著40条, 主干15, 交叉12, 直接OX0, EF=44.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.88e-25 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 7.88e-25 | 56.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.64e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 8.64e-17 | 41.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.07e-15 | 22.7 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.07e-15 | 18.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.31e-15 | 42.5 |
| Hepatitis B | 交叉调控 | 1.06e-14 | 46.4 |
| Chemical carcinogenesis | 交叉调控 | 1.25e-14 | 34.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.93e-14 | 25.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.54e-14 | 23.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.07e-13 | 49.2 |
| MicroRNAs in cancer | 交叉调控 | 2.00e-13 | 26.5 |
| Human cytomegalovirus infection | 其他 | 2.25e-13 | 33.4 |
| HIF-1 signaling pathway | 信号主干通路 | 6.27e-13 | 56.4 |
| Thyroid hormone signaling pathway | 其他 | 1.52e-12 | 50.8 |
| Pancreatic cancer | 信号主干通路 | 2.22e-12 | 71.9 |
| Estrogen signaling pathway | 交叉调控 | 3.92e-12 | 44.9 |
| Rap1 signaling pathway | 信号主干通路 | 3.92e-12 | 32.5 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 6.76e-12 | 61.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.4_K50</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=32.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.36e-32 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 3.10e-27 | 42.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.54e-22 | 24.3 |
| MAPK signaling pathway | 信号主干通路 | 3.65e-19 | 25.1 |
| Hepatitis B | 交叉调控 | 1.56e-18 | 38.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.59e-18 | 20.3 |
| Prostate cancer | 交叉调控 | 3.69e-18 | 54.9 |
| Central carbon metabolism in cancer | 交叉调控 | 3.96e-18 | 70.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.43e-17 | 31.9 |
| MicroRNAs in cancer | 交叉调控 | 1.55e-17 | 22.5 |
| Rap1 signaling pathway | 信号主干通路 | 4.15e-17 | 29.3 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 1.04e-16 | 16.5 |
| Ras signaling pathway | 信号主干通路 | 1.70e-16 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.17e-16 | 49.2 |
| Lipid and atherosclerosis | 交叉调控 | 1.87e-15 | 26.7 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 8.86e-15 | 21.7 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.10e-14 | 46.4 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.13e-14 | 24.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.34e-14 | 15.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.66e-14 | 27.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.5_K10</b> — 显著40条, 主干12, 交叉13, 直接OX0, EF=90.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.57e-13 | 80.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.23e-10 | 66.7 |
| Pathways in cancer | 交叉调控 | 4.70e-10 | 30.9 |
| Estrogen signaling pathway | 交叉调控 | 6.65e-10 | 89.8 |
| Hepatitis B | 交叉调控 | 1.47e-09 | 75.9 |
| Chemical carcinogenesis | 交叉调控 | 1.29e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.64e-08 | 26.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.39e-08 | 73.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.91e-07 | 29.8 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 6.42e-07 | 210.3 |
| Human cytomegalovirus infection | 其他 | 7.26e-07 | 45.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.45e-06 | 63.7 |
| IL-17 signaling pathway | 交叉调控 | 1.46e-06 | 87.2 |
| Prostate cancer | 交叉调控 | 1.49e-06 | 84.5 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 1.54e-06 | 19.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.54e-06 | 82.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.60e-06 | 30.4 |
| HIF-1 signaling pathway | 信号主干通路 | 1.99e-06 | 75.2 |
| Toxoplasmosis | 其他 | 2.05e-06 | 73.2 |
| MicroRNAs in cancer | 交叉调控 | 2.05e-06 | 33.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.5_K20</b> — 显著40条, 主干12, 交叉14, 直接OX0, EF=52.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.72e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 2.24e-18 | 29.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.55e-16 | 58.4 |
| Chemical carcinogenesis | 交叉调控 | 3.72e-15 | 47.2 |
| Lipid and atherosclerosis | 交叉调控 | 9.72e-14 | 47.7 |
| Hepatitis B | 交叉调控 | 5.38e-13 | 56.9 |
| Prostate cancer | 交叉调控 | 6.76e-13 | 84.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.61e-13 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.37e-12 | 75.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.56e-12 | 21.5 |
| MicroRNAs in cancer | 交叉调控 | 1.93e-12 | 33.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.98e-12 | 25.5 |
| Human cytomegalovirus infection | 其他 | 5.83e-12 | 41.0 |
| Estrogen signaling pathway | 交叉调控 | 6.68e-12 | 59.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.94e-12 | 59.0 |
| Pancreatic cancer | 信号主干通路 | 7.93e-12 | 94.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.06e-11 | 27.9 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 2.32e-11 | 80.6 |
| TNF signaling pathway | 交叉调控 | 1.13e-10 | 64.1 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 1.70e-10 | 27.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.5_K30</b> — 显著40条, 主干14, 交叉13, 直接OX0, EF=43.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.79e-25 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 7.79e-25 | 56.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.23e-19 | 82.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.06e-18 | 44.5 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.64e-17 | 46.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.83e-17 | 24.1 |
| Hepatitis B | 交叉调控 | 1.27e-16 | 50.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.03e-15 | 18.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.61e-15 | 54.1 |
| Chemical carcinogenesis | 交叉调控 | 1.09e-14 | 34.3 |
| Pancreatic cancer | 信号主干通路 | 2.57e-14 | 80.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.57e-14 | 25.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.27e-14 | 23.7 |
| Human immunodeficiency virus 1 infection | 其他 | 1.15e-13 | 35.5 |
| Human cytomegalovirus infection | 其他 | 2.04e-13 | 33.4 |
| C-type lectin receptor signaling pathway | 交叉调控 | 3.71e-13 | 59.1 |
| Th17 cell differentiation | 交叉调控 | 4.48e-13 | 57.5 |
| HIF-1 signaling pathway | 信号主干通路 | 4.96e-13 | 56.4 |
| TNF signaling pathway | 交叉调控 | 5.63e-13 | 54.9 |
| Toxoplasmosis | 其他 | 5.63e-13 | 54.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.5_K50</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=33.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 2.26e-34 | 23.9 |
| Proteoglycans in cancer | 交叉调控 | 3.10e-27 | 42.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.54e-22 | 24.3 |
| MAPK signaling pathway | 信号主干通路 | 8.62e-21 | 26.5 |
| Hepatitis B | 交叉调控 | 1.49e-18 | 38.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.49e-18 | 31.2 |
| Prostate cancer | 交叉调控 | 3.16e-18 | 54.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.45e-18 | 20.3 |
| Central carbon metabolism in cancer | 交叉调控 | 3.47e-18 | 70.3 |
| Ras signaling pathway | 信号主干通路 | 4.99e-18 | 28.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.60e-16 | 49.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 3.71e-16 | 23.0 |
| MicroRNAs in cancer | 交叉调控 | 4.26e-16 | 21.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.41e-16 | 29.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 4.75e-16 | 25.9 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.53e-16 | 29.7 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 8.79e-16 | 23.9 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 8.79e-16 | 15.7 |
| Lipid and atherosclerosis | 交叉调控 | 1.87e-15 | 26.7 |
| positive regulation of intracellular signal transductio | 其他 | 4.48e-15 | 14.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.7_K10</b> — 显著40条, 主干8, 交叉10, 直接OX0, EF=96.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.43e-13 | 80.0 |
| Chemical carcinogenesis | 交叉调控 | 2.54e-10 | 60.0 |
| Pathways in cancer | 交叉调控 | 4.60e-10 | 30.9 |
| Lipid and atherosclerosis | 交叉调控 | 9.99e-09 | 57.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.99e-08 | 26.4 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.99e-08 | 76.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.08e-08 | 59.4 |
| Estrogen signaling pathway | 交叉调控 | 9.41e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.00e-07 | 29.8 |
| Hepatitis B | 交叉调控 | 1.82e-07 | 63.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.90e-07 | 26.3 |
| positive regulation of pri-miRNA transcription by RNA p | 其他 | 1.90e-07 | 241.2 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 2.89e-07 | 210.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.77e-07 | 53.1 |
| regulation of pri-miRNA transcription by RNA polymerase | 其他 | 4.57e-07 | 182.2 |
| Human cytomegalovirus infection | 其他 | 7.10e-07 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.07e-06 | 30.4 |
| positive regulation of cellular biosynthetic process (G | 其他 | 1.07e-06 | 56.9 |
| negative regulation of production of miRNAs involved in | 其他 | 1.41e-06 | 512.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.84e-06 | 82.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.7_K20</b> — 显著40条, 主干10, 交叉16, 直接OX1, EF=61.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.72e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 1.85e-20 | 30.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.76e-17 | 102.5 |
| Hepatitis B | 交叉调控 | 4.76e-17 | 69.6 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.73e-16 | 58.4 |
| Lipid and atherosclerosis | 交叉调控 | 7.62e-16 | 52.4 |
| Prolactin signaling pathway | 交叉调控 | 4.48e-14 | 117.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 5.61e-14 | 27.6 |
| Estrogen signaling pathway | 交叉调控 | 8.66e-14 | 67.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.80e-14 | 66.4 |
| Chemical carcinogenesis | 交叉调控 | 1.41e-13 | 42.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.29e-13 | 30.4 |
| IL-17 signaling pathway | 交叉调控 | 3.32e-13 | 87.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 4.99e-13 | 21.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 7.01e-13 | 78.8 |
| TNF signaling pathway | 交叉调控 | 1.10e-12 | 73.2 |
| Toxoplasmosis | 其他 | 1.10e-12 | 73.2 |
| Relaxin signaling pathway | 交叉调控 | 3.27e-12 | 63.6 |
| Human cytomegalovirus infection | 其他 | 4.01e-12 | 41.0 |
| Pancreatic cancer | 信号主干通路 | 6.53e-12 | 94.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.7_K30</b> — 显著40条, 主干15, 交叉13, 直接OX2, EF=38.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.11e-24 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 9.39e-23 | 53.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.10e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 8.10e-17 | 41.3 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 8.70e-17 | 27.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 7.30e-16 | 22.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.23e-15 | 42.5 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.78e-15 | 54.1 |
| MAPK signaling pathway | 信号主干通路 | 2.81e-15 | 30.2 |
| Hepatitis B | 交叉调控 | 7.47e-15 | 46.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.70e-14 | 17.6 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.09e-13 | 32.4 |
| TNF signaling pathway | 交叉调控 | 8.30e-13 | 54.9 |
| Toxoplasmosis | 其他 | 8.30e-13 | 54.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 8.30e-13 | 23.2 |
| Prolactin signaling pathway | 交叉调控 | 1.22e-12 | 78.1 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.16e-12 | 27.1 |
| Pancreatic cancer | 信号主干通路 | 2.24e-12 | 71.9 |
| positive regulation of intracellular signal transductio | 其他 | 2.38e-12 | 17.5 |
| Focal adhesion | 交叉调控 | 2.66e-12 | 34.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.7_K50</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=37.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.97e-40 | 26.3 |
| MAPK signaling pathway | 信号主干通路 | 1.13e-25 | 30.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.13e-25 | 26.6 |
| Proteoglycans in cancer | 交叉调控 | 1.13e-25 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 9.92e-24 | 37.1 |
| Ras signaling pathway | 信号主干通路 | 5.77e-23 | 33.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.91e-23 | 23.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 7.53e-21 | 30.8 |
| Focal adhesion | 交叉调控 | 6.07e-19 | 32.6 |
| Hepatitis B | 交叉调控 | 9.49e-19 | 38.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.41e-18 | 30.5 |
| Prostate cancer | 交叉调控 | 2.15e-18 | 54.9 |
| Central carbon metabolism in cancer | 交叉调控 | 2.45e-18 | 70.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.72e-18 | 53.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 6.62e-18 | 24.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 8.57e-18 | 31.9 |
| MicroRNAs in cancer | 交叉调控 | 9.68e-18 | 22.5 |
| Prolactin signaling pathway | 交叉调控 | 1.87e-16 | 64.4 |
| Pancreatic cancer | 信号主干通路 | 4.58e-16 | 59.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 5.25e-16 | 29.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.9_K10</b> — 显著40条, 主干9, 交叉12, 直接OX0, EF=102.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Lipid and atherosclerosis | 交叉调控 | 1.08e-12 | 76.3 |
| Chemical carcinogenesis | 交叉调控 | 1.27e-12 | 68.6 |
| Pathways in cancer | 交叉调控 | 2.79e-12 | 34.7 |
| Estrogen signaling pathway | 交叉调控 | 2.79e-12 | 104.7 |
| Hepatitis B | 交叉调控 | 7.36e-12 | 88.6 |
| Prolactin signaling pathway | 交叉调控 | 7.97e-12 | 175.7 |
| Colorectal cancer | 信号主干通路 | 2.44e-11 | 143.0 |
| Proteoglycans in cancer | 交叉调控 | 2.44e-11 | 70.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.81e-11 | 123.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.37e-09 | 63.7 |
| Pancreatic cancer | 信号主干通路 | 2.44e-09 | 134.9 |
| Focal adhesion | 交叉调控 | 2.52e-09 | 61.2 |
| ErbB signaling pathway | 其他 | 3.66e-09 | 120.6 |
| Human cytomegalovirus infection | 其他 | 4.27e-09 | 54.7 |
| Prostate cancer | 交叉调控 | 6.22e-09 | 105.7 |
| Shigellosis | 其他 | 6.40e-09 | 50.0 |
| Salmonella infection | 其他 | 6.48e-09 | 49.4 |
| HIF-1 signaling pathway | 信号主干通路 | 9.37e-09 | 94.0 |
| Toxoplasmosis | 其他 | 1.02e-08 | 91.5 |
| Relaxin signaling pathway | 交叉调控 | 1.94e-08 | 79.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.9_K20</b> — 显著40条, 主干13, 交叉11, 直接OX1, EF=55.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 4.40e-18 | 29.0 |
| Proteoglycans in cancer | 交叉调控 | 9.68e-18 | 60.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.16e-17 | 57.2 |
| Hepatitis B | 交叉调控 | 6.72e-15 | 63.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.37e-15 | 92.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.66e-14 | 53.1 |
| Prolactin signaling pathway | 交叉调控 | 4.40e-14 | 117.1 |
| Human immunodeficiency virus 1 infection | 其他 | 5.18e-14 | 48.3 |
| Chemical carcinogenesis | 交叉调控 | 1.54e-13 | 42.9 |
| Toxoplasmosis | 其他 | 1.52e-12 | 73.2 |
| Focal adhesion | 交叉调控 | 2.06e-12 | 45.9 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 3.82e-12 | 25.4 |
| Relaxin signaling pathway | 交叉调控 | 4.02e-12 | 63.6 |
| Ras signaling pathway | 信号主干通路 | 6.33e-12 | 39.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.33e-12 | 59.0 |
| Pancreatic cancer | 信号主干通路 | 7.28e-12 | 94.4 |
| Shigellosis | 其他 | 8.83e-12 | 37.5 |
| ErbB signaling pathway | 其他 | 1.44e-11 | 84.4 |
| Colorectal cancer | 信号主干通路 | 1.48e-11 | 83.4 |
| GnRH signaling pathway | 其他 | 2.47e-11 | 77.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.9_K30</b> — 显著40条, 主干13, 交叉12, 直接OX3, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.14e-24 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 1.04e-18 | 46.7 |
| Lipid and atherosclerosis | 交叉调控 | 1.37e-18 | 44.5 |
| Hepatitis B | 交叉调控 | 1.92e-18 | 54.8 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.12e-17 | 75.2 |
| Prolactin signaling pathway | 交叉调控 | 1.22e-16 | 97.6 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 1.61e-15 | 22.5 |
| MAPK signaling pathway | 信号主干通路 | 2.87e-15 | 30.2 |
| Toxoplasmosis | 其他 | 1.26e-14 | 61.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.80e-14 | 38.9 |
| Focal adhesion | 交叉调控 | 6.52e-14 | 37.4 |
| Epstein-Barr virus infection | 其他 | 6.52e-14 | 37.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.70e-14 | 49.2 |
| Human immunodeficiency virus 1 infection | 其他 | 9.42e-14 | 35.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.11e-13 | 23.7 |
| MicroRNAs in cancer | 交叉调控 | 1.37e-13 | 26.5 |
| Choline metabolism in cancer | 信号主干通路 | 1.79e-13 | 62.8 |
| Ras signaling pathway | 信号主干通路 | 2.08e-13 | 32.4 |
| Chemical carcinogenesis | 交叉调控 | 2.71e-13 | 31.5 |
| Thyroid hormone signaling pathway | 其他 | 1.05e-12 | 50.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.03_S0.9_K50</b> — 显著40条, 主干15, 交叉12, 直接OX1, EF=44.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.30e-32 | 23.2 |
| Hepatitis B | 交叉调控 | 8.88e-22 | 43.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.31e-21 | 22.3 |
| Proteoglycans in cancer | 交叉调控 | 3.65e-20 | 34.0 |
| MicroRNAs in cancer | 交叉调控 | 9.05e-19 | 23.8 |
| MAPK signaling pathway | 信号主干通路 | 1.08e-17 | 23.7 |
| Epstein-Barr virus infection | 其他 | 3.64e-17 | 30.4 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.09e-16 | 49.2 |
| Progesterone-mediated oocyte maturation | 其他 | 3.09e-16 | 49.2 |
| Focal adhesion | 交叉调控 | 9.93e-16 | 28.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.66e-15 | 27.3 |
| Lipid and atherosclerosis | 交叉调控 | 2.10e-15 | 26.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.04e-15 | 18.5 |
| Prostate cancer | 交叉调控 | 9.24e-15 | 46.5 |
| Prolactin signaling pathway | 交叉调控 | 1.79e-14 | 58.6 |
| Human papillomavirus infection | 其他 | 2.39e-14 | 18.6 |
| Ras signaling pathway | 信号主干通路 | 1.44e-13 | 23.0 |
| Chemical carcinogenesis | 交叉调控 | 1.99e-13 | 22.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 3.88e-13 | 32.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.14e-13 | 25.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.4_K10</b> — 显著40条, 主干11, 交叉14, 直接OX0, EF=100.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.71e-13 | 80.0 |
| Pathways in cancer | 交叉调控 | 7.20e-10 | 30.9 |
| Hepatitis B | 交叉调控 | 2.51e-09 | 75.9 |
| Lipid and atherosclerosis | 交叉调控 | 1.04e-08 | 57.2 |
| Human cytomegalovirus infection | 其他 | 1.10e-08 | 54.7 |
| Chemical carcinogenesis | 交叉调控 | 1.32e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.55e-08 | 26.4 |
| MicroRNAs in cancer | 交叉调控 | 5.38e-08 | 39.7 |
| Estrogen signaling pathway | 交叉调控 | 6.14e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.85e-07 | 29.8 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.06e-07 | 53.1 |
| IL-17 signaling pathway | 交叉调控 | 1.34e-06 | 87.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.57e-06 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.95e-06 | 75.2 |
| TNF signaling pathway | 交叉调控 | 1.95e-06 | 73.2 |
| Toxoplasmosis | 其他 | 1.95e-06 | 73.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.33e-06 | 59.0 |
| negative regulation of production of miRNAs involved in | 其他 | 4.77e-06 | 512.5 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 6.39e-06 | 25.4 |
| negative regulation of apoptotic signaling pathway (GO: | 信号主干通路 | 6.39e-06 | 105.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.4_K20</b> — 显著40条, 主干17, 交叉13, 直接OX0, EF=43.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 4.46e-18 | 29.0 |
| Proteoglycans in cancer | 交叉调控 | 9.79e-18 | 60.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.61e-13 | 47.7 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.50e-12 | 21.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.90e-12 | 25.5 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.96e-12 | 47.8 |
| Chemical carcinogenesis | 交叉调控 | 1.84e-11 | 38.6 |
| Pancreatic cancer | 信号主干通路 | 1.84e-11 | 94.4 |
| Hepatitis B | 交叉调控 | 4.44e-11 | 50.6 |
| Prostate cancer | 交叉调控 | 8.01e-11 | 74.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.86e-11 | 71.8 |
| MicroRNAs in cancer | 交叉调控 | 1.14e-10 | 29.8 |
| Th17 cell differentiation | 交叉调控 | 1.18e-10 | 67.1 |
| HIF-1 signaling pathway | 信号主干通路 | 1.23e-10 | 65.8 |
| Toxoplasmosis | 其他 | 1.38e-10 | 64.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.68e-10 | 26.1 |
| Human cytomegalovirus infection | 其他 | 2.91e-10 | 36.4 |
| Estrogen signaling pathway | 交叉调控 | 4.69e-10 | 52.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.89e-10 | 51.6 |
| JAK-STAT signaling pathway | 信号主干通路 | 1.36e-09 | 44.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.4_K30</b> — 显著40条, 主干13, 交叉13, 直接OX0, EF=44.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.25e-26 | 28.3 |
| Proteoglycans in cancer | 交叉调控 | 7.75e-25 | 56.7 |
| Hepatitis B | 交叉调控 | 2.40e-20 | 59.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.65e-19 | 82.0 |
| Lipid and atherosclerosis | 交叉调控 | 8.43e-19 | 44.5 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.36e-17 | 46.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-15 | 22.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.60e-15 | 54.1 |
| Chemical carcinogenesis | 交叉调控 | 1.08e-14 | 34.3 |
| Pancreatic cancer | 信号主干通路 | 2.55e-14 | 80.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.55e-14 | 25.1 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 3.92e-14 | 17.6 |
| Epstein-Barr virus infection | 其他 | 6.69e-14 | 37.2 |
| Human immunodeficiency virus 1 infection | 其他 | 1.05e-13 | 35.5 |
| Human cytomegalovirus infection | 其他 | 1.79e-13 | 33.4 |
| Prostate cancer | 交叉调控 | 1.79e-13 | 63.4 |
| C-type lectin receptor signaling pathway | 交叉调控 | 3.20e-13 | 59.1 |
| Th17 cell differentiation | 交叉调控 | 3.90e-13 | 57.5 |
| HIF-1 signaling pathway | 信号主干通路 | 4.35e-13 | 56.4 |
| Toxoplasmosis | 其他 | 5.01e-13 | 54.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.4_K50</b> — 显著40条, 主干17, 交叉13, 直接OX0, EF=31.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 2.24e-34 | 23.9 |
| Proteoglycans in cancer | 交叉调控 | 3.07e-27 | 42.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.04e-20 | 23.2 |
| Hepatitis B | 交叉调控 | 3.14e-20 | 40.5 |
| MAPK signaling pathway | 信号主干通路 | 1.12e-17 | 23.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.89e-17 | 31.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.18e-16 | 19.3 |
| Prostate cancer | 交叉调控 | 2.51e-16 | 50.7 |
| Central carbon metabolism in cancer | 交叉调控 | 3.17e-16 | 64.4 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.17e-16 | 49.2 |
| MicroRNAs in cancer | 交叉调控 | 4.64e-16 | 21.2 |
| Epstein-Barr virus infection | 其他 | 9.08e-16 | 28.4 |
| Rap1 signaling pathway | 信号主干通路 | 1.44e-15 | 27.3 |
| Lipid and atherosclerosis | 交叉调控 | 1.85e-15 | 26.7 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.76e-15 | 15.7 |
| Ras signaling pathway | 信号主干通路 | 5.01e-15 | 24.7 |
| Chemical carcinogenesis | 交叉调控 | 7.09e-15 | 24.0 |
| Prolactin signaling pathway | 交叉调控 | 1.63e-14 | 58.6 |
| Gastric cancer | 交叉调控 | 2.10e-14 | 33.0 |
| HIF-1 signaling pathway | 信号主干通路 | 2.62e-14 | 41.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.5_K10</b> — 显著40条, 主干12, 交叉13, 直接OX0, EF=90.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.57e-13 | 80.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.23e-10 | 66.7 |
| Pathways in cancer | 交叉调控 | 4.70e-10 | 30.9 |
| Estrogen signaling pathway | 交叉调控 | 6.65e-10 | 89.8 |
| Hepatitis B | 交叉调控 | 1.47e-09 | 75.9 |
| Chemical carcinogenesis | 交叉调控 | 1.29e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.64e-08 | 26.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.39e-08 | 73.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.91e-07 | 29.8 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 6.42e-07 | 210.3 |
| Human cytomegalovirus infection | 其他 | 7.26e-07 | 45.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.45e-06 | 63.7 |
| IL-17 signaling pathway | 交叉调控 | 1.46e-06 | 87.2 |
| Prostate cancer | 交叉调控 | 1.49e-06 | 84.5 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 1.54e-06 | 19.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.54e-06 | 82.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.60e-06 | 30.4 |
| HIF-1 signaling pathway | 信号主干通路 | 1.99e-06 | 75.2 |
| Toxoplasmosis | 其他 | 2.05e-06 | 73.2 |
| MicroRNAs in cancer | 交叉调控 | 2.05e-06 | 33.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.5_K20</b> — 显著40条, 主干10, 交叉14, 直接OX0, EF=55.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 2.21e-24 | 75.0 |
| Pathways in cancer | 交叉调控 | 2.11e-16 | 27.0 |
| Chemical carcinogenesis | 交叉调控 | 4.96e-15 | 47.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.07e-14 | 53.1 |
| Lipid and atherosclerosis | 交叉调控 | 9.72e-14 | 47.7 |
| Human cytomegalovirus infection | 其他 | 1.28e-13 | 45.6 |
| Hepatitis B | 交叉调控 | 4.61e-13 | 56.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.92e-13 | 30.4 |
| Prostate cancer | 交叉调控 | 5.92e-13 | 84.5 |
| HIF-1 signaling pathway | 信号主干通路 | 1.37e-12 | 75.2 |
| MicroRNAs in cancer | 交叉调控 | 1.93e-12 | 33.1 |
| Estrogen signaling pathway | 交叉调控 | 7.29e-12 | 59.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.51e-12 | 59.0 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 2.67e-11 | 80.6 |
| IL-17 signaling pathway | 交叉调控 | 3.68e-11 | 76.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.35e-11 | 71.8 |
| TNF signaling pathway | 交叉调控 | 1.13e-10 | 64.1 |
| VEGF signaling pathway | 其他 | 1.94e-10 | 104.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.09e-10 | 26.1 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 2.93e-10 | 27.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.5_K30</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=39.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.19e-26 | 28.3 |
| Proteoglycans in cancer | 交叉调控 | 7.38e-25 | 56.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.59e-17 | 46.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.85e-17 | 24.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.32e-17 | 27.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.26e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 6.48e-17 | 41.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 6.89e-16 | 18.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.70e-16 | 27.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.52e-15 | 54.1 |
| Hepatitis B | 交叉调控 | 7.47e-15 | 46.4 |
| Central carbon metabolism in cancer | 交叉调控 | 1.24e-14 | 87.9 |
| Pancreatic cancer | 信号主干通路 | 2.43e-14 | 80.9 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 5.47e-14 | 29.3 |
| positive regulation of intracellular signal transductio | 其他 | 8.25e-14 | 18.8 |
| MAPK signaling pathway | 信号主干通路 | 9.04e-14 | 27.9 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.05e-13 | 32.4 |
| Chemical carcinogenesis | 交叉调控 | 3.51e-13 | 31.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 3.51e-13 | 59.1 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 3.96e-13 | 19.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.5_K50</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=29.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.43e-36 | 24.7 |
| Proteoglycans in cancer | 交叉调控 | 3.09e-27 | 42.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 8.57e-21 | 23.2 |
| MAPK signaling pathway | 信号主干通路 | 8.57e-21 | 26.5 |
| Hepatitis B | 交叉调控 | 2.52e-20 | 40.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.58e-18 | 53.3 |
| Ras signaling pathway | 信号主干通路 | 6.39e-18 | 28.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.43e-17 | 31.9 |
| Rap1 signaling pathway | 信号主干通路 | 4.58e-17 | 29.3 |
| Lipid and atherosclerosis | 交叉调控 | 5.89e-17 | 28.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.17e-16 | 19.3 |
| Prostate cancer | 交叉调控 | 1.60e-16 | 50.7 |
| Central carbon metabolism in cancer | 交叉调控 | 2.21e-16 | 64.4 |
| Prolactin signaling pathway | 交叉调控 | 2.21e-16 | 64.4 |
| MicroRNAs in cancer | 交叉调控 | 3.33e-16 | 21.2 |
| Pancreatic cancer | 信号主干通路 | 5.01e-16 | 59.3 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.72e-15 | 15.7 |
| Chemical carcinogenesis | 交叉调控 | 6.68e-15 | 24.0 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 8.77e-15 | 21.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.01e-15 | 35.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.7_K10</b> — 显著40条, 主干12, 交叉11, 直接OX0, EF=98.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 7.42e-13 | 80.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.39e-10 | 66.7 |
| Chemical carcinogenesis | 交叉调控 | 1.95e-10 | 60.0 |
| Pathways in cancer | 交叉调控 | 3.98e-10 | 30.9 |
| Estrogen signaling pathway | 交叉调控 | 6.00e-10 | 89.8 |
| Hepatitis B | 交叉调控 | 1.39e-09 | 75.9 |
| Human cytomegalovirus infection | 其他 | 8.66e-09 | 54.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.38e-08 | 102.5 |
| HIF-1 signaling pathway | 信号主干通路 | 1.90e-08 | 94.0 |
| Toxoplasmosis | 其他 | 1.96e-08 | 91.5 |
| Hepatitis C | 信号主干通路 | 9.80e-08 | 65.3 |
| Tuberculosis | 其他 | 1.79e-07 | 56.9 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.34e-07 | 53.1 |
| Prolactin signaling pathway | 交叉调控 | 3.21e-07 | 117.1 |
| Shigellosis | 其他 | 6.40e-07 | 41.7 |
| Colorectal cancer | 信号主干通路 | 6.40e-07 | 95.3 |
| Salmonella infection | 其他 | 6.40e-07 | 41.2 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 6.62e-07 | 92.1 |
| IL-17 signaling pathway | 交叉调控 | 7.82e-07 | 87.2 |
| Prostate cancer | 交叉调控 | 8.43e-07 | 84.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.7_K20</b> — 显著40条, 主干10, 交叉16, 直接OX1, EF=61.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.72e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 1.85e-20 | 30.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.76e-17 | 102.5 |
| Hepatitis B | 交叉调控 | 4.76e-17 | 69.6 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.73e-16 | 58.4 |
| Lipid and atherosclerosis | 交叉调控 | 7.62e-16 | 52.4 |
| Prolactin signaling pathway | 交叉调控 | 4.48e-14 | 117.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 5.61e-14 | 27.6 |
| Estrogen signaling pathway | 交叉调控 | 8.66e-14 | 67.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.80e-14 | 66.4 |
| Chemical carcinogenesis | 交叉调控 | 1.41e-13 | 42.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.29e-13 | 30.4 |
| IL-17 signaling pathway | 交叉调控 | 3.32e-13 | 87.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 4.99e-13 | 21.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 7.01e-13 | 78.8 |
| TNF signaling pathway | 交叉调控 | 1.10e-12 | 73.2 |
| Toxoplasmosis | 其他 | 1.10e-12 | 73.2 |
| Relaxin signaling pathway | 交叉调控 | 3.27e-12 | 63.6 |
| Human cytomegalovirus infection | 其他 | 4.01e-12 | 41.0 |
| Pancreatic cancer | 信号主干通路 | 6.53e-12 | 94.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.7_K30</b> — 显著40条, 主干14, 交叉13, 直接OX1, EF=38.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.38e-25 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 7.38e-25 | 56.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.59e-17 | 46.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.26e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 6.48e-17 | 41.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.50e-15 | 22.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.78e-15 | 54.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.40e-15 | 25.4 |
| Hepatitis B | 交叉调控 | 8.54e-15 | 46.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.78e-14 | 17.6 |
| Pancreatic cancer | 信号主干通路 | 3.04e-14 | 80.9 |
| MAPK signaling pathway | 信号主干通路 | 1.11e-13 | 27.9 |
| MicroRNAs in cancer | 交叉调控 | 1.88e-13 | 26.5 |
| Chemical carcinogenesis | 交叉调控 | 4.10e-13 | 31.5 |
| TNF signaling pathway | 交叉调控 | 6.52e-13 | 54.9 |
| Toxoplasmosis | 其他 | 6.52e-13 | 54.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.52e-13 | 23.2 |
| Central carbon metabolism in cancer | 交叉调控 | 9.17e-13 | 78.1 |
| Prolactin signaling pathway | 交叉调控 | 9.17e-13 | 78.1 |
| Epstein-Barr virus infection | 其他 | 2.34e-12 | 33.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.7_K50</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=30.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.91e-40 | 26.3 |
| Proteoglycans in cancer | 交叉调控 | 2.98e-27 | 42.0 |
| MAPK signaling pathway | 信号主干通路 | 1.83e-22 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.83e-22 | 24.3 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.50e-21 | 22.3 |
| Hepatitis B | 交叉调控 | 2.44e-20 | 40.5 |
| Ras signaling pathway | 信号主干通路 | 1.57e-19 | 30.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.80e-18 | 30.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.04e-18 | 53.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.22e-17 | 31.9 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.72e-17 | 27.5 |
| Rap1 signaling pathway | 信号主干通路 | 3.98e-17 | 29.3 |
| Prostate cancer | 交叉调控 | 1.55e-16 | 50.7 |
| Central carbon metabolism in cancer | 交叉调控 | 2.13e-16 | 64.4 |
| Prolactin signaling pathway | 交叉调控 | 2.13e-16 | 64.4 |
| MicroRNAs in cancer | 交叉调控 | 3.21e-16 | 21.2 |
| Pancreatic cancer | 信号主干通路 | 4.83e-16 | 59.3 |
| Focal adhesion | 交叉调控 | 5.65e-16 | 28.6 |
| Epstein-Barr virus infection | 其他 | 5.70e-16 | 28.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.04e-15 | 16.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.9_K10</b> — 显著40条, 主干9, 交叉12, 直接OX0, EF=102.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Lipid and atherosclerosis | 交叉调控 | 1.08e-12 | 76.3 |
| Chemical carcinogenesis | 交叉调控 | 1.27e-12 | 68.6 |
| Pathways in cancer | 交叉调控 | 2.79e-12 | 34.7 |
| Estrogen signaling pathway | 交叉调控 | 2.79e-12 | 104.7 |
| Hepatitis B | 交叉调控 | 7.36e-12 | 88.6 |
| Prolactin signaling pathway | 交叉调控 | 7.97e-12 | 175.7 |
| Colorectal cancer | 信号主干通路 | 2.44e-11 | 143.0 |
| Proteoglycans in cancer | 交叉调控 | 2.44e-11 | 70.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.81e-11 | 123.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.37e-09 | 63.7 |
| Pancreatic cancer | 信号主干通路 | 2.44e-09 | 134.9 |
| Focal adhesion | 交叉调控 | 2.52e-09 | 61.2 |
| ErbB signaling pathway | 其他 | 3.66e-09 | 120.6 |
| Human cytomegalovirus infection | 其他 | 4.27e-09 | 54.7 |
| Prostate cancer | 交叉调控 | 6.22e-09 | 105.7 |
| Shigellosis | 其他 | 6.40e-09 | 50.0 |
| Salmonella infection | 其他 | 6.48e-09 | 49.4 |
| HIF-1 signaling pathway | 信号主干通路 | 9.37e-09 | 94.0 |
| Toxoplasmosis | 其他 | 1.02e-08 | 91.5 |
| Relaxin signaling pathway | 交叉调控 | 1.94e-08 | 79.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.9_K20</b> — 显著40条, 主干10, 交叉13, 直接OX1, EF=57.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 4.40e-18 | 29.0 |
| Proteoglycans in cancer | 交叉调控 | 9.68e-18 | 60.0 |
| Prolactin signaling pathway | 交叉调控 | 4.26e-16 | 131.8 |
| Lipid and atherosclerosis | 交叉调控 | 1.12e-15 | 52.4 |
| Chemical carcinogenesis | 交叉调控 | 2.92e-15 | 47.2 |
| Hepatitis B | 交叉调控 | 4.48e-15 | 63.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.26e-15 | 92.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.00e-14 | 53.1 |
| Focal adhesion | 交叉调控 | 2.52e-12 | 45.9 |
| Human immunodeficiency virus 1 infection | 其他 | 3.68e-12 | 43.5 |
| Relaxin signaling pathway | 交叉调控 | 4.39e-12 | 63.6 |
| Estrogen signaling pathway | 交叉调控 | 6.56e-12 | 59.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.82e-12 | 59.0 |
| ErbB signaling pathway | 其他 | 1.75e-11 | 84.4 |
| Colorectal cancer | 信号主干通路 | 1.78e-11 | 83.4 |
| GnRH signaling pathway | 其他 | 2.93e-11 | 77.2 |
| MAPK signaling pathway | 信号主干通路 | 3.90e-11 | 31.4 |
| Choline metabolism in cancer | 信号主干通路 | 3.90e-11 | 73.2 |
| Progesterone-mediated oocyte maturation | 其他 | 4.15e-11 | 71.8 |
| Th17 cell differentiation | 交叉调控 | 6.40e-11 | 67.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.9_K30</b> — 显著40条, 主干10, 交叉11, 直接OX1, EF=41.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 9.08e-23 | 25.7 |
| Hepatitis B | 交叉调控 | 3.49e-20 | 59.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.16e-19 | 82.0 |
| Proteoglycans in cancer | 交叉调控 | 5.16e-19 | 46.7 |
| Progesterone-mediated oocyte maturation | 其他 | 5.09e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 5.49e-17 | 41.3 |
| Prolactin signaling pathway | 交叉调控 | 1.04e-16 | 97.6 |
| Epstein-Barr virus infection | 其他 | 1.36e-15 | 40.6 |
| Chemical carcinogenesis | 交叉调控 | 9.29e-15 | 34.3 |
| Toxoplasmosis | 其他 | 1.00e-14 | 61.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.91e-14 | 38.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.66e-14 | 49.2 |
| MAPK signaling pathway | 信号主干通路 | 7.78e-14 | 27.9 |
| Human immunodeficiency virus 1 infection | 其他 | 8.70e-14 | 35.5 |
| Pancreatic cancer | 信号主干通路 | 1.98e-12 | 71.9 |
| Relaxin signaling pathway | 交叉调控 | 2.11e-12 | 47.7 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 2.58e-12 | 19.7 |
| Estrogen signaling pathway | 交叉调控 | 3.45e-12 | 44.9 |
| Coronavirus disease | 其他 | 8.53e-12 | 29.5 |
| GnRH signaling pathway | 其他 | 8.53e-12 | 58.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.05_S0.9_K50</b> — 显著40条, 主干13, 交叉12, 直接OX0, EF=33.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.18e-36 | 24.7 |
| Proteoglycans in cancer | 交叉调控 | 2.14e-25 | 40.0 |
| Hepatitis B | 交叉调控 | 7.68e-24 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.47e-23 | 23.3 |
| MAPK signaling pathway | 信号主干通路 | 7.94e-21 | 26.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.19e-19 | 22.0 |
| MicroRNAs in cancer | 交叉调控 | 5.85e-19 | 23.8 |
| Progesterone-mediated oocyte maturation | 其他 | 4.42e-18 | 53.3 |
| Epstein-Barr virus infection | 其他 | 2.65e-17 | 30.4 |
| Rap1 signaling pathway | 信号主干通路 | 4.24e-17 | 29.3 |
| Prostate cancer | 交叉调控 | 1.57e-16 | 50.7 |
| Ras signaling pathway | 信号主干通路 | 1.57e-16 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.93e-16 | 49.2 |
| Prolactin signaling pathway | 交叉调控 | 1.93e-16 | 64.4 |
| Chemical carcinogenesis | 交叉调控 | 1.93e-16 | 25.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.25e-16 | 29.7 |
| Pancreatic cancer | 信号主干通路 | 4.35e-16 | 59.3 |
| Focal adhesion | 交叉调控 | 5.10e-16 | 28.6 |
| Human papillomavirus infection | 其他 | 6.78e-16 | 19.8 |
| Lipid and atherosclerosis | 交叉调控 | 1.18e-15 | 26.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.4_K10</b> — 显著40条, 主干11, 交叉14, 直接OX0, EF=100.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.71e-13 | 80.0 |
| Pathways in cancer | 交叉调控 | 7.20e-10 | 30.9 |
| Hepatitis B | 交叉调控 | 2.51e-09 | 75.9 |
| Lipid and atherosclerosis | 交叉调控 | 1.04e-08 | 57.2 |
| Human cytomegalovirus infection | 其他 | 1.10e-08 | 54.7 |
| Chemical carcinogenesis | 交叉调控 | 1.32e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.55e-08 | 26.4 |
| MicroRNAs in cancer | 交叉调控 | 5.38e-08 | 39.7 |
| Estrogen signaling pathway | 交叉调控 | 6.14e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.85e-07 | 29.8 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.06e-07 | 53.1 |
| IL-17 signaling pathway | 交叉调控 | 1.34e-06 | 87.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.57e-06 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.95e-06 | 75.2 |
| TNF signaling pathway | 交叉调控 | 1.95e-06 | 73.2 |
| Toxoplasmosis | 其他 | 1.95e-06 | 73.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.33e-06 | 59.0 |
| negative regulation of production of miRNAs involved in | 其他 | 4.77e-06 | 512.5 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 6.39e-06 | 25.4 |
| negative regulation of apoptotic signaling pathway (GO: | 信号主干通路 | 6.39e-06 | 105.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.4_K20</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=40.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.07e-19 | 65.0 |
| Pathways in cancer | 交叉调控 | 2.03e-18 | 29.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.33e-14 | 23.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.79e-14 | 27.6 |
| Lipid and atherosclerosis | 交叉调控 | 1.29e-11 | 42.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.01e-11 | 59.0 |
| Chemical carcinogenesis | 交叉调控 | 2.01e-11 | 38.6 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.94e-10 | 42.5 |
| positive regulation of intracellular signal transductio | 其他 | 2.71e-10 | 20.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.85e-10 | 36.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.33e-10 | 25.4 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.89e-10 | 26.1 |
| positive regulation of cellular biosynthetic process (G | 其他 | 7.51e-10 | 45.6 |
| Estrogen signaling pathway | 交叉调控 | 8.56e-10 | 52.4 |
| positive regulation of protein modification process (GO | 其他 | 2.59e-09 | 38.3 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 4.41e-09 | 24.9 |
| MicroRNAs in cancer | 交叉调控 | 5.69e-09 | 26.5 |
| Prostate cancer | 交叉调控 | 6.08e-09 | 63.4 |
| Bladder cancer | 交叉调控 | 6.08e-09 | 125.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.43e-09 | 61.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.4_K30</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=45.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.14e-28 | 29.6 |
| Proteoglycans in cancer | 交叉调控 | 7.75e-25 | 56.7 |
| Hepatitis B | 交叉调控 | 2.40e-20 | 59.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.65e-19 | 82.0 |
| Lipid and atherosclerosis | 交叉调控 | 8.43e-19 | 44.5 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.36e-17 | 46.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.79e-17 | 24.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.93e-16 | 27.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.01e-15 | 18.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.40e-15 | 54.1 |
| Chemical carcinogenesis | 交叉调控 | 9.60e-15 | 34.3 |
| Pancreatic cancer | 信号主干通路 | 2.55e-14 | 80.9 |
| Epstein-Barr virus infection | 其他 | 6.69e-14 | 37.2 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.16e-14 | 23.7 |
| Human cytomegalovirus infection | 其他 | 1.93e-13 | 33.4 |
| Prostate cancer | 交叉调控 | 1.93e-13 | 63.4 |
| C-type lectin receptor signaling pathway | 交叉调控 | 3.42e-13 | 59.1 |
| Th17 cell differentiation | 交叉调控 | 4.16e-13 | 57.5 |
| HIF-1 signaling pathway | 信号主干通路 | 4.63e-13 | 56.4 |
| Toxoplasmosis | 其他 | 5.29e-13 | 54.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.4_K50</b> — 显著40条, 主干17, 交叉12, 直接OX0, EF=31.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 2.25e-34 | 23.9 |
| Proteoglycans in cancer | 交叉调控 | 3.09e-27 | 42.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.04e-20 | 23.2 |
| Hepatitis B | 交叉调控 | 3.16e-20 | 40.5 |
| MAPK signaling pathway | 信号主干通路 | 1.13e-17 | 23.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.90e-17 | 31.9 |
| Lipid and atherosclerosis | 交叉调控 | 8.42e-17 | 28.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.17e-16 | 19.3 |
| Prostate cancer | 交叉调控 | 2.21e-16 | 50.7 |
| Central carbon metabolism in cancer | 交叉调控 | 2.87e-16 | 64.4 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.87e-16 | 49.2 |
| Rap1 signaling pathway | 信号主干通路 | 1.58e-15 | 27.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.59e-15 | 16.2 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 1.81e-15 | 15.7 |
| Ras signaling pathway | 信号主干通路 | 5.87e-15 | 24.7 |
| Chemical carcinogenesis | 交叉调控 | 8.22e-15 | 24.0 |
| MicroRNAs in cancer | 交叉调控 | 1.02e-14 | 19.8 |
| Prolactin signaling pathway | 交叉调控 | 1.75e-14 | 58.6 |
| Gastric cancer | 交叉调控 | 2.24e-14 | 33.0 |
| Epstein-Barr virus infection | 其他 | 2.34e-14 | 26.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.5_K10</b> — 显著40条, 主干12, 交叉13, 直接OX0, EF=90.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.57e-13 | 80.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.23e-10 | 66.7 |
| Pathways in cancer | 交叉调控 | 4.70e-10 | 30.9 |
| Estrogen signaling pathway | 交叉调控 | 6.65e-10 | 89.8 |
| Hepatitis B | 交叉调控 | 1.47e-09 | 75.9 |
| Chemical carcinogenesis | 交叉调控 | 1.29e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.64e-08 | 26.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.39e-08 | 73.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.91e-07 | 29.8 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 6.42e-07 | 210.3 |
| Human cytomegalovirus infection | 其他 | 7.26e-07 | 45.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.45e-06 | 63.7 |
| IL-17 signaling pathway | 交叉调控 | 1.46e-06 | 87.2 |
| Prostate cancer | 交叉调控 | 1.49e-06 | 84.5 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 1.54e-06 | 19.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.54e-06 | 82.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.60e-06 | 30.4 |
| HIF-1 signaling pathway | 信号主干通路 | 1.99e-06 | 75.2 |
| Toxoplasmosis | 其他 | 2.05e-06 | 73.2 |
| MicroRNAs in cancer | 交叉调控 | 2.05e-06 | 33.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.5_K20</b> — 显著40条, 主干10, 交叉14, 直接OX0, EF=56.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 2.22e-24 | 75.0 |
| Pathways in cancer | 交叉调控 | 2.13e-16 | 27.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.58e-16 | 58.4 |
| Lipid and atherosclerosis | 交叉调控 | 1.15e-15 | 52.4 |
| Chemical carcinogenesis | 交叉调控 | 2.99e-15 | 47.2 |
| Human cytomegalovirus infection | 其他 | 1.14e-13 | 45.6 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.14e-13 | 66.4 |
| Hepatitis B | 交叉调控 | 4.06e-13 | 56.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.78e-13 | 30.4 |
| Prostate cancer | 交叉调控 | 5.29e-13 | 84.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.12e-13 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.13e-12 | 75.2 |
| TNF signaling pathway | 交叉调控 | 1.30e-12 | 73.2 |
| MicroRNAs in cancer | 交叉调控 | 1.49e-12 | 33.1 |
| Estrogen signaling pathway | 交叉调控 | 5.76e-12 | 59.9 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 2.33e-11 | 80.6 |
| IL-17 signaling pathway | 交叉调控 | 3.24e-11 | 76.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 4.06e-11 | 19.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 7.33e-11 | 23.4 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 1.71e-10 | 27.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.5_K30</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=39.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.19e-26 | 28.3 |
| Proteoglycans in cancer | 交叉调控 | 7.38e-25 | 56.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.59e-17 | 46.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.85e-17 | 24.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.32e-17 | 27.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.26e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 6.48e-17 | 41.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 6.89e-16 | 18.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.70e-16 | 27.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.52e-15 | 54.1 |
| Hepatitis B | 交叉调控 | 7.47e-15 | 46.4 |
| Central carbon metabolism in cancer | 交叉调控 | 1.24e-14 | 87.9 |
| Pancreatic cancer | 信号主干通路 | 2.43e-14 | 80.9 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 5.47e-14 | 29.3 |
| positive regulation of intracellular signal transductio | 其他 | 8.25e-14 | 18.8 |
| MAPK signaling pathway | 信号主干通路 | 9.04e-14 | 27.9 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.05e-13 | 32.4 |
| Chemical carcinogenesis | 交叉调控 | 3.51e-13 | 31.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 3.51e-13 | 59.1 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 3.96e-13 | 19.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.5_K50</b> — 显著40条, 主干14, 交叉13, 直接OX0, EF=32.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 6.10e-40 | 26.3 |
| Proteoglycans in cancer | 交叉调控 | 3.07e-27 | 42.0 |
| Hepatitis B | 交叉调控 | 8.26e-24 | 45.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.89e-22 | 24.3 |
| MAPK signaling pathway | 信号主干通路 | 6.83e-21 | 26.5 |
| Lipid and atherosclerosis | 交叉调控 | 2.16e-18 | 30.5 |
| Prostate cancer | 交叉调控 | 3.13e-18 | 54.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.16e-18 | 53.3 |
| Ras signaling pathway | 信号主干通路 | 4.94e-18 | 28.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.14e-17 | 31.9 |
| Gastric cancer | 交叉调控 | 1.18e-17 | 38.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.75e-17 | 19.3 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 9.75e-17 | 16.5 |
| Central carbon metabolism in cancer | 交叉调控 | 2.20e-16 | 64.4 |
| Prolactin signaling pathway | 交叉调控 | 2.20e-16 | 64.4 |
| MicroRNAs in cancer | 交叉调控 | 3.31e-16 | 21.2 |
| Pancreatic cancer | 信号主干通路 | 4.98e-16 | 59.3 |
| Epstein-Barr virus infection | 其他 | 6.24e-16 | 28.4 |
| Rap1 signaling pathway | 信号主干通路 | 1.02e-15 | 27.3 |
| Human T-cell leukemia virus 1 infection | 其他 | 1.73e-15 | 26.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.7_K10</b> — 显著40条, 主干8, 交叉10, 直接OX0, EF=96.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.43e-13 | 80.0 |
| Chemical carcinogenesis | 交叉调控 | 2.54e-10 | 60.0 |
| Pathways in cancer | 交叉调控 | 4.60e-10 | 30.9 |
| Lipid and atherosclerosis | 交叉调控 | 9.99e-09 | 57.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.99e-08 | 26.4 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.99e-08 | 76.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.08e-08 | 59.4 |
| Estrogen signaling pathway | 交叉调控 | 9.41e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.00e-07 | 29.8 |
| Hepatitis B | 交叉调控 | 1.82e-07 | 63.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.90e-07 | 26.3 |
| positive regulation of pri-miRNA transcription by RNA p | 其他 | 1.90e-07 | 241.2 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 2.89e-07 | 210.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.77e-07 | 53.1 |
| regulation of pri-miRNA transcription by RNA polymerase | 其他 | 4.57e-07 | 182.2 |
| Human cytomegalovirus infection | 其他 | 7.10e-07 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.07e-06 | 30.4 |
| positive regulation of cellular biosynthetic process (G | 其他 | 1.07e-06 | 56.9 |
| negative regulation of production of miRNAs involved in | 其他 | 1.41e-06 | 512.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.84e-06 | 82.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.7_K20</b> — 显著40条, 主干10, 交叉16, 直接OX1, EF=61.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.72e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 1.85e-20 | 30.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.76e-17 | 102.5 |
| Hepatitis B | 交叉调控 | 4.76e-17 | 69.6 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.73e-16 | 58.4 |
| Lipid and atherosclerosis | 交叉调控 | 7.62e-16 | 52.4 |
| Prolactin signaling pathway | 交叉调控 | 4.48e-14 | 117.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 5.61e-14 | 27.6 |
| Estrogen signaling pathway | 交叉调控 | 8.66e-14 | 67.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.80e-14 | 66.4 |
| Chemical carcinogenesis | 交叉调控 | 1.41e-13 | 42.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.29e-13 | 30.4 |
| IL-17 signaling pathway | 交叉调控 | 3.32e-13 | 87.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 4.99e-13 | 21.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 7.01e-13 | 78.8 |
| TNF signaling pathway | 交叉调控 | 1.10e-12 | 73.2 |
| Toxoplasmosis | 其他 | 1.10e-12 | 73.2 |
| Relaxin signaling pathway | 交叉调控 | 3.27e-12 | 63.6 |
| Human cytomegalovirus infection | 其他 | 4.01e-12 | 41.0 |
| Pancreatic cancer | 信号主干通路 | 6.53e-12 | 94.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.7_K30</b> — 显著40条, 主干13, 交叉15, 直接OX2, EF=39.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.38e-25 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 7.38e-25 | 56.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.59e-17 | 46.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.93e-17 | 24.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.26e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 6.48e-17 | 41.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.06e-15 | 18.7 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.59e-15 | 25.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.78e-15 | 54.1 |
| Hepatitis B | 交叉调控 | 8.54e-15 | 46.4 |
| Pancreatic cancer | 信号主干通路 | 2.70e-14 | 80.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.70e-14 | 25.1 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 5.63e-14 | 29.3 |
| MAPK signaling pathway | 信号主干通路 | 9.95e-14 | 27.9 |
| Chemical carcinogenesis | 交叉调控 | 3.80e-13 | 31.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 3.80e-13 | 59.1 |
| Th17 cell differentiation | 交叉调控 | 4.57e-13 | 57.5 |
| TNF signaling pathway | 交叉调控 | 6.04e-13 | 54.9 |
| Toxoplasmosis | 其他 | 6.04e-13 | 54.9 |
| Central carbon metabolism in cancer | 交叉调控 | 8.63e-13 | 78.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.7_K50</b> — 显著40条, 主干14, 交叉13, 直接OX0, EF=31.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.85e-40 | 26.3 |
| Proteoglycans in cancer | 交叉调控 | 2.95e-27 | 42.0 |
| MAPK signaling pathway | 信号主干通路 | 4.05e-24 | 29.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.05e-24 | 25.5 |
| Hepatitis B | 交叉调控 | 4.76e-24 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.23e-20 | 21.3 |
| Ras signaling pathway | 信号主干通路 | 1.56e-19 | 30.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.21e-18 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 1.56e-18 | 30.5 |
| Prostate cancer | 交叉调控 | 2.34e-18 | 54.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.20e-18 | 53.3 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 3.33e-18 | 17.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 9.92e-18 | 31.9 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.09e-17 | 27.5 |
| Progesterone-mediated oocyte maturation | 其他 | 1.96e-16 | 49.2 |
| Central carbon metabolism in cancer | 交叉调控 | 1.96e-16 | 64.4 |
| Prolactin signaling pathway | 交叉调控 | 1.96e-16 | 64.4 |
| MicroRNAs in cancer | 交叉调控 | 2.97e-16 | 21.2 |
| Pancreatic cancer | 信号主干通路 | 4.49e-16 | 59.3 |
| Focal adhesion | 交叉调控 | 5.26e-16 | 28.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.9_K10</b> — 显著40条, 主干9, 交叉12, 直接OX0, EF=102.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Lipid and atherosclerosis | 交叉调控 | 1.08e-12 | 76.3 |
| Chemical carcinogenesis | 交叉调控 | 1.27e-12 | 68.6 |
| Pathways in cancer | 交叉调控 | 2.79e-12 | 34.7 |
| Estrogen signaling pathway | 交叉调控 | 2.79e-12 | 104.7 |
| Hepatitis B | 交叉调控 | 7.36e-12 | 88.6 |
| Prolactin signaling pathway | 交叉调控 | 7.97e-12 | 175.7 |
| Colorectal cancer | 信号主干通路 | 2.44e-11 | 143.0 |
| Proteoglycans in cancer | 交叉调控 | 2.44e-11 | 70.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.81e-11 | 123.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.37e-09 | 63.7 |
| Pancreatic cancer | 信号主干通路 | 2.44e-09 | 134.9 |
| Focal adhesion | 交叉调控 | 2.52e-09 | 61.2 |
| ErbB signaling pathway | 其他 | 3.66e-09 | 120.6 |
| Human cytomegalovirus infection | 其他 | 4.27e-09 | 54.7 |
| Prostate cancer | 交叉调控 | 6.22e-09 | 105.7 |
| Shigellosis | 其他 | 6.40e-09 | 50.0 |
| Salmonella infection | 其他 | 6.48e-09 | 49.4 |
| HIF-1 signaling pathway | 信号主干通路 | 9.37e-09 | 94.0 |
| Toxoplasmosis | 其他 | 1.02e-08 | 91.5 |
| Relaxin signaling pathway | 交叉调控 | 1.94e-08 | 79.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.9_K20</b> — 显著40条, 主干9, 交叉12, 直接OX1, EF=55.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.23e-19 | 65.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.84e-17 | 57.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.95e-17 | 102.5 |
| Hepatitis B | 交叉调控 | 4.95e-17 | 69.6 |
| Pathways in cancer | 交叉调控 | 8.80e-17 | 27.0 |
| Prolactin signaling pathway | 交叉调控 | 2.26e-16 | 131.8 |
| Chemical carcinogenesis | 交叉调控 | 2.21e-15 | 47.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.12e-14 | 53.1 |
| Human immunodeficiency virus 1 infection | 其他 | 4.87e-14 | 48.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.24e-14 | 66.4 |
| MAPK signaling pathway | 信号主干通路 | 1.07e-12 | 34.9 |
| Toxoplasmosis | 其他 | 1.34e-12 | 73.2 |
| Focal adhesion | 交叉调控 | 1.85e-12 | 45.9 |
| Sphingolipid signaling pathway | 其他 | 1.89e-12 | 68.9 |
| Relaxin signaling pathway | 交叉调控 | 3.22e-12 | 63.6 |
| Fc epsilon RI signaling pathway | 其他 | 3.22e-12 | 105.5 |
| Human cytomegalovirus infection | 其他 | 3.93e-12 | 41.0 |
| Estrogen signaling pathway | 交叉调控 | 4.63e-12 | 59.9 |
| Coronavirus disease | 其他 | 4.63e-12 | 39.8 |
| Shigellosis | 其他 | 7.48e-12 | 37.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.9_K30</b> — 显著40条, 主干10, 交叉12, 直接OX1, EF=38.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 9.08e-23 | 25.7 |
| Proteoglycans in cancer | 交叉调控 | 1.06e-20 | 50.0 |
| Hepatitis B | 交叉调控 | 2.33e-20 | 59.1 |
| Progesterone-mediated oocyte maturation | 其他 | 4.50e-19 | 82.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.09e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 5.49e-17 | 41.3 |
| MAPK signaling pathway | 信号主干通路 | 2.86e-15 | 30.2 |
| Chemical carcinogenesis | 交叉调控 | 1.05e-14 | 34.3 |
| Prolactin signaling pathway | 交叉调控 | 1.26e-14 | 87.9 |
| Epstein-Barr virus infection | 其他 | 7.13e-14 | 37.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.36e-14 | 49.2 |
| Human immunodeficiency virus 1 infection | 其他 | 1.02e-13 | 35.5 |
| Toxoplasmosis | 其他 | 7.08e-13 | 54.9 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.83e-12 | 35.4 |
| Relaxin signaling pathway | 交叉调控 | 2.25e-12 | 47.7 |
| Focal adhesion | 交叉调控 | 2.41e-12 | 34.0 |
| Estrogen signaling pathway | 交叉调控 | 3.45e-12 | 44.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 8.27e-12 | 22.0 |
| GnRH signaling pathway | 其他 | 8.53e-12 | 58.8 |
| Ras signaling pathway | 信号主干通路 | 8.53e-12 | 29.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.08_S0.9_K50</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=45.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.28e-32 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 2.18e-25 | 40.0 |
| Hepatitis B | 交叉调控 | 7.81e-24 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.56e-23 | 23.3 |
| MAPK signaling pathway | 信号主干通路 | 8.06e-21 | 26.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.23e-19 | 22.0 |
| Progesterone-mediated oocyte maturation | 其他 | 5.24e-18 | 53.3 |
| Chemical carcinogenesis | 交叉调控 | 9.71e-18 | 27.4 |
| MicroRNAs in cancer | 交叉调控 | 1.63e-17 | 22.5 |
| Rap1 signaling pathway | 信号主干通路 | 4.31e-17 | 29.3 |
| Prostate cancer | 交叉调控 | 1.59e-16 | 50.7 |
| Ras signaling pathway | 信号主干通路 | 1.59e-16 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.03e-16 | 49.2 |
| Prolactin signaling pathway | 交叉调控 | 2.08e-16 | 64.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.54e-16 | 29.7 |
| Pancreatic cancer | 信号主干通路 | 4.71e-16 | 59.3 |
| Focal adhesion | 交叉调控 | 5.50e-16 | 28.6 |
| Lipid and atherosclerosis | 交叉调控 | 1.34e-15 | 26.7 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.49e-15 | 15.7 |
| Alzheimer disease | 交叉调控 | 3.83e-15 | 17.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.4_K10</b> — 显著40条, 主干11, 交叉14, 直接OX0, EF=100.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.71e-13 | 80.0 |
| Pathways in cancer | 交叉调控 | 7.20e-10 | 30.9 |
| Hepatitis B | 交叉调控 | 2.51e-09 | 75.9 |
| Lipid and atherosclerosis | 交叉调控 | 1.04e-08 | 57.2 |
| Human cytomegalovirus infection | 其他 | 1.10e-08 | 54.7 |
| Chemical carcinogenesis | 交叉调控 | 1.32e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.55e-08 | 26.4 |
| MicroRNAs in cancer | 交叉调控 | 5.38e-08 | 39.7 |
| Estrogen signaling pathway | 交叉调控 | 6.14e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.85e-07 | 29.8 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.06e-07 | 53.1 |
| IL-17 signaling pathway | 交叉调控 | 1.34e-06 | 87.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.57e-06 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.95e-06 | 75.2 |
| TNF signaling pathway | 交叉调控 | 1.95e-06 | 73.2 |
| Toxoplasmosis | 其他 | 1.95e-06 | 73.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.33e-06 | 59.0 |
| negative regulation of production of miRNAs involved in | 其他 | 4.77e-06 | 512.5 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 6.39e-06 | 25.4 |
| negative regulation of apoptotic signaling pathway (GO: | 信号主干通路 | 6.39e-06 | 105.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.4_K20</b> — 显著40条, 主干13, 交叉15, 直接OX0, EF=47.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.17e-19 | 65.0 |
| Pathways in cancer | 交叉调控 | 2.23e-18 | 29.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.48e-12 | 21.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.88e-12 | 25.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.02e-12 | 82.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.96e-12 | 47.8 |
| Lipid and atherosclerosis | 交叉调控 | 8.45e-12 | 42.9 |
| Estrogen signaling pathway | 交叉调控 | 1.28e-11 | 59.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.28e-11 | 59.0 |
| Chemical carcinogenesis | 交叉调控 | 1.38e-11 | 38.6 |
| Hepatitis B | 交叉调控 | 3.30e-11 | 50.6 |
| Bladder cancer | 交叉调控 | 3.30e-11 | 150.0 |
| Prostate cancer | 交叉调控 | 5.83e-11 | 74.0 |
| MicroRNAs in cancer | 交叉调控 | 9.53e-11 | 29.8 |
| Th17 cell differentiation | 交叉调控 | 9.95e-11 | 67.1 |
| HIF-1 signaling pathway | 信号主干通路 | 1.06e-10 | 65.8 |
| Human cytomegalovirus infection | 其他 | 2.91e-10 | 36.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.51e-10 | 25.4 |
| Pancreatic cancer | 信号主干通路 | 9.83e-10 | 80.9 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 2.44e-09 | 69.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.4_K30</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=43.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.58e-24 | 56.7 |
| Pathways in cancer | 交叉调控 | 4.77e-23 | 25.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.45e-17 | 24.1 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.45e-17 | 19.8 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.64e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 8.64e-17 | 41.3 |
| Chemical carcinogenesis | 交叉调控 | 2.79e-16 | 37.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.09e-15 | 42.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.63e-15 | 25.4 |
| Hepatitis B | 交叉调控 | 9.11e-15 | 46.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.85e-14 | 39.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.93e-14 | 25.1 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.07e-13 | 49.2 |
| Human cytomegalovirus infection | 其他 | 2.47e-13 | 33.4 |
| HIF-1 signaling pathway | 信号主干通路 | 6.84e-13 | 56.4 |
| Pancreatic cancer | 信号主干通路 | 2.59e-12 | 71.9 |
| positive regulation of intracellular signal transductio | 其他 | 2.99e-12 | 17.5 |
| Estrogen signaling pathway | 交叉调控 | 4.48e-12 | 44.9 |
| Rap1 signaling pathway | 信号主干通路 | 4.48e-12 | 32.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.33e-12 | 42.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.4_K50</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=34.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.90e-29 | 21.6 |
| Proteoglycans in cancer | 交叉调控 | 2.35e-25 | 40.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.06e-20 | 23.2 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.03e-19 | 21.3 |
| MAPK signaling pathway | 信号主干通路 | 1.43e-17 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.96e-17 | 31.7 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 7.12e-17 | 16.5 |
| Hepatitis B | 交叉调控 | 8.83e-17 | 35.4 |
| Prostate cancer | 交叉调控 | 2.99e-16 | 50.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.64e-16 | 49.2 |
| Central carbon metabolism in cancer | 交叉调控 | 3.64e-16 | 64.4 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 4.27e-16 | 80.5 |
| MicroRNAs in cancer | 交叉调控 | 5.25e-16 | 21.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 5.35e-16 | 29.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.62e-16 | 16.2 |
| Rap1 signaling pathway | 信号主干通路 | 1.60e-15 | 27.3 |
| Lipid and atherosclerosis | 交叉调控 | 2.05e-15 | 26.7 |
| peptidyl-tyrosine phosphorylation (GO:0018108) | 信号主干通路 | 3.40e-15 | 65.4 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 3.90e-15 | 21.7 |
| positive regulation of intracellular signal transductio | 其他 | 4.19e-15 | 14.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.5_K10</b> — 显著40条, 主干12, 交叉13, 直接OX0, EF=90.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.57e-13 | 80.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.23e-10 | 66.7 |
| Pathways in cancer | 交叉调控 | 4.70e-10 | 30.9 |
| Estrogen signaling pathway | 交叉调控 | 6.65e-10 | 89.8 |
| Hepatitis B | 交叉调控 | 1.47e-09 | 75.9 |
| Chemical carcinogenesis | 交叉调控 | 1.29e-08 | 51.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.64e-08 | 26.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.39e-08 | 73.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.91e-07 | 29.8 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 6.42e-07 | 210.3 |
| Human cytomegalovirus infection | 其他 | 7.26e-07 | 45.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.45e-06 | 63.7 |
| IL-17 signaling pathway | 交叉调控 | 1.46e-06 | 87.2 |
| Prostate cancer | 交叉调控 | 1.49e-06 | 84.5 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 1.54e-06 | 19.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.54e-06 | 82.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.60e-06 | 30.4 |
| HIF-1 signaling pathway | 信号主干通路 | 1.99e-06 | 75.2 |
| Toxoplasmosis | 其他 | 2.05e-06 | 73.2 |
| MicroRNAs in cancer | 交叉调控 | 2.05e-06 | 33.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.5_K20</b> — 显著40条, 主干12, 交叉14, 直接OX0, EF=52.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.72e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 2.24e-18 | 29.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.55e-16 | 58.4 |
| Chemical carcinogenesis | 交叉调控 | 3.72e-15 | 47.2 |
| Lipid and atherosclerosis | 交叉调控 | 9.72e-14 | 47.7 |
| Hepatitis B | 交叉调控 | 5.38e-13 | 56.9 |
| Prostate cancer | 交叉调控 | 6.76e-13 | 84.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.61e-13 | 82.0 |
| HIF-1 signaling pathway | 信号主干通路 | 1.37e-12 | 75.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.56e-12 | 21.5 |
| MicroRNAs in cancer | 交叉调控 | 1.93e-12 | 33.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.98e-12 | 25.5 |
| Human cytomegalovirus infection | 其他 | 5.83e-12 | 41.0 |
| Estrogen signaling pathway | 交叉调控 | 6.68e-12 | 59.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.94e-12 | 59.0 |
| Pancreatic cancer | 信号主干通路 | 7.93e-12 | 94.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.06e-11 | 27.9 |
| PD-L1 expression and PD-1 checkpoint pathway in cancer | 信号主干通路 | 2.32e-11 | 80.6 |
| TNF signaling pathway | 交叉调控 | 1.13e-10 | 64.1 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 1.70e-10 | 27.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.5_K30</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=45.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.57e-24 | 56.7 |
| Pathways in cancer | 交叉调控 | 4.74e-23 | 25.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.35e-17 | 24.1 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.35e-17 | 19.8 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.60e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 8.60e-17 | 41.3 |
| Chemical carcinogenesis | 交叉调控 | 2.77e-16 | 37.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.09e-15 | 42.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.57e-15 | 25.4 |
| Hepatitis B | 交叉调控 | 9.06e-15 | 46.4 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.91e-14 | 25.1 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.07e-13 | 49.2 |
| Human cytomegalovirus infection | 其他 | 2.46e-13 | 33.4 |
| HIF-1 signaling pathway | 信号主干通路 | 6.80e-13 | 56.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.75e-12 | 36.3 |
| Pancreatic cancer | 信号主干通路 | 2.58e-12 | 71.9 |
| positive regulation of intracellular signal transductio | 其他 | 2.87e-12 | 17.5 |
| Estrogen signaling pathway | 交叉调控 | 4.71e-12 | 44.9 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.11e-12 | 42.4 |
| MicroRNAs in cancer | 交叉调控 | 5.93e-12 | 24.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.5_K50</b> — 显著40条, 主干15, 交叉14, 直接OX0, EF=34.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 2.27e-34 | 23.9 |
| Proteoglycans in cancer | 交叉调控 | 3.12e-27 | 42.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.55e-22 | 24.3 |
| MAPK signaling pathway | 信号主干通路 | 8.66e-21 | 26.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.95e-20 | 21.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 5.52e-19 | 33.7 |
| Hepatitis B | 交叉调控 | 1.50e-18 | 38.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.50e-18 | 31.2 |
| Prostate cancer | 交叉调控 | 3.18e-18 | 54.9 |
| Central carbon metabolism in cancer | 交叉调控 | 3.49e-18 | 70.3 |
| Ras signaling pathway | 信号主干通路 | 5.02e-18 | 28.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 6.80e-18 | 24.4 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 5.15e-17 | 16.5 |
| positive regulation of intracellular signal transductio | 其他 | 2.45e-16 | 15.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.45e-16 | 25.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.61e-16 | 49.2 |
| Chemical carcinogenesis | 交叉调控 | 2.68e-16 | 25.7 |
| MicroRNAs in cancer | 交叉调控 | 3.92e-16 | 21.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.09e-16 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.62e-16 | 35.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.7_K10</b> — 显著40条, 主干8, 交叉10, 直接OX0, EF=96.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.43e-13 | 80.0 |
| Chemical carcinogenesis | 交叉调控 | 2.54e-10 | 60.0 |
| Pathways in cancer | 交叉调控 | 4.60e-10 | 30.9 |
| Lipid and atherosclerosis | 交叉调控 | 9.99e-09 | 57.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.99e-08 | 26.4 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.99e-08 | 76.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.08e-08 | 59.4 |
| Estrogen signaling pathway | 交叉调控 | 9.41e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.00e-07 | 29.8 |
| Hepatitis B | 交叉调控 | 1.82e-07 | 63.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.90e-07 | 26.3 |
| positive regulation of pri-miRNA transcription by RNA p | 其他 | 1.90e-07 | 241.2 |
| regulation of nitric-oxide synthase activity (GO:005099 | 其他 | 2.89e-07 | 210.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.77e-07 | 53.1 |
| regulation of pri-miRNA transcription by RNA polymerase | 其他 | 4.57e-07 | 182.2 |
| Human cytomegalovirus infection | 其他 | 7.10e-07 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.07e-06 | 30.4 |
| positive regulation of cellular biosynthetic process (G | 其他 | 1.07e-06 | 56.9 |
| negative regulation of production of miRNAs involved in | 其他 | 1.41e-06 | 512.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.84e-06 | 82.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.7_K20</b> — 显著40条, 主干10, 交叉16, 直接OX1, EF=61.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.72e-22 | 70.0 |
| Pathways in cancer | 交叉调控 | 1.85e-20 | 30.9 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.76e-17 | 102.5 |
| Hepatitis B | 交叉调控 | 4.76e-17 | 69.6 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.73e-16 | 58.4 |
| Lipid and atherosclerosis | 交叉调控 | 7.62e-16 | 52.4 |
| Prolactin signaling pathway | 交叉调控 | 4.48e-14 | 117.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 5.61e-14 | 27.6 |
| Estrogen signaling pathway | 交叉调控 | 8.66e-14 | 67.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.80e-14 | 66.4 |
| Chemical carcinogenesis | 交叉调控 | 1.41e-13 | 42.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.29e-13 | 30.4 |
| IL-17 signaling pathway | 交叉调控 | 3.32e-13 | 87.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 4.99e-13 | 21.5 |
| C-type lectin receptor signaling pathway | 交叉调控 | 7.01e-13 | 78.8 |
| TNF signaling pathway | 交叉调控 | 1.10e-12 | 73.2 |
| Toxoplasmosis | 其他 | 1.10e-12 | 73.2 |
| Relaxin signaling pathway | 交叉调控 | 3.27e-12 | 63.6 |
| Human cytomegalovirus infection | 其他 | 4.01e-12 | 41.0 |
| Pancreatic cancer | 信号主干通路 | 6.53e-12 | 94.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.7_K30</b> — 显著40条, 主干14, 交叉15, 直接OX1, EF=39.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 8.93e-23 | 25.7 |
| Proteoglycans in cancer | 交叉调控 | 9.39e-23 | 53.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.90e-17 | 24.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.41e-17 | 27.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.10e-17 | 75.2 |
| Lipid and atherosclerosis | 交叉调控 | 8.10e-17 | 41.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.03e-16 | 18.7 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.23e-15 | 42.5 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.78e-15 | 54.1 |
| Hepatitis B | 交叉调控 | 8.54e-15 | 46.4 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 5.58e-14 | 29.3 |
| positive regulation of intracellular signal transductio | 其他 | 7.12e-14 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 7.12e-14 | 46.7 |
| MAPK signaling pathway | 信号主干通路 | 1.24e-13 | 27.9 |
| TNF signaling pathway | 交叉调控 | 9.06e-13 | 54.9 |
| Toxoplasmosis | 其他 | 9.06e-13 | 54.9 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.01e-12 | 36.3 |
| Prolactin signaling pathway | 交叉调控 | 1.33e-12 | 78.1 |
| Pancreatic cancer | 信号主干通路 | 2.43e-12 | 71.9 |
| Focal adhesion | 交叉调控 | 2.84e-12 | 34.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.7_K50</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=34.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 4.67e-38 | 25.5 |
| Proteoglycans in cancer | 交叉调控 | 2.26e-25 | 40.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.51e-24 | 25.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 6.01e-23 | 23.3 |
| MAPK signaling pathway | 信号主干通路 | 1.81e-22 | 27.9 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 9.61e-21 | 35.7 |
| Rap1 signaling pathway | 信号主干通路 | 3.40e-20 | 33.2 |
| Prostate cancer | 交叉调控 | 4.00e-20 | 59.2 |
| Ras signaling pathway | 信号主干通路 | 1.36e-19 | 30.0 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.64e-19 | 25.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 2.29e-19 | 40.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 7.84e-19 | 17.9 |
| Lipid and atherosclerosis | 交叉调控 | 1.59e-18 | 30.5 |
| Central carbon metabolism in cancer | 交叉调控 | 2.99e-18 | 70.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.82e-18 | 27.5 |
| MicroRNAs in cancer | 交叉调控 | 1.36e-17 | 22.5 |
| Focal adhesion | 交叉调控 | 1.88e-17 | 30.6 |
| Hepatitis B | 交叉调控 | 3.53e-17 | 35.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.04e-16 | 13.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.68e-16 | 15.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.9_K10</b> — 显著40条, 主干9, 交叉12, 直接OX0, EF=102.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Lipid and atherosclerosis | 交叉调控 | 1.08e-12 | 76.3 |
| Chemical carcinogenesis | 交叉调控 | 1.27e-12 | 68.6 |
| Pathways in cancer | 交叉调控 | 2.79e-12 | 34.7 |
| Estrogen signaling pathway | 交叉调控 | 2.79e-12 | 104.7 |
| Hepatitis B | 交叉调控 | 7.36e-12 | 88.6 |
| Prolactin signaling pathway | 交叉调控 | 7.97e-12 | 175.7 |
| Colorectal cancer | 信号主干通路 | 2.44e-11 | 143.0 |
| Proteoglycans in cancer | 交叉调控 | 2.44e-11 | 70.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.81e-11 | 123.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.37e-09 | 63.7 |
| Pancreatic cancer | 信号主干通路 | 2.44e-09 | 134.9 |
| Focal adhesion | 交叉调控 | 2.52e-09 | 61.2 |
| ErbB signaling pathway | 其他 | 3.66e-09 | 120.6 |
| Human cytomegalovirus infection | 其他 | 4.27e-09 | 54.7 |
| Prostate cancer | 交叉调控 | 6.22e-09 | 105.7 |
| Shigellosis | 其他 | 6.40e-09 | 50.0 |
| Salmonella infection | 其他 | 6.48e-09 | 49.4 |
| HIF-1 signaling pathway | 信号主干通路 | 9.37e-09 | 94.0 |
| Toxoplasmosis | 其他 | 1.02e-08 | 91.5 |
| Relaxin signaling pathway | 交叉调控 | 1.94e-08 | 79.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.9_K20</b> — 显著40条, 主干13, 交叉11, 直接OX1, EF=55.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 4.40e-18 | 29.0 |
| Proteoglycans in cancer | 交叉调控 | 9.68e-18 | 60.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.16e-17 | 57.2 |
| Hepatitis B | 交叉调控 | 6.72e-15 | 63.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.37e-15 | 92.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.66e-14 | 53.1 |
| Prolactin signaling pathway | 交叉调控 | 4.40e-14 | 117.1 |
| Human immunodeficiency virus 1 infection | 其他 | 5.18e-14 | 48.3 |
| Chemical carcinogenesis | 交叉调控 | 1.54e-13 | 42.9 |
| Toxoplasmosis | 其他 | 1.52e-12 | 73.2 |
| Focal adhesion | 交叉调控 | 2.06e-12 | 45.9 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 3.82e-12 | 25.4 |
| Relaxin signaling pathway | 交叉调控 | 4.02e-12 | 63.6 |
| Ras signaling pathway | 信号主干通路 | 6.33e-12 | 39.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.33e-12 | 59.0 |
| Pancreatic cancer | 信号主干通路 | 7.28e-12 | 94.4 |
| Shigellosis | 其他 | 8.83e-12 | 37.5 |
| ErbB signaling pathway | 其他 | 1.44e-11 | 84.4 |
| Colorectal cancer | 信号主干通路 | 1.48e-11 | 83.4 |
| GnRH signaling pathway | 其他 | 2.47e-11 | 77.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.9_K30</b> — 显著40条, 主干13, 交叉12, 直接OX3, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.14e-24 | 27.0 |
| Proteoglycans in cancer | 交叉调控 | 1.04e-18 | 46.7 |
| Lipid and atherosclerosis | 交叉调控 | 1.37e-18 | 44.5 |
| Hepatitis B | 交叉调控 | 1.92e-18 | 54.8 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 5.12e-17 | 75.2 |
| Prolactin signaling pathway | 交叉调控 | 1.22e-16 | 97.6 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 1.61e-15 | 22.5 |
| MAPK signaling pathway | 信号主干通路 | 2.87e-15 | 30.2 |
| Toxoplasmosis | 其他 | 1.26e-14 | 61.0 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 4.80e-14 | 38.9 |
| Focal adhesion | 交叉调控 | 6.52e-14 | 37.4 |
| Epstein-Barr virus infection | 其他 | 6.52e-14 | 37.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 7.70e-14 | 49.2 |
| Human immunodeficiency virus 1 infection | 其他 | 9.42e-14 | 35.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.11e-13 | 23.7 |
| MicroRNAs in cancer | 交叉调控 | 1.37e-13 | 26.5 |
| Choline metabolism in cancer | 信号主干通路 | 1.79e-13 | 62.8 |
| Ras signaling pathway | 信号主干通路 | 2.08e-13 | 32.4 |
| Chemical carcinogenesis | 交叉调控 | 2.71e-13 | 31.5 |
| Thyroid hormone signaling pathway | 其他 | 1.05e-12 | 50.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.0_S0.9_K50</b> — 显著40条, 主干15, 交叉12, 直接OX1, EF=27.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.30e-32 | 23.2 |
| Hepatitis B | 交叉调控 | 8.88e-22 | 43.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.38e-21 | 22.3 |
| Proteoglycans in cancer | 交叉调控 | 3.65e-20 | 34.0 |
| MicroRNAs in cancer | 交叉调控 | 9.02e-19 | 23.8 |
| Epstein-Barr virus infection | 其他 | 9.02e-19 | 32.5 |
| MAPK signaling pathway | 信号主干通路 | 9.00e-18 | 23.7 |
| Focal adhesion | 交叉调控 | 2.89e-17 | 30.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.56e-16 | 19.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.47e-16 | 49.2 |
| Progesterone-mediated oocyte maturation | 其他 | 2.47e-16 | 49.2 |
| Rap1 signaling pathway | 信号主干通路 | 1.51e-15 | 27.3 |
| Lipid and atherosclerosis | 交叉调控 | 1.93e-15 | 26.7 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.53e-15 | 15.7 |
| Prostate cancer | 交叉调控 | 9.24e-15 | 46.5 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.50e-14 | 27.6 |
| Prolactin signaling pathway | 交叉调控 | 1.67e-14 | 58.6 |
| Human papillomavirus infection | 其他 | 2.24e-14 | 18.6 |
| Viral carcinogenesis | 信号主干通路 | 2.39e-14 | 26.3 |
| Ras signaling pathway | 信号主干通路 | 1.28e-13 | 23.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.4_K10</b> — 显著40条, 主干6, 交叉15, 直接OX0, EF=89.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 5.30e-06 | 361.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 5.43e-05 | 19.8 |
| bicarbonate transport (GO:0015701) | 其他 | 8.94e-05 | 175.7 |
| response to peptide (GO:1901652) | 其他 | 8.94e-05 | 157.7 |
| Pathways in cancer | 交叉调控 | 1.27e-04 | 19.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.49e-04 | 21.3 |
| MicroRNAs in cancer | 交叉调控 | 3.09e-04 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.09e-04 | 61.5 |
| Estrogen signaling pathway | 交叉调控 | 5.52e-04 | 44.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.52e-04 | 44.2 |
| Hepatitis B | 交叉调控 | 7.46e-04 | 38.0 |
| cellular response to UV-A (GO:0071492) | 其他 | 7.55e-04 | 372.7 |
| organic anion transport (GO:0015711) | 其他 | 7.55e-04 | 62.1 |
| response to UV-A (GO:0070141) | 其他 | 9.59e-04 | 292.9 |
| Proteoglycans in cancer | 交叉调控 | 1.31e-03 | 30.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.34e-03 | 28.6 |
| Bladder cancer | 交叉调控 | 1.50e-03 | 100.0 |
| Chemical carcinogenesis | 交叉调控 | 1.50e-03 | 25.7 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 1.58e-03 | 215.8 |
| positive regulation of vascular associated smooth muscl | 其他 | 1.72e-03 | 195.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.4_K20</b> — 显著25条, 主干2, 交叉10, 直接OX0, EF=91.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 5.61e-16 | 422.1 |
| bicarbonate transport (GO:0015701) | 其他 | 8.21e-13 | 205.0 |
| organic anion transport (GO:0015711) | 其他 | 8.76e-10 | 72.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.27e-04 | 41.0 |
| MicroRNAs in cancer | 交叉调控 | 3.51e-04 | 16.5 |
| response to peptide (GO:1901652) | 其他 | 1.03e-03 | 78.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.15e-03 | 9.9 |
| regulation of anion transport (GO:0044070) | 其他 | 2.15e-03 | 256.2 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 2.43e-03 | 47.3 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 2.47e-03 | 205.0 |
| cellular response to UV-A (GO:0071492) | 其他 | 2.64e-03 | 186.4 |
| Pathways in cancer | 交叉调控 | 3.42e-03 | 9.7 |
| Serotonergic synapse | 交叉调控 | 3.55e-03 | 27.2 |
| response to UV-A (GO:0070141) | 其他 | 3.71e-03 | 146.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 3.71e-03 | 10.6 |
| Proximal tubule bicarbonate reclamation | 其他 | 3.76e-03 | 89.1 |
| Estrogen signaling pathway | 交叉调控 | 4.08e-03 | 22.4 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.08e-03 | 22.1 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 5.31e-03 | 113.9 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 5.44e-03 | 107.9 |
| ... 还有 5 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.4_K30</b> — 显著30条, 主干3, 交叉5, 直接OX0, EF=96.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.64e-20 | 361.8 |
| bicarbonate transport (GO:0015701) | 其他 | 3.84e-16 | 175.7 |
| organic anion transport (GO:0015711) | 其他 | 4.43e-12 | 62.1 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 1.57e-07 | 273.3 |
| Serotonergic synapse | 交叉调控 | 8.17e-07 | 36.3 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 4.17e-06 | 410.0 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 4.55e-06 | 52.6 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 9.71e-06 | 292.9 |
| Gap junction | 交叉调控 | 2.89e-04 | 31.1 |
| Pathways in cancer | 交叉调控 | 2.89e-04 | 9.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.15e-04 | 27.3 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 5.37e-04 | 82.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 5.97e-04 | 75.9 |
| MicroRNAs in cancer | 交叉调控 | 1.62e-03 | 11.0 |
| response to peptide (GO:1901652) | 其他 | 1.64e-03 | 52.6 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 1.84e-03 | 48.8 |
| Neuroactive ligand-receptor interaction | 其他 | 2.17e-03 | 10.0 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 2.44e-03 | 20.1 |
| regulation of anion transport (GO:0044070) | 其他 | 2.77e-03 | 170.8 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 3.31e-03 | 36.0 |
| ... 还有 10 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.4_K50</b> — 显著31条, 主干3, 交叉6, 直接OX0, EF=75.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 9.80e-22 | 294.1 |
| bicarbonate transport (GO:0015701) | 其他 | 4.77e-17 | 142.9 |
| organic anion transport (GO:0015711) | 其他 | 2.33e-14 | 55.6 |
| Serotonergic synapse | 交叉调控 | 1.87e-07 | 31.0 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 6.47e-07 | 200.0 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 6.47e-07 | 46.2 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 1.08e-05 | 300.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 3.13e-05 | 214.3 |
| Neuroactive ligand-receptor interaction | 其他 | 2.29e-04 | 10.3 |
| Pathways in cancer | 交叉调控 | 3.31e-04 | 7.5 |
| Proximal tubule bicarbonate reclamation | 其他 | 3.52e-04 | 65.2 |
| Gap junction | 交叉调控 | 6.67e-04 | 22.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 9.45e-04 | 20.0 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 1.72e-03 | 60.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 1.91e-03 | 55.6 |
| Calcium signaling pathway | 信号主干通路 | 2.01e-03 | 10.4 |
| chemical synaptic transmission (GO:0007268) | 其他 | 2.63e-03 | 9.8 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 2.76e-03 | 200.0 |
| cellular metal ion homeostasis (GO:0006875) | 其他 | 3.70e-03 | 19.2 |
| response to peptide (GO:1901652) | 其他 | 3.91e-03 | 38.5 |
| ... 还有 11 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.5_K10</b> — 显著40条, 主干6, 交叉15, 直接OX0, EF=98.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 6.43e-09 | 482.4 |
| bicarbonate transport (GO:0015701) | 其他 | 5.08e-07 | 234.3 |
| organic anion transport (GO:0015711) | 其他 | 1.71e-05 | 82.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.71e-05 | 19.8 |
| response to peptide (GO:1901652) | 其他 | 6.31e-05 | 157.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.13e-04 | 21.3 |
| Pathways in cancer | 交叉调控 | 1.26e-04 | 19.3 |
| MicroRNAs in cancer | 交叉调控 | 3.05e-04 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.05e-04 | 61.5 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-04 | 44.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-04 | 44.2 |
| cellular response to UV-A (GO:0071492) | 其他 | 6.38e-04 | 372.7 |
| Proximal tubule bicarbonate reclamation | 其他 | 6.45e-04 | 178.3 |
| Hepatitis B | 交叉调控 | 6.45e-04 | 38.0 |
| response to UV-A (GO:0070141) | 其他 | 9.04e-04 | 292.9 |
| Proteoglycans in cancer | 交叉调控 | 1.15e-03 | 30.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.20e-03 | 28.6 |
| Bladder cancer | 交叉调控 | 1.36e-03 | 100.0 |
| Chemical carcinogenesis | 交叉调控 | 1.36e-03 | 25.7 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 1.48e-03 | 215.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.5_K20</b> — 显著32条, 主干3, 交叉10, 直接OX0, EF=102.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 4.99e-13 | 361.8 |
| bicarbonate transport (GO:0015701) | 其他 | 3.10e-10 | 175.7 |
| organic anion transport (GO:0015711) | 其他 | 1.03e-07 | 62.1 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 1.52e-05 | 307.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.25e-04 | 41.0 |
| Serotonergic synapse | 交叉调控 | 1.35e-04 | 36.3 |
| MicroRNAs in cancer | 交叉调控 | 2.58e-04 | 16.5 |
| response to peptide (GO:1901652) | 其他 | 8.46e-04 | 78.8 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 8.46e-04 | 410.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 1.33e-03 | 292.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.33e-03 | 9.9 |
| regulation of anion transport (GO:0044070) | 其他 | 1.33e-03 | 256.2 |
| negative regulation of autophagy (GO:0010507) | 信号主干通路 | 1.33e-03 | 52.1 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 1.60e-03 | 47.3 |
| negative regulation of cellular catabolic process (GO:0 | 其他 | 1.74e-03 | 44.6 |
| cellular response to UV-A (GO:0071492) | 其他 | 1.93e-03 | 186.4 |
| Pathways in cancer | 交叉调控 | 2.68e-03 | 9.7 |
| response to UV-A (GO:0070141) | 其他 | 2.91e-03 | 146.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.91e-03 | 10.6 |
| Proximal tubule bicarbonate reclamation | 其他 | 3.68e-03 | 89.1 |
| ... 还有 12 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.5_K30</b> — 显著33条, 主干4, 交叉7, 直接OX0, EF=85.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.62e-20 | 361.8 |
| bicarbonate transport (GO:0015701) | 其他 | 3.96e-16 | 175.7 |
| organic anion transport (GO:0015711) | 其他 | 4.58e-12 | 62.1 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 1.62e-07 | 273.3 |
| Serotonergic synapse | 交叉调控 | 8.09e-07 | 36.3 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 4.31e-06 | 410.0 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 4.70e-06 | 52.6 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 1.00e-05 | 292.9 |
| Pathways in cancer | 交叉调控 | 2.92e-05 | 10.3 |
| Gap junction | 交叉调控 | 2.35e-04 | 31.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.12e-04 | 27.3 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 5.55e-04 | 82.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 6.16e-04 | 75.9 |
| MicroRNAs in cancer | 交叉调控 | 1.60e-03 | 11.0 |
| response to peptide (GO:1901652) | 其他 | 1.69e-03 | 52.6 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 1.91e-03 | 48.8 |
| Neuroactive ligand-receptor interaction | 其他 | 2.15e-03 | 10.0 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 2.52e-03 | 20.1 |
| regulation of anion transport (GO:0044070) | 其他 | 2.86e-03 | 170.8 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 3.42e-03 | 36.0 |
| ... 还有 13 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.5_K50</b> — 显著31条, 主干2, 交叉6, 直接OX0, EF=78.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 7.36e-22 | 301.5 |
| bicarbonate transport (GO:0015701) | 其他 | 3.59e-17 | 146.4 |
| organic anion transport (GO:0015711) | 其他 | 1.70e-14 | 56.9 |
| Serotonergic synapse | 交叉调控 | 1.55e-07 | 31.7 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 5.50e-07 | 205.0 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 5.50e-07 | 47.3 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 9.93e-06 | 307.5 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 2.89e-05 | 219.6 |
| Neuroactive ligand-receptor interaction | 其他 | 1.92e-04 | 10.5 |
| Pathways in cancer | 交叉调控 | 2.70e-04 | 7.7 |
| Proximal tubule bicarbonate reclamation | 其他 | 3.24e-04 | 66.8 |
| Gap junction | 交叉调控 | 6.00e-04 | 23.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 8.50e-04 | 20.5 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 1.59e-03 | 61.5 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 1.76e-03 | 56.9 |
| Calcium signaling pathway | 信号主干通路 | 1.77e-03 | 10.7 |
| chemical synaptic transmission (GO:0007268) | 其他 | 2.26e-03 | 10.0 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 2.61e-03 | 205.0 |
| cellular metal ion homeostasis (GO:0006875) | 其他 | 3.34e-03 | 19.7 |
| response to peptide (GO:1901652) | 其他 | 3.61e-03 | 39.4 |
| ... 还有 11 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.7_K10</b> — 显著40条, 主干6, 交叉16, 直接OX0, EF=132.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 3.75e-06 | 615.0 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 3.75e-06 | 126.2 |
| Serotonergic synapse | 交叉调控 | 1.70e-05 | 72.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.97e-05 | 19.8 |
| response to peptide (GO:1901652) | 其他 | 7.27e-05 | 157.7 |
| Pathways in cancer | 交叉调控 | 1.26e-04 | 19.3 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.30e-04 | 21.3 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 1.34e-04 | 820.0 |
| negative regulation of autophagy (GO:0010507) | 信号主干通路 | 1.47e-04 | 104.2 |
| negative regulation of cellular catabolic process (GO:0 | 其他 | 1.87e-04 | 89.1 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 1.87e-04 | 585.7 |
| MicroRNAs in cancer | 交叉调控 | 2.77e-04 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.77e-04 | 61.5 |
| Neuroactive ligand-receptor interaction | 其他 | 2.77e-04 | 24.0 |
| cellular response to UV-A (GO:0071492) | 其他 | 4.41e-04 | 372.7 |
| Estrogen signaling pathway | 交叉调控 | 4.67e-04 | 44.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 4.67e-04 | 44.2 |
| Hepatitis B | 交叉调控 | 6.45e-04 | 38.0 |
| response to UV-A (GO:0070141) | 其他 | 6.62e-04 | 292.9 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 1.05e-03 | 45.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.7_K20</b> — 显著39条, 主干4, 交叉12, 直接OX0, EF=105.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 5.61e-16 | 422.1 |
| bicarbonate transport (GO:0015701) | 其他 | 8.19e-13 | 205.0 |
| organic anion transport (GO:0015711) | 其他 | 8.74e-10 | 72.5 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 2.05e-08 | 410.0 |
| Serotonergic synapse | 交叉调控 | 4.91e-08 | 54.4 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 4.68e-07 | 78.8 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 6.92e-07 | 615.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 2.01e-06 | 439.3 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 1.12e-04 | 123.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 1.25e-04 | 113.9 |
| Pathways in cancer | 交叉调控 | 2.63e-04 | 11.6 |
| MicroRNAs in cancer | 交叉调控 | 2.63e-04 | 16.5 |
| Neuroactive ligand-receptor interaction | 其他 | 3.34e-04 | 15.0 |
| response to peptide (GO:1901652) | 其他 | 3.44e-04 | 78.8 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 3.53e-04 | 73.2 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 3.53e-04 | 30.1 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.76e-04 | 9.9 |
| regulation of anion transport (GO:0044070) | 其他 | 7.76e-04 | 256.2 |
| negative regulation of autophagy (GO:0010507) | 信号主干通路 | 7.76e-04 | 52.1 |
| negative regulation of cellular catabolic process (GO:0 | 其他 | 1.16e-03 | 44.6 |
| ... 还有 19 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.7_K30</b> — 显著32条, 主干3, 交叉7, 直接OX0, EF=94.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.50e-20 | 361.8 |
| bicarbonate transport (GO:0015701) | 其他 | 3.75e-16 | 175.7 |
| organic anion transport (GO:0015711) | 其他 | 4.33e-12 | 62.1 |
| Serotonergic synapse | 交叉调控 | 1.42e-08 | 42.3 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 9.28e-08 | 63.1 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 1.15e-07 | 273.3 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 3.26e-06 | 410.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 9.48e-06 | 292.9 |
| Gap junction | 交叉调控 | 2.52e-04 | 31.1 |
| Neuroactive ligand-receptor interaction | 其他 | 2.52e-04 | 12.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.98e-04 | 27.3 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 5.24e-04 | 82.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 5.82e-04 | 75.9 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 1.29e-03 | 273.3 |
| response to peptide (GO:1901652) | 其他 | 1.40e-03 | 52.6 |
| anterograde trans-synaptic signaling (GO:0098916) | 其他 | 1.40e-03 | 14.0 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 1.50e-03 | 48.8 |
| MicroRNAs in cancer | 交叉调控 | 1.53e-03 | 11.0 |
| Pathways in cancer | 交叉调控 | 1.71e-03 | 7.7 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 2.01e-03 | 20.1 |
| ... 还有 12 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.7_K50</b> — 显著31条, 主干3, 交叉5, 直接OX0, EF=84.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.65e-19 | 293.3 |
| bicarbonate transport (GO:0015701) | 其他 | 3.91e-15 | 142.5 |
| organic anion transport (GO:0015711) | 其他 | 5.59e-13 | 56.0 |
| Serotonergic synapse | 交叉调控 | 8.52e-08 | 34.3 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 3.34e-07 | 51.1 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 3.34e-07 | 221.6 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 7.53e-06 | 332.4 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 2.19e-05 | 237.5 |
| Neuroactive ligand-receptor interaction | 其他 | 1.09e-04 | 11.4 |
| Pathways in cancer | 交叉调控 | 1.43e-04 | 8.3 |
| Proximal tubule bicarbonate reclamation | 其他 | 2.52e-04 | 72.3 |
| Gap junction | 交叉调控 | 4.31e-04 | 25.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 6.11e-04 | 22.2 |
| Calcium signaling pathway | 信号主干通路 | 1.19e-03 | 11.5 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 1.21e-03 | 66.5 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 1.34e-03 | 61.6 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 2.39e-03 | 221.6 |
| cellular metal ion homeostasis (GO:0006875) | 其他 | 2.59e-03 | 21.3 |
| response to peptide (GO:1901652) | 其他 | 3.00e-03 | 42.6 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 3.44e-03 | 39.6 |
| ... 还有 11 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.9_K10</b> — 显著40条, 主干7, 交叉15, 直接OX0, EF=143.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 2.36e-09 | 820.0 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 1.38e-08 | 157.7 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 1.07e-07 | 1230.0 |
| Serotonergic synapse | 交叉调控 | 1.10e-07 | 90.7 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 2.82e-07 | 878.6 |
| Neuroactive ligand-receptor interaction | 其他 | 1.40e-05 | 30.1 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 1.47e-05 | 246.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 1.56e-05 | 227.8 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 2.13e-05 | 60.3 |
| response to peptide (GO:1901652) | 其他 | 3.65e-05 | 157.7 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 4.07e-05 | 146.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.50e-05 | 21.3 |
| Pathways in cancer | 交叉调控 | 8.28e-05 | 19.3 |
| negative regulation of autophagy (GO:0010507) | 信号主干通路 | 9.38e-05 | 104.2 |
| negative regulation of cellular catabolic process (GO:0 | 其他 | 1.28e-04 | 89.1 |
| anterograde trans-synaptic signaling (GO:0098916) | 其他 | 1.28e-04 | 33.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.61e-04 | 16.5 |
| Gap junction | 交叉调控 | 1.72e-04 | 69.9 |
| MicroRNAs in cancer | 交叉调控 | 1.72e-04 | 26.5 |
| Inflammatory mediator regulation of TRP channels | 信号主干通路 | 1.72e-04 | 62.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.9_K20</b> — 显著40条, 主干4, 交叉12, 直接OX0, EF=103.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 8.09e-08 | 410.0 |
| Serotonergic synapse | 交叉调控 | 1.16e-07 | 54.4 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 1.23e-06 | 78.8 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 1.51e-06 | 615.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 3.97e-06 | 439.3 |
| Pathways in cancer | 交叉调控 | 2.87e-05 | 13.5 |
| Gap junction | 交叉调控 | 6.00e-05 | 46.6 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 2.06e-04 | 123.0 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 2.18e-04 | 113.9 |
| cellular metal ion homeostasis (GO:0006875) | 其他 | 2.38e-04 | 39.4 |
| response to peptide (GO:1901652) | 其他 | 4.63e-04 | 78.8 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 4.63e-04 | 30.1 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 4.63e-04 | 410.0 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 4.63e-04 | 73.2 |
| Neuroactive ligand-receptor interaction | 其他 | 4.92e-04 | 15.0 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 9.89e-04 | 53.9 |
| icosanoid metabolic process (GO:0006690) | 其他 | 9.89e-04 | 53.9 |
| negative regulation of autophagy (GO:0010507) | 信号主干通路 | 1.02e-03 | 52.1 |
| kynurenine metabolic process (GO:0070189) | 其他 | 1.21e-03 | 227.8 |
| indolalkylamine catabolic process (GO:0046218) | 其他 | 1.32e-03 | 205.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.9_K30</b> — 显著33条, 主干4, 交叉8, 直接OX0, EF=91.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Serotonergic synapse | 交叉调控 | 8.28e-09 | 50.8 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 9.27e-08 | 75.7 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 1.17e-07 | 328.0 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 3.39e-06 | 492.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 8.89e-06 | 351.4 |
| Pathways in cancer | 交叉调控 | 1.05e-05 | 12.4 |
| Neuroactive ligand-receptor interaction | 其他 | 1.27e-04 | 14.4 |
| Gap junction | 交叉调控 | 1.27e-04 | 37.3 |
| Calcium signaling pathway | 信号主干通路 | 2.56e-04 | 17.1 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 4.61e-04 | 98.4 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 4.87e-04 | 91.1 |
| cellular metal ion homeostasis (GO:0006875) | 其他 | 6.77e-04 | 31.5 |
| anterograde trans-synaptic signaling (GO:0098916) | 其他 | 8.34e-04 | 16.8 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 9.05e-04 | 328.0 |
| response to peptide (GO:1901652) | 其他 | 9.05e-04 | 63.1 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 1.03e-03 | 58.6 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 1.14e-03 | 24.1 |
| chemical synaptic transmission (GO:0007268) | 其他 | 1.52e-03 | 13.4 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.87e-03 | 43.2 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.87e-03 | 43.2 |
| ... 还有 13 条通路 ... | | | |

</details>

<details>
<summary><b>P0.15_S0.9_K50</b> — 显著33条, 主干4, 交叉8, 直接OX0, EF=91.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Serotonergic synapse | 交叉调控 | 8.28e-09 | 50.8 |
| G protein-coupled receptor signaling pathway, coupled t | 其他 | 9.27e-08 | 75.7 |
| G protein-coupled serotonin receptor signaling pathway  | 其他 | 1.17e-07 | 328.0 |
| regulation of phosphatidylinositol biosynthetic process | 其他 | 3.39e-06 | 492.0 |
| positive regulation of phospholipid biosynthetic proces | 其他 | 8.89e-06 | 351.4 |
| Pathways in cancer | 交叉调控 | 1.05e-05 | 12.4 |
| Neuroactive ligand-receptor interaction | 其他 | 1.27e-04 | 14.4 |
| Gap junction | 交叉调控 | 1.27e-04 | 37.3 |
| Calcium signaling pathway | 信号主干通路 | 2.56e-04 | 17.1 |
| negative regulation of sequestering of calcium ion (GO: | 其他 | 4.61e-04 | 98.4 |
| release of sequestered calcium ion into cytosol (GO:005 | 其他 | 4.87e-04 | 91.1 |
| cellular metal ion homeostasis (GO:0006875) | 其他 | 6.77e-04 | 31.5 |
| anterograde trans-synaptic signaling (GO:0098916) | 其他 | 8.34e-04 | 16.8 |
| tryptophan catabolic process to kynurenine (GO:0019441) | 其他 | 9.05e-04 | 328.0 |
| response to peptide (GO:1901652) | 其他 | 9.05e-04 | 63.1 |
| calcium ion transmembrane import into cytosol (GO:00975 | 其他 | 1.03e-03 | 58.6 |
| cellular calcium ion homeostasis (GO:0006874) | 其他 | 1.14e-03 | 24.1 |
| chemical synaptic transmission (GO:0007268) | 其他 | 1.52e-03 | 13.4 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.87e-03 | 43.2 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.87e-03 | 43.2 |
| ... 还有 13 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K10</b> — 显著40条, 主干11, 交叉15, 直接OX0, EF=87.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.25e-10 | 70.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.25e-10 | 66.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.82e-10 | 88.5 |
| Hepatitis B | 交叉调控 | 1.87e-09 | 75.9 |
| Human cytomegalovirus infection | 其他 | 1.09e-08 | 54.7 |
| IL-17 signaling pathway | 交叉调控 | 1.12e-08 | 109.0 |
| Chemical carcinogenesis | 交叉调控 | 1.12e-08 | 51.5 |
| Pathways in cancer | 交叉调控 | 1.78e-08 | 27.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.65e-08 | 26.4 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 3.72e-08 | 22.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.72e-08 | 35.5 |
| MicroRNAs in cancer | 交叉调控 | 4.16e-08 | 39.7 |
| Estrogen signaling pathway | 交叉调控 | 4.88e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 8.00e-08 | 29.8 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 8.00e-08 | 29.6 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 9.67e-07 | 33.2 |
| negative regulation of programmed cell death (GO:004306 | 其他 | 9.72e-07 | 32.3 |
| Prostate cancer | 交叉调控 | 1.38e-06 | 84.5 |
| regulation of mitochondrial depolarization (GO:0051900) | 其他 | 1.39e-06 | 559.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.43e-06 | 82.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K20</b> — 显著40条, 主干10, 交叉14, 直接OX0, EF=52.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.20e-19 | 65.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.80e-17 | 57.2 |
| Pathways in cancer | 交叉调控 | 1.43e-16 | 27.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.45e-15 | 73.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.63e-15 | 92.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.34e-14 | 23.1 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.76e-14 | 53.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.80e-14 | 27.6 |
| Human cytomegalovirus infection | 其他 | 1.12e-13 | 45.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.52e-13 | 30.4 |
| Chemical carcinogenesis | 交叉调控 | 1.80e-13 | 42.9 |
| Hepatitis B | 交叉调控 | 3.65e-13 | 56.9 |
| TNF signaling pathway | 交叉调控 | 1.57e-12 | 73.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.72e-12 | 29.0 |
| Estrogen signaling pathway | 交叉调控 | 6.79e-12 | 59.9 |
| Shigellosis | 其他 | 1.12e-11 | 37.5 |
| IL-17 signaling pathway | 交叉调控 | 3.74e-11 | 76.3 |
| Prostate cancer | 交叉调控 | 4.37e-11 | 74.0 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 5.62e-11 | 23.2 |
| C-type lectin receptor signaling pathway | 交叉调控 | 6.75e-11 | 69.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K30</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=41.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.14e-28 | 29.6 |
| Proteoglycans in cancer | 交叉调控 | 7.75e-25 | 56.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.81e-21 | 32.8 |
| Hepatitis B | 交叉调控 | 1.80e-20 | 59.1 |
| Lipid and atherosclerosis | 交叉调控 | 8.43e-19 | 44.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.38e-17 | 75.2 |
| MAPK signaling pathway | 信号主干通路 | 5.10e-17 | 32.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.94e-17 | 27.1 |
| Ras signaling pathway | 信号主干通路 | 1.16e-16 | 38.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 7.16e-16 | 42.5 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.12e-15 | 54.1 |
| Rap1 signaling pathway | 信号主干通路 | 1.64e-15 | 39.0 |
| Prostate cancer | 交叉调控 | 1.94e-15 | 70.4 |
| Human cytomegalovirus infection | 其他 | 3.20e-15 | 36.4 |
| Chemical carcinogenesis | 交叉调控 | 6.17e-15 | 34.3 |
| TNF signaling pathway | 交叉调控 | 6.89e-15 | 61.0 |
| Central carbon metabolism in cancer | 交叉调控 | 7.32e-15 | 87.9 |
| Relaxin signaling pathway | 交叉调控 | 2.60e-14 | 53.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.98e-14 | 21.3 |
| Apoptosis | 信号主干通路 | 6.56e-14 | 48.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K50</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=27.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 2.24e-34 | 23.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.78e-22 | 24.3 |
| Hepatitis B | 交叉调控 | 4.96e-22 | 43.0 |
| Proteoglycans in cancer | 交叉调控 | 4.96e-22 | 36.0 |
| Lipid and atherosclerosis | 交叉调控 | 2.59e-18 | 30.5 |
| Prostate cancer | 交叉调控 | 3.66e-18 | 54.9 |
| MAPK signaling pathway | 信号主干通路 | 8.03e-18 | 23.7 |
| Rap1 signaling pathway | 信号主干通路 | 5.13e-17 | 29.3 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.12e-16 | 19.3 |
| Ras signaling pathway | 信号主干通路 | 2.06e-16 | 26.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.57e-16 | 49.2 |
| Gastric cancer | 交叉调控 | 6.54e-16 | 35.8 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.60e-15 | 15.7 |
| Human T-cell leukemia virus 1 infection | 其他 | 2.60e-15 | 26.2 |
| Relaxin signaling pathway | 交叉调控 | 4.69e-15 | 38.1 |
| MicroRNAs in cancer | 交叉调控 | 1.01e-14 | 19.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.02e-14 | 35.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.36e-14 | 27.6 |
| Epstein-Barr virus infection | 其他 | 2.33e-14 | 26.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.54e-14 | 15.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K10</b> — 显著40条, 主干11, 交叉15, 直接OX0, EF=87.7</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.25e-10 | 70.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.25e-10 | 66.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.82e-10 | 88.5 |
| Hepatitis B | 交叉调控 | 1.87e-09 | 75.9 |
| Human cytomegalovirus infection | 其他 | 1.09e-08 | 54.7 |
| IL-17 signaling pathway | 交叉调控 | 1.12e-08 | 109.0 |
| Chemical carcinogenesis | 交叉调控 | 1.12e-08 | 51.5 |
| Pathways in cancer | 交叉调控 | 1.78e-08 | 27.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.65e-08 | 26.4 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 3.72e-08 | 22.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.72e-08 | 35.5 |
| MicroRNAs in cancer | 交叉调控 | 4.16e-08 | 39.7 |
| Estrogen signaling pathway | 交叉调控 | 4.88e-08 | 74.8 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 8.00e-08 | 29.8 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 8.00e-08 | 29.6 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 9.67e-07 | 33.2 |
| negative regulation of programmed cell death (GO:004306 | 其他 | 9.72e-07 | 32.3 |
| Prostate cancer | 交叉调控 | 1.38e-06 | 84.5 |
| regulation of mitochondrial depolarization (GO:0051900) | 其他 | 1.39e-06 | 559.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.43e-06 | 82.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K20</b> — 显著40条, 主干13, 交叉13, 直接OX0, EF=49.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.75e-22 | 70.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.78e-17 | 57.2 |
| Pathways in cancer | 交叉调控 | 1.42e-16 | 27.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.43e-15 | 73.7 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.61e-15 | 33.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.54e-15 | 92.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.73e-14 | 53.1 |
| Human cytomegalovirus infection | 其他 | 1.11e-13 | 45.6 |
| Chemical carcinogenesis | 交叉调控 | 1.78e-13 | 42.9 |
| Hepatitis B | 交叉调控 | 3.61e-13 | 56.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.48e-13 | 21.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.26e-12 | 25.5 |
| TNF signaling pathway | 交叉调控 | 1.55e-12 | 73.2 |
| MicroRNAs in cancer | 交叉调控 | 1.76e-12 | 33.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.09e-12 | 29.0 |
| Estrogen signaling pathway | 交叉调控 | 6.20e-12 | 59.9 |
| Shigellosis | 其他 | 1.03e-11 | 37.5 |
| IL-17 signaling pathway | 交叉调控 | 3.45e-11 | 76.3 |
| Prostate cancer | 交叉调控 | 4.05e-11 | 74.0 |
| Toxoplasmosis | 其他 | 1.07e-10 | 64.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K30</b> — 显著40条, 主干15, 交叉12, 直接OX0, EF=37.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.26e-26 | 28.3 |
| Proteoglycans in cancer | 交叉调控 | 7.79e-25 | 56.7 |
| Hepatitis B | 交叉调控 | 2.42e-20 | 59.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.63e-19 | 30.9 |
| Lipid and atherosclerosis | 交叉调控 | 8.47e-19 | 44.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.40e-17 | 75.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 9.26e-16 | 42.5 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.41e-15 | 54.1 |
| MAPK signaling pathway | 信号主干通路 | 2.08e-15 | 30.2 |
| Human immunodeficiency virus 1 infection | 其他 | 2.08e-15 | 38.7 |
| Prostate cancer | 交叉调控 | 2.13e-15 | 70.4 |
| MicroRNAs in cancer | 交叉调控 | 3.22e-15 | 28.7 |
| Human cytomegalovirus infection | 其他 | 3.22e-15 | 36.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.25e-15 | 25.4 |
| Ras signaling pathway | 信号主干通路 | 4.33e-15 | 35.3 |
| Chemical carcinogenesis | 交叉调控 | 5.79e-15 | 34.3 |
| TNF signaling pathway | 交叉调控 | 6.50e-15 | 61.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.17e-14 | 21.3 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 2.17e-14 | 21.1 |
| Relaxin signaling pathway | 交叉调控 | 2.62e-14 | 53.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K50</b> — 显著40条, 主干14, 交叉12, 直接OX0, EF=28.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.33e-32 | 23.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 8.30e-24 | 25.5 |
| Hepatitis B | 交叉调控 | 4.89e-22 | 43.0 |
| Proteoglycans in cancer | 交叉调控 | 4.89e-22 | 36.0 |
| Rap1 signaling pathway | 信号主干通路 | 3.42e-20 | 33.2 |
| Ras signaling pathway | 信号主干通路 | 1.60e-19 | 30.0 |
| MAPK signaling pathway | 信号主干通路 | 2.04e-19 | 25.1 |
| Lipid and atherosclerosis | 交叉调控 | 1.60e-18 | 30.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.23e-18 | 20.3 |
| Prostate cancer | 交叉调控 | 1.92e-16 | 50.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.31e-16 | 49.2 |
| Progesterone-mediated oocyte maturation | 其他 | 2.31e-16 | 49.2 |
| MicroRNAs in cancer | 交叉调控 | 3.81e-16 | 21.2 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 3.97e-16 | 29.7 |
| Epstein-Barr virus infection | 其他 | 7.03e-16 | 28.4 |
| Human T-cell leukemia virus 1 infection | 其他 | 2.05e-15 | 26.2 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.47e-15 | 15.7 |
| Relaxin signaling pathway | 交叉调控 | 3.75e-15 | 38.1 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 8.83e-15 | 35.4 |
| Gastric cancer | 交叉调控 | 1.95e-14 | 33.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K10</b> — 显著40条, 主干10, 交叉14, 直接OX0, EF=77.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.25e-10 | 70.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.25e-10 | 66.7 |
| Chemical carcinogenesis | 交叉调控 | 1.76e-10 | 60.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.91e-10 | 29.7 |
| Pathways in cancer | 交叉调控 | 3.57e-10 | 30.9 |
| Hepatitis B | 交叉调控 | 1.49e-09 | 75.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.79e-09 | 34.0 |
| Human cytomegalovirus infection | 其他 | 9.08e-09 | 54.7 |
| IL-17 signaling pathway | 交叉调控 | 1.03e-08 | 109.0 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.82e-08 | 35.5 |
| MicroRNAs in cancer | 交叉调控 | 4.67e-08 | 39.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.25e-08 | 73.7 |
| Estrogen signaling pathway | 交叉调控 | 5.25e-08 | 74.8 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.49e-07 | 53.1 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 1.20e-06 | 33.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.20e-06 | 63.7 |
| Prostate cancer | 交叉调控 | 1.26e-06 | 84.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.32e-06 | 82.0 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 1.32e-06 | 19.3 |
| positive regulation of cellular biosynthetic process (G | 其他 | 1.50e-06 | 56.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K20</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=52.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 1.13e-19 | 65.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.13e-19 | 62.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.13e-17 | 81.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.77e-17 | 102.5 |
| Pathways in cancer | 交叉调控 | 8.56e-17 | 27.0 |
| Hepatitis B | 交叉调控 | 4.61e-15 | 63.3 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.40e-15 | 33.0 |
| TNF signaling pathway | 交叉调控 | 1.56e-14 | 82.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.06e-14 | 53.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.70e-14 | 27.6 |
| IL-17 signaling pathway | 交叉调控 | 4.11e-13 | 87.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 4.80e-13 | 21.5 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 9.82e-13 | 25.4 |
| Toxoplasmosis | 其他 | 1.56e-12 | 73.2 |
| Human immunodeficiency virus 1 infection | 其他 | 3.44e-12 | 43.5 |
| Relaxin signaling pathway | 交叉调控 | 4.14e-12 | 63.6 |
| Human cytomegalovirus infection | 其他 | 4.99e-12 | 41.0 |
| Estrogen signaling pathway | 交叉调控 | 5.79e-12 | 59.9 |
| Chemical carcinogenesis | 交叉调控 | 7.47e-12 | 38.6 |
| Shigellosis | 其他 | 9.09e-12 | 37.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K30</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=37.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.15e-28 | 29.6 |
| Proteoglycans in cancer | 交叉调控 | 9.97e-23 | 53.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.88e-21 | 32.8 |
| Lipid and atherosclerosis | 交叉调控 | 1.15e-20 | 47.7 |
| Hepatitis B | 交叉调控 | 1.60e-18 | 54.8 |
| Rap1 signaling pathway | 信号主干通路 | 3.79e-17 | 42.3 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.79e-17 | 75.2 |
| MAPK signaling pathway | 信号主干通路 | 4.51e-17 | 32.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.84e-17 | 27.1 |
| Ras signaling pathway | 信号主干通路 | 1.04e-16 | 38.3 |
| Relaxin signaling pathway | 交叉调控 | 4.86e-16 | 58.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 5.92e-16 | 42.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.59e-16 | 22.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.44e-16 | 54.1 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.18e-15 | 31.6 |
| Human immunodeficiency virus 1 infection | 其他 | 1.57e-15 | 38.7 |
| Human cytomegalovirus infection | 其他 | 3.00e-15 | 36.4 |
| Chemical carcinogenesis | 交叉调控 | 5.82e-15 | 34.3 |
| TNF signaling pathway | 交叉调控 | 6.53e-15 | 61.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.83e-14 | 17.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K50</b> — 显著40条, 主干14, 交叉12, 直接OX0, EF=34.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.30e-32 | 23.2 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.04e-23 | 23.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.63e-22 | 24.3 |
| Hepatitis B | 交叉调控 | 4.77e-22 | 43.0 |
| Proteoglycans in cancer | 交叉调控 | 4.77e-22 | 36.0 |
| MAPK signaling pathway | 信号主干通路 | 6.55e-21 | 26.5 |
| Rap1 signaling pathway | 信号主干通路 | 2.77e-20 | 33.2 |
| Ras signaling pathway | 信号主干通路 | 1.33e-19 | 30.0 |
| Lipid and atherosclerosis | 交叉调控 | 1.56e-18 | 30.5 |
| Progesterone-mediated oocyte maturation | 其他 | 3.55e-18 | 53.3 |
| Prostate cancer | 交叉调控 | 1.69e-16 | 50.7 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 2.25e-16 | 49.2 |
| MicroRNAs in cancer | 交叉调控 | 3.71e-16 | 21.2 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.26e-15 | 15.7 |
| Human T-cell leukemia virus 1 infection | 其他 | 2.30e-15 | 26.2 |
| Relaxin signaling pathway | 交叉调控 | 4.18e-15 | 38.1 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.76e-15 | 35.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.31e-14 | 27.6 |
| Epstein-Barr virus infection | 其他 | 2.24e-14 | 26.4 |
| TNF signaling pathway | 交叉调控 | 3.43e-14 | 40.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.9_K10</b> — 显著40条, 主干11, 交叉11, 直接OX0, EF=99.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Lipid and atherosclerosis | 交叉调控 | 1.08e-12 | 76.3 |
| Hepatitis B | 交叉调控 | 1.84e-11 | 88.6 |
| Colorectal cancer | 信号主干通路 | 5.70e-11 | 143.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.08e-10 | 123.0 |
| Chemical carcinogenesis | 交叉调控 | 1.16e-10 | 60.0 |
| Pathways in cancer | 交叉调控 | 2.62e-10 | 30.9 |
| Estrogen signaling pathway | 交叉调控 | 4.23e-10 | 89.8 |
| Prolactin signaling pathway | 交叉调控 | 2.20e-09 | 146.4 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 2.62e-09 | 63.7 |
| Pancreatic cancer | 信号主干通路 | 2.62e-09 | 134.9 |
| Focal adhesion | 交叉调控 | 2.62e-09 | 61.2 |
| Viral carcinogenesis | 信号主干通路 | 2.62e-09 | 60.6 |
| Proteoglycans in cancer | 交叉调控 | 2.62e-09 | 60.0 |
| Human immunodeficiency virus 1 infection | 其他 | 2.98e-09 | 58.0 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 3.00e-09 | 33.8 |
| ErbB signaling pathway | 其他 | 3.17e-09 | 120.6 |
| Human cytomegalovirus infection | 其他 | 3.74e-09 | 54.7 |
| Prostate cancer | 交叉调控 | 5.48e-09 | 105.7 |
| Shigellosis | 其他 | 5.51e-09 | 50.0 |
| Progesterone-mediated oocyte maturation | 其他 | 5.51e-09 | 102.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.9_K20</b> — 显著40条, 主干14, 交叉9, 直接OX0, EF=55.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Hepatitis B | 交叉调控 | 1.13e-18 | 75.9 |
| Lipid and atherosclerosis | 交叉调控 | 1.84e-17 | 57.2 |
| Progesterone-mediated oocyte maturation | 其他 | 6.54e-17 | 102.5 |
| Pathways in cancer | 交叉调控 | 1.10e-16 | 27.0 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 7.80e-15 | 92.2 |
| Proteoglycans in cancer | 交叉调控 | 5.20e-14 | 50.0 |
| Human immunodeficiency virus 1 infection | 其他 | 6.26e-14 | 48.3 |
| Chemical carcinogenesis | 交叉调控 | 1.84e-13 | 42.9 |
| MAPK signaling pathway | 信号主干通路 | 1.31e-12 | 34.9 |
| Toxoplasmosis | 其他 | 1.61e-12 | 73.2 |
| Focal adhesion | 交叉调控 | 2.09e-12 | 45.9 |
| Epstein-Barr virus infection | 其他 | 2.09e-12 | 45.7 |
| negative regulation of apoptotic process (GO:0043066) | 信号主干通路 | 3.61e-12 | 25.4 |
| Ras signaling pathway | 信号主干通路 | 6.70e-12 | 39.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.70e-12 | 59.0 |
| Apoptosis | 信号主干通路 | 7.22e-12 | 57.7 |
| Pancreatic cancer | 信号主干通路 | 7.22e-12 | 94.4 |
| Shigellosis | 其他 | 8.80e-12 | 37.5 |
| ErbB signaling pathway | 其他 | 1.44e-11 | 84.4 |
| Colorectal cancer | 信号主干通路 | 1.49e-11 | 83.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.9_K30</b> — 显著40条, 主干16, 交叉10, 直接OX0, EF=42.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Hepatitis B | 交叉调控 | 2.35e-26 | 71.7 |
| Pathways in cancer | 交叉调控 | 4.57e-23 | 25.7 |
| Progesterone-mediated oocyte maturation | 其他 | 3.70e-21 | 88.8 |
| Proteoglycans in cancer | 交叉调控 | 5.34e-21 | 50.0 |
| Lipid and atherosclerosis | 交叉调控 | 8.20e-19 | 44.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.26e-17 | 75.2 |
| Relaxin signaling pathway | 交叉调控 | 6.69e-16 | 58.3 |
| Epstein-Barr virus infection | 其他 | 1.22e-15 | 40.6 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.22e-15 | 54.1 |
| Human immunodeficiency virus 1 infection | 其他 | 1.83e-15 | 38.7 |
| MAPK signaling pathway | 信号主干通路 | 1.83e-15 | 30.2 |
| MicroRNAs in cancer | 交叉调控 | 3.34e-15 | 28.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.73e-14 | 25.1 |
| Focal adhesion | 交叉调控 | 4.85e-14 | 37.4 |
| Rap1 signaling pathway | 信号主干通路 | 7.35e-14 | 35.8 |
| GnRH signaling pathway | 其他 | 1.03e-13 | 66.1 |
| Human T-cell leukemia virus 1 infection | 其他 | 1.03e-13 | 34.3 |
| Prostate cancer | 交叉调控 | 1.36e-13 | 63.4 |
| Cellular senescence | 信号主干通路 | 1.58e-13 | 43.8 |
| Ras signaling pathway | 信号主干通路 | 1.66e-13 | 32.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.9_K50</b> — 显著40条, 主干15, 交叉12, 直接OX0, EF=45.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.63e-29 | 21.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.16e-24 | 24.4 |
| Hepatitis B | 交叉调控 | 8.78e-22 | 43.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.85e-21 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 2.71e-20 | 34.0 |
| MAPK signaling pathway | 信号主干通路 | 2.75e-19 | 25.1 |
| Progesterone-mediated oocyte maturation | 其他 | 5.27e-18 | 53.3 |
| Ras signaling pathway | 信号主干通路 | 6.04e-18 | 28.3 |
| Rap1 signaling pathway | 信号主干通路 | 4.88e-17 | 29.3 |
| Lipid and atherosclerosis | 交叉调控 | 6.19e-17 | 28.6 |
| protein phosphorylation (GO:0006468) | 信号主干通路 | 2.37e-15 | 15.7 |
| Relaxin signaling pathway | 交叉调控 | 5.79e-15 | 38.1 |
| Prostate cancer | 交叉调控 | 1.08e-14 | 46.5 |
| MicroRNAs in cancer | 交叉调控 | 1.12e-14 | 19.8 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 1.29e-14 | 45.1 |
| Fc epsilon RI signaling pathway | 其他 | 1.30e-14 | 60.3 |
| Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 1.38e-14 | 27.6 |
| Epstein-Barr virus infection | 其他 | 2.35e-14 | 26.4 |
| Human immunodeficiency virus 1 infection | 其他 | 4.16e-14 | 25.1 |
| Chemical carcinogenesis | 交叉调控 | 1.86e-13 | 22.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.2_S0.4_K10</b> — 显著31条, 主干2, 交叉10, 直接OX0, EF=117.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 3.62e-12 | 602.9 |
| bicarbonate transport (GO:0015701) | 其他 | 7.50e-10 | 292.9 |
| organic anion transport (GO:0015711) | 其他 | 8.15e-08 | 103.5 |
| response to peptide (GO:1901652) | 其他 | 6.68e-05 | 157.7 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 3.87e-04 | 16.5 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 4.40e-04 | 61.5 |
| cellular response to UV-A (GO:0071492) | 其他 | 6.07e-04 | 372.7 |
| response to UV-A (GO:0070141) | 其他 | 8.37e-04 | 292.9 |
| Proximal tubule bicarbonate reclamation | 其他 | 1.17e-03 | 178.3 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 1.35e-03 | 215.8 |
| Pathways in cancer | 交叉调控 | 1.41e-03 | 15.4 |
| positive regulation of vascular associated smooth muscl | 其他 | 1.45e-03 | 195.2 |
| Proteoglycans in cancer | 交叉调控 | 1.50e-03 | 30.0 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.70e-03 | 17.0 |
| Bladder cancer | 交叉调控 | 1.89e-03 | 100.0 |
| endodermal cell differentiation (GO:0035987) | 其他 | 2.72e-03 | 128.1 |
| endoderm formation (GO:0001706) | 其他 | 3.04e-03 | 113.9 |
| regulation of vascular associated smooth muscle cell pr | 其他 | 3.04e-03 | 110.8 |
| MicroRNAs in cancer | 交叉调控 | 3.62e-03 | 19.8 |
| response to amyloid-beta (GO:1904645) | 其他 | 3.98e-03 | 93.2 |
| ... 还有 11 条通路 ... | | | |

</details>

<details>
<summary><b>P0.2_S0.4_K20</b> — 显著20条, 主干3, 交叉2, 直接OX0, EF=116.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.31e-22 | 542.6 |
| bicarbonate transport (GO:0015701) | 其他 | 2.81e-18 | 263.6 |
| organic anion transport (GO:0015711) | 其他 | 3.33e-14 | 93.2 |
| Pathways in cancer | 交叉调控 | 2.05e-05 | 13.5 |
| response to peptide (GO:1901652) | 其他 | 8.60e-04 | 78.8 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.93e-03 | 12.7 |
| cellular response to UV-A (GO:0071492) | 其他 | 2.93e-03 | 186.4 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 2.93e-03 | 186.4 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.47e-03 | 30.8 |
| response to UV-A (GO:0070141) | 其他 | 4.14e-03 | 146.4 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 4.18e-03 | 136.7 |
| Proximal tubule bicarbonate reclamation | 其他 | 4.75e-03 | 89.1 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 5.44e-03 | 107.9 |
| prostaglandin biosynthetic process (GO:0001516) | 其他 | 5.44e-03 | 107.9 |
| positive regulation of vascular associated smooth muscl | 其他 | 6.06e-03 | 97.6 |
| cellular response to hypoxia (GO:0071456) | 信号主干通路 | 7.77e-03 | 23.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.77e-03 | 8.3 |
| regulation of cellular pH (GO:0030641) | 其他 | 9.78e-03 | 66.1 |
| prostaglandin metabolic process (GO:0006693) | 其他 | 9.78e-03 | 66.1 |
| endodermal cell differentiation (GO:0035987) | 其他 | 9.78e-03 | 64.1 |

</details>

<details>
<summary><b>P0.2_S0.4_K30</b> — 显著12条, 主干0, 交叉3, 直接OX0, EF=97.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.20e-23 | 402.0 |
| bicarbonate transport (GO:0015701) | 其他 | 9.70e-19 | 195.2 |
| organic anion transport (GO:0015711) | 其他 | 3.88e-14 | 69.0 |
| Proximal tubule bicarbonate reclamation | 其他 | 2.83e-04 | 89.1 |
| Pathways in cancer | 交叉调控 | 3.82e-04 | 9.0 |
| response to peptide (GO:1901652) | 其他 | 4.54e-03 | 52.6 |
| Proteoglycans in cancer | 交叉调控 | 6.38e-03 | 13.3 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 8.55e-03 | 36.0 |
| icosanoid metabolic process (GO:0006690) | 其他 | 8.55e-03 | 36.0 |
| cellular response to UV-A (GO:0071492) | 其他 | 8.59e-03 | 124.2 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 8.59e-03 | 124.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 9.55e-03 | 20.5 |

</details>

<details>
<summary><b>P0.2_S0.4_K50</b> — 显著7条, 主干0, 交叉2, 直接OX0, EF=112.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 4.82e-23 | 376.8 |
| bicarbonate transport (GO:0015701) | 其他 | 2.13e-18 | 183.0 |
| organic anion transport (GO:0015711) | 其他 | 7.59e-16 | 71.2 |
| Proximal tubule bicarbonate reclamation | 其他 | 3.51e-04 | 83.6 |
| Pathways in cancer | 交叉调控 | 6.14e-04 | 8.4 |
| response to peptide (GO:1901652) | 其他 | 5.67e-03 | 49.3 |
| Proteoglycans in cancer | 交叉调控 | 8.39e-03 | 12.5 |

</details>

<details>
<summary><b>P0.2_S0.5_K10</b> — 显著33条, 主干2, 交叉10, 直接OX0, EF=123.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 1.81e-15 | 723.5 |
| bicarbonate transport (GO:0015701) | 其他 | 8.82e-13 | 351.4 |
| organic anion transport (GO:0015711) | 其他 | 3.01e-10 | 124.2 |
| response to peptide (GO:1901652) | 其他 | 6.27e-05 | 157.7 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 3.63e-04 | 16.5 |
| cellular response to UV-A (GO:0071492) | 其他 | 5.70e-04 | 372.7 |
| response to UV-A (GO:0070141) | 其他 | 7.86e-04 | 292.9 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 1.26e-03 | 215.8 |
| positive regulation of vascular associated smooth muscl | 其他 | 1.36e-03 | 195.2 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.60e-03 | 17.0 |
| Proximal tubule bicarbonate reclamation | 其他 | 1.75e-03 | 178.3 |
| Pathways in cancer | 交叉调控 | 1.87e-03 | 15.4 |
| Proteoglycans in cancer | 交叉调控 | 1.87e-03 | 30.0 |
| Bladder cancer | 交叉调控 | 2.26e-03 | 100.0 |
| endodermal cell differentiation (GO:0035987) | 其他 | 2.32e-03 | 128.1 |
| regulation of cellular pH (GO:0030641) | 其他 | 2.32e-03 | 132.3 |
| endoderm formation (GO:0001706) | 其他 | 2.45e-03 | 113.9 |
| regulation of intracellular pH (GO:0051453) | 其他 | 2.45e-03 | 110.8 |
| regulation of vascular associated smooth muscle cell pr | 其他 | 2.45e-03 | 110.8 |
| response to amyloid-beta (GO:1904645) | 其他 | 3.24e-03 | 93.2 |
| ... 还有 13 条通路 ... | | | |

</details>

<details>
<summary><b>P0.2_S0.5_K20</b> — 显著13条, 主干1, 交叉2, 直接OX0, EF=129.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.39e-22 | 542.6 |
| bicarbonate transport (GO:0015701) | 其他 | 3.20e-18 | 263.6 |
| organic anion transport (GO:0015711) | 其他 | 3.80e-14 | 93.2 |
| Pathways in cancer | 交叉调控 | 2.13e-05 | 13.5 |
| response to peptide (GO:1901652) | 其他 | 9.80e-04 | 78.8 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.60e-03 | 30.8 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.81e-03 | 12.7 |
| cellular response to UV-A (GO:0071492) | 其他 | 4.00e-03 | 186.4 |
| Proximal tubule bicarbonate reclamation | 其他 | 4.92e-03 | 89.1 |
| response to UV-A (GO:0070141) | 其他 | 5.51e-03 | 146.4 |
| Retrograde endocannabinoid signaling | 其他 | 6.85e-03 | 20.8 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 8.84e-03 | 107.9 |
| positive regulation of vascular associated smooth muscl | 其他 | 9.49e-03 | 97.6 |

</details>

<details>
<summary><b>P0.2_S0.5_K30</b> — 显著12条, 主干0, 交叉3, 直接OX0, EF=97.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.22e-23 | 402.0 |
| bicarbonate transport (GO:0015701) | 其他 | 9.77e-19 | 195.2 |
| organic anion transport (GO:0015711) | 其他 | 3.91e-14 | 69.0 |
| Proximal tubule bicarbonate reclamation | 其他 | 2.85e-04 | 89.1 |
| Pathways in cancer | 交叉调控 | 3.85e-04 | 9.0 |
| response to peptide (GO:1901652) | 其他 | 4.57e-03 | 52.6 |
| Proteoglycans in cancer | 交叉调控 | 6.44e-03 | 13.3 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 8.62e-03 | 36.0 |
| icosanoid metabolic process (GO:0006690) | 其他 | 8.62e-03 | 36.0 |
| cellular response to UV-A (GO:0071492) | 其他 | 8.66e-03 | 124.2 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 8.66e-03 | 124.2 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 9.64e-03 | 20.5 |

</details>

<details>
<summary><b>P0.2_S0.5_K50</b> — 显著11条, 主干0, 交叉2, 直接OX0, EF=101.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 3.28e-23 | 389.0 |
| bicarbonate transport (GO:0015701) | 其他 | 1.46e-18 | 188.9 |
| organic anion transport (GO:0015711) | 其他 | 4.96e-16 | 73.5 |
| Proximal tubule bicarbonate reclamation | 其他 | 3.16e-04 | 86.3 |
| Pathways in cancer | 交叉调控 | 4.86e-04 | 8.7 |
| response to peptide (GO:1901652) | 其他 | 5.10e-03 | 50.9 |
| Proteoglycans in cancer | 交叉调控 | 7.33e-03 | 12.9 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 9.34e-03 | 34.8 |
| icosanoid metabolic process (GO:0006690) | 其他 | 9.34e-03 | 34.8 |
| cellular response to UV-A (GO:0071492) | 其他 | 9.34e-03 | 120.2 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 9.34e-03 | 120.2 |

</details>

<details>
<summary><b>P0.2_S0.7_K10</b> — 显著24条, 主干3, 交叉4, 直接OX0, EF=125.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 5.16e-09 | 482.4 |
| bicarbonate transport (GO:0015701) | 其他 | 4.13e-07 | 234.3 |
| organic anion transport (GO:0015711) | 其他 | 1.46e-05 | 82.8 |
| response to peptide (GO:1901652) | 其他 | 6.84e-05 | 157.7 |
| Pathways in cancer | 交叉调控 | 1.01e-04 | 19.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 3.96e-04 | 16.5 |
| cellular response to UV-A (GO:0071492) | 其他 | 6.22e-04 | 372.7 |
| response to UV-A (GO:0070141) | 其他 | 8.57e-04 | 292.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.13e-03 | 20.3 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 1.21e-03 | 215.8 |
| positive regulation of vascular associated smooth muscl | 其他 | 1.32e-03 | 195.2 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.57e-03 | 17.0 |
| endodermal cell differentiation (GO:0035987) | 其他 | 2.54e-03 | 128.1 |
| Proteoglycans in cancer | 交叉调控 | 2.78e-03 | 30.0 |
| endoderm formation (GO:0001706) | 其他 | 2.88e-03 | 113.9 |
| regulation of vascular associated smooth muscle cell pr | 其他 | 2.88e-03 | 110.8 |
| Bladder cancer | 交叉调控 | 3.15e-03 | 100.0 |
| response to amyloid-beta (GO:1904645) | 其他 | 3.79e-03 | 93.2 |
| positive regulation of smooth muscle cell proliferation | 其他 | 3.87e-03 | 89.1 |
| MicroRNAs in cancer | 交叉调控 | 5.63e-03 | 19.8 |
| ... 还有 4 条通路 ... | | | |

</details>

<details>
<summary><b>P0.2_S0.7_K20</b> — 显著26条, 主干3, 交叉3, 直接OX0, EF=100.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 2.08e-22 | 542.6 |
| bicarbonate transport (GO:0015701) | 其他 | 2.52e-18 | 263.6 |
| organic anion transport (GO:0015711) | 其他 | 2.99e-14 | 93.2 |
| Pathways in cancer | 交叉调控 | 3.46e-04 | 11.6 |
| response to peptide (GO:1901652) | 其他 | 7.71e-04 | 78.8 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.46e-03 | 53.9 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.46e-03 | 53.9 |
| cellular response to UV-A (GO:0071492) | 其他 | 2.25e-03 | 186.4 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 2.25e-03 | 186.4 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 3.12e-03 | 30.8 |
| response to UV-A (GO:0070141) | 其他 | 3.25e-03 | 146.4 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 3.33e-03 | 136.7 |
| Proximal tubule bicarbonate reclamation | 其他 | 4.27e-03 | 89.1 |
| prostaglandin biosynthetic process (GO:0001516) | 其他 | 4.43e-03 | 107.9 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 4.43e-03 | 107.9 |
| positive regulation of vascular associated smooth muscl | 其他 | 4.98e-03 | 97.6 |
| Retrograde endocannabinoid signaling | 其他 | 5.94e-03 | 20.8 |
| cellular response to hypoxia (GO:0071456) | 信号主干通路 | 6.47e-03 | 23.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 6.47e-03 | 8.3 |
| regulation of cellular pH (GO:0030641) | 其他 | 8.22e-03 | 66.1 |
| ... 还有 6 条通路 ... | | | |

</details>

<details>
<summary><b>P0.2_S0.7_K30</b> — 显著15条, 主干1, 交叉2, 直接OX0, EF=99.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 1.24e-20 | 387.6 |
| bicarbonate transport (GO:0015701) | 其他 | 1.65e-16 | 188.3 |
| organic anion transport (GO:0015711) | 其他 | 1.65e-14 | 74.0 |
| Proximal tubule bicarbonate reclamation | 其他 | 2.24e-04 | 95.5 |
| Pathways in cancer | 交叉调控 | 2.28e-04 | 9.7 |
| response to peptide (GO:1901652) | 其他 | 3.54e-03 | 56.3 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 6.68e-03 | 38.5 |
| icosanoid metabolic process (GO:0006690) | 其他 | 6.68e-03 | 38.5 |
| cellular response to UV-A (GO:0071492) | 其他 | 7.21e-03 | 133.1 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 7.21e-03 | 133.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 9.52e-03 | 22.0 |
| response to UV-A (GO:0070141) | 其他 | 9.60e-03 | 104.6 |
| acetyl-CoA metabolic process (GO:0006084) | 其他 | 9.60e-03 | 97.6 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 9.60e-03 | 97.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.77e-03 | 9.1 |

</details>

<details>
<summary><b>P0.2_S0.7_K50</b> — 显著15条, 主干1, 交叉2, 直接OX0, EF=99.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Nitrogen metabolism | 其他 | 1.24e-20 | 387.6 |
| bicarbonate transport (GO:0015701) | 其他 | 1.65e-16 | 188.3 |
| organic anion transport (GO:0015711) | 其他 | 1.65e-14 | 74.0 |
| Proximal tubule bicarbonate reclamation | 其他 | 2.24e-04 | 95.5 |
| Pathways in cancer | 交叉调控 | 2.28e-04 | 9.7 |
| response to peptide (GO:1901652) | 其他 | 3.54e-03 | 56.3 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 6.68e-03 | 38.5 |
| icosanoid metabolic process (GO:0006690) | 其他 | 6.68e-03 | 38.5 |
| cellular response to UV-A (GO:0071492) | 其他 | 7.21e-03 | 133.1 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 7.21e-03 | 133.1 |
| AGE-RAGE signaling pathway in diabetic complications | 交叉调控 | 9.52e-03 | 22.0 |
| response to UV-A (GO:0070141) | 其他 | 9.60e-03 | 104.6 |
| acetyl-CoA metabolic process (GO:0006084) | 其他 | 9.60e-03 | 97.6 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 9.60e-03 | 97.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.77e-03 | 9.1 |

</details>

<details>
<summary><b>P0.2_S0.9_K10</b> — 显著36条, 主干4, 交叉11, 直接OX0, EF=102.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 4.64e-06 | 23.2 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.19e-04 | 25.4 |
| response to peptide (GO:1901652) | 其他 | 1.19e-04 | 157.7 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.90e-04 | 107.9 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.90e-04 | 107.9 |
| cellular response to UV-A (GO:0071492) | 其他 | 6.01e-04 | 372.7 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 6.01e-04 | 372.7 |
| response to UV-A (GO:0070141) | 其他 | 8.51e-04 | 292.9 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 8.59e-04 | 273.3 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 1.12e-03 | 215.8 |
| prostaglandin biosynthetic process (GO:0001516) | 其他 | 1.12e-03 | 215.8 |
| positive regulation of vascular associated smooth muscl | 其他 | 1.25e-03 | 195.2 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.52e-03 | 17.0 |
| Retrograde endocannabinoid signaling | 其他 | 1.63e-03 | 41.6 |
| Proteoglycans in cancer | 交叉调控 | 2.30e-03 | 30.0 |
| Rap1 signaling pathway | 信号主干通路 | 2.30e-03 | 29.3 |
| prostaglandin metabolic process (GO:0006693) | 其他 | 2.31e-03 | 132.3 |
| endodermal cell differentiation (GO:0035987) | 其他 | 2.31e-03 | 128.1 |
| positive regulation of cell motility (GO:2000147) | 其他 | 2.45e-03 | 27.8 |
| regulation of vascular associated smooth muscle cell pr | 其他 | 2.45e-03 | 110.8 |
| ... 还有 16 条通路 ... | | | |

</details>

<details>
<summary><b>P0.2_S0.9_K20</b> — 显著17条, 主干2, 交叉1, 直接OX0, EF=103.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.44e-05 | 15.9 |
| response to peptide (GO:1901652) | 其他 | 1.84e-03 | 92.8 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.84e-03 | 63.5 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.84e-03 | 63.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.84e-03 | 14.9 |
| cellular response to UV-A (GO:0071492) | 其他 | 2.76e-03 | 219.3 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 2.76e-03 | 219.3 |
| response to UV-A (GO:0070141) | 其他 | 3.50e-03 | 172.3 |
| acetyl-CoA metabolic process (GO:0006084) | 其他 | 3.50e-03 | 160.8 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 3.50e-03 | 160.8 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 4.65e-03 | 126.9 |
| prostaglandin biosynthetic process (GO:0001516) | 其他 | 4.65e-03 | 126.9 |
| positive regulation of vascular associated smooth muscl | 其他 | 5.07e-03 | 114.8 |
| fatty acid metabolic process (GO:0006631) | 其他 | 5.07e-03 | 29.2 |
| cellular response to hypoxia (GO:0071456) | 信号主干通路 | 5.54e-03 | 27.6 |
| prostaglandin metabolic process (GO:0006693) | 其他 | 9.22e-03 | 77.8 |
| endodermal cell differentiation (GO:0035987) | 其他 | 9.22e-03 | 75.4 |

</details>

<details>
<summary><b>P0.2_S0.9_K30</b> — 显著17条, 主干2, 交叉1, 直接OX0, EF=103.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.44e-05 | 15.9 |
| response to peptide (GO:1901652) | 其他 | 1.84e-03 | 92.8 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.84e-03 | 63.5 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.84e-03 | 63.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.84e-03 | 14.9 |
| cellular response to UV-A (GO:0071492) | 其他 | 2.76e-03 | 219.3 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 2.76e-03 | 219.3 |
| response to UV-A (GO:0070141) | 其他 | 3.50e-03 | 172.3 |
| acetyl-CoA metabolic process (GO:0006084) | 其他 | 3.50e-03 | 160.8 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 3.50e-03 | 160.8 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 4.65e-03 | 126.9 |
| prostaglandin biosynthetic process (GO:0001516) | 其他 | 4.65e-03 | 126.9 |
| positive regulation of vascular associated smooth muscl | 其他 | 5.07e-03 | 114.8 |
| fatty acid metabolic process (GO:0006631) | 其他 | 5.07e-03 | 29.2 |
| cellular response to hypoxia (GO:0071456) | 信号主干通路 | 5.54e-03 | 27.6 |
| prostaglandin metabolic process (GO:0006693) | 其他 | 9.22e-03 | 77.8 |
| endodermal cell differentiation (GO:0035987) | 其他 | 9.22e-03 | 75.4 |

</details>

<details>
<summary><b>P0.2_S0.9_K50</b> — 显著17条, 主干2, 交叉1, 直接OX0, EF=103.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.44e-05 | 15.9 |
| response to peptide (GO:1901652) | 其他 | 1.84e-03 | 92.8 |
| arachidonic acid metabolic process (GO:0019369) | 其他 | 1.84e-03 | 63.5 |
| icosanoid metabolic process (GO:0006690) | 其他 | 1.84e-03 | 63.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.84e-03 | 14.9 |
| cellular response to UV-A (GO:0071492) | 其他 | 2.76e-03 | 219.3 |
| cyclooxygenase pathway (GO:0019371) | 其他 | 2.76e-03 | 219.3 |
| response to UV-A (GO:0070141) | 其他 | 3.50e-03 | 172.3 |
| acetyl-CoA metabolic process (GO:0006084) | 其他 | 3.50e-03 | 160.8 |
| prostanoid biosynthetic process (GO:0046457) | 其他 | 3.50e-03 | 160.8 |
| regulation of neuroinflammatory response (GO:0150077) | 其他 | 4.65e-03 | 126.9 |
| prostaglandin biosynthetic process (GO:0001516) | 其他 | 4.65e-03 | 126.9 |
| positive regulation of vascular associated smooth muscl | 其他 | 5.07e-03 | 114.8 |
| fatty acid metabolic process (GO:0006631) | 其他 | 5.07e-03 | 29.2 |
| cellular response to hypoxia (GO:0071456) | 信号主干通路 | 5.54e-03 | 27.6 |
| prostaglandin metabolic process (GO:0006693) | 其他 | 9.22e-03 | 77.8 |
| endodermal cell differentiation (GO:0035987) | 其他 | 9.22e-03 | 75.4 |

</details>

## 5. 靶点稳健性分析

| 象限 | 说明 | 靶点数 |
|------|------|--------|
| Q1: 候选稳健 | 高频率 + 高相对连接度 | 50 |
| Q2: 低频高重要 | 低频率 + 高相对连接度 | 12 |
| Q3: 稳健但边缘 | 高频率 + 低相对连接度 | 14 |
| Q4: 可忽略 | 低频率 + 低相对连接度 | 48 |

### 5.1 Q1 候选稳健靶点列表

| 排名 | 基因 | 频率 | 出现次数 | 平均Degree(原始) | 平均相对连接度 | 平均排名 | 平均中介中心性 |
|------|------|------|---------|----------------|-------------------|---------|--------------|
| 1 | STAT3 | 100.0% | 112 | 39.2 | 0.883 | 2.9 | 0.0582 |
| 2 | BCL2 | 85.7% | 96 | 39.0 | 0.740 | 5.6 | 0.0433 |
| 3 | MMP9 | 85.7% | 96 | 31.9 | 0.681 | 8.9 | 0.0523 |
| 4 | PTGS2 | 79.5% | 89 | 30.3 | 0.690 | 9.3 | 0.0683 |
| 5 | SRC | 71.4% | 80 | 58.0 | 0.931 | 2.2 | 0.1301 |
| 6 | EGFR | 71.4% | 80 | 52.5 | 0.852 | 3.1 | 0.0677 |
| 7 | CASP3 | 71.4% | 80 | 42.8 | 0.657 | 7.5 | 0.0426 |
| 8 | MMP2 | 71.4% | 80 | 18.5 | 0.517 | 14.8 | 0.0112 |
| 9 | HSP90AA1 | 67.9% | 76 | 38.9 | 0.684 | 7.1 | 0.0473 |
| 10 | TNF | 64.3% | 72 | 59.6 | 0.829 | 5.0 | 0.0945 |
| 11 | MAPK1 | 58.0% | 65 | 28.7 | 0.556 | 11.5 | 0.0334 |
| 12 | ESR1 | 57.1% | 64 | 46.5 | 0.696 | 6.9 | 0.0569 |
| 13 | FGF2 | 55.4% | 62 | 28.9 | 0.513 | 13.4 | 0.0194 |
| 14 | PIK3CA | 54.5% | 61 | 25.7 | 0.511 | 14.5 | 0.0241 |
| 15 | CDK1 | 54.5% | 61 | 12.4 | 0.313 | 22.3 | 0.0140 |
| 16 | HIF1A | 53.6% | 60 | 45.9 | 0.631 | 8.8 | 0.0330 |
| 17 | KDR | 51.8% | 58 | 31.1 | 0.475 | 15.2 | 0.0171 |
| 18 | MAPK14 | 47.3% | 53 | 24.5 | 0.421 | 18.9 | 0.0126 |
| 19 | MAPK8 | 46.4% | 52 | 20.0 | 0.469 | 16.7 | 0.0202 |
| 20 | APP | 42.9% | 48 | 30.8 | 0.466 | 16.4 | 0.0529 |
| 21 | BCL2L1 | 42.9% | 48 | 25.4 | 0.405 | 20.3 | 0.0054 |
| 22 | MTOR | 36.6% | 41 | 36.3 | 0.446 | 19.3 | 0.0138 |
| 23 | ERBB2 | 35.7% | 40 | 43.0 | 0.554 | 14.0 | 0.0264 |
| 24 | ICAM1 | 35.7% | 40 | 35.9 | 0.437 | 18.7 | 0.0075 |
| 25 | PARP1 | 35.7% | 40 | 26.1 | 0.377 | 23.8 | 0.0176 |
| 26 | STAT1 | 34.8% | 39 | 30.4 | 0.438 | 19.4 | 0.0114 |
| 27 | PRKCA | 34.8% | 39 | 20.1 | 0.349 | 24.0 | 0.0387 |
| 28 | MET | 29.5% | 33 | 19.8 | 0.363 | 26.4 | 0.0036 |
| 29 | IL2 | 28.6% | 32 | 36.4 | 0.425 | 20.9 | 0.0153 |
| 30 | KIT | 26.8% | 30 | 25.8 | 0.349 | 27.2 | 0.0131 |
| 31 | CA9 | 25.0% | 28 | 5.1 | 0.435 | 13.4 | 0.0990 |
| 32 | CDK2 | 24.1% | 27 | 19.9 | 0.325 | 30.2 | 0.0040 |
| 33 | CCNA2 | 23.2% | 26 | 19.2 | 0.319 | 31.2 | 0.0032 |
| 34 | CA1 | 21.4% | 24 | 5.7 | 0.658 | 5.3 | 0.0666 |
| 35 | MAPT | 21.4% | 24 | 17.0 | 0.313 | 29.1 | 0.0195 |
| 36 | CA2 | 20.5% | 23 | 4.5 | 0.463 | 7.4 | 0.0415 |
| 37 | HTR1A | 20.5% | 23 | 7.0 | 0.371 | 20.6 | 0.0551 |
| 38 | HTR2A | 19.6% | 22 | 9.0 | 0.519 | 20.8 | 0.0178 |
| 39 | HDAC1 | 19.6% | 22 | 19.2 | 0.365 | 21.8 | 0.0139 |
| 40 | MAPK9 | 19.6% | 22 | 9.6 | 0.357 | 22.5 | 0.0033 |
| 41 | FAAH | 19.6% | 22 | 1.8 | 0.306 | 15.5 | 0.0382 |
| 42 | CA4 | 18.8% | 21 | 4.5 | 0.471 | 10.5 | 0.0362 |
| 43 | CA3 | 17.9% | 20 | 3.3 | 0.417 | 11.1 | 0.0034 |
| 44 | ESR2 | 17.9% | 20 | 12.9 | 0.296 | 30.8 | 0.0046 |
| 45 | CA7 | 17.0% | 19 | 2.7 | 0.300 | 14.3 | 0.0005 |
| 46 | CA12 | 17.0% | 19 | 2.7 | 0.297 | 13.1 | 0.0035 |
| 47 | NTRK1 | 16.1% | 18 | 23.1 | 0.300 | 35.5 | 0.0119 |
| 48 | ABCB1 | 14.3% | 16 | 13.9 | 0.301 | 27.8 | 0.0237 |
| 49 | CCNA1 | 14.3% | 16 | 11.3 | 0.294 | 32.9 | 0.0009 |
| 50 | IDO1 | 13.4% | 15 | 6.6 | 0.307 | 24.7 | 0.0210 |

## 6. 通路排名稳定性 (Top 25)

| 排名 | 通路 | 类别 | 频率 | 平均排名 | 排名波动 | 平均logP | 平均EF |
|------|------|------|------|---------|---------|---------|-------|
| 1 | Nitrogen metabolism | 其他 | 20.5% | 1.0 | 0.0 | 17.6 | 426.3 |
| 2 | bicarbonate transport (GO:0015701) | 其他 | 20.5% | 2.0 | 0.2 | 13.9 | 207.1 |
| 3 | organic anion transport (GO:0015711) | 其他 | 20.5% | 3.5 | 2.1 | 11.0 | 75.1 |
| 3 | Pathways in cancer | 交叉调控 | 100.0% | 3.6 | 3.6 | 16.4 | 23.5 |
| 4 | G protein-coupled serotonin receptor signaling pat | 其他 | 12.5% | 4.1 | 2.4 | 6.4 | 347.9 |
| 5 | G protein-coupled receptor signaling pathway, coup | 其他 | 12.5% | 5.4 | 3.3 | 5.8 | 71.5 |
| 5 | Proteoglycans in cancer | 交叉调控 | 87.5% | 5.4 | 8.2 | 16.3 | 51.8 |
| 6 | regulation of phosphatidylinositol biosynthetic pr | 其他 | 11.6% | 6.0 | 1.7 | 5.2 | 526.5 |
| 7 | positive regulation of phospholipid biosynthetic p | 其他 | 11.6% | 7.7 | 2.4 | 4.9 | 376.0 |
| 7 | Lipid and atherosclerosis | 交叉调控 | 78.6% | 7.9 | 7.8 | 13.4 | 44.9 |
| 8 | Hepatitis B | 交叉调控 | 80.4% | 8.3 | 7.7 | 13.5 | 54.3 |
| 8 | response to peptide (GO:1901652) | 其他 | 28.6% | 8.5 | 5.5 | 3.1 | 87.9 |
| 8 | Serotonergic synapse | 交叉调控 | 15.2% | 8.6 | 10.0 | 5.6 | 44.6 |
| 8 | MAPK signaling pathway | 信号主干通路 | 34.8% | 8.7 | 5.1 | 16.5 | 28.0 |
| 9 | PI3K-Akt signaling pathway | 信号主干通路 | 41.1% | 9.7 | 7.2 | 16.2 | 25.1 |
| 10 | Progesterone-mediated oocyte maturation | 其他 | 13.4% | 10.5 | 6.6 | 15.6 | 66.3 |
| 10 | AGE-RAGE signaling pathway in diabetic complicatio | 交叉调控 | 94.6% | 10.6 | 6.2 | 10.7 | 66.7 |
| 10 | acetyl-CoA metabolic process (GO:0006084) | 其他 | 4.5% | 10.6 | 2.2 | 2.3 | 135.5 |
| 11 | Neuroactive ligand-receptor interaction | 其他 | 10.7% | 11.1 | 3.7 | 3.6 | 14.8 |
| 11 | prostanoid biosynthetic process (GO:0046457) | 其他 | 7.1% | 11.2 | 1.9 | 2.4 | 153.0 |
| 11 | Kaposi sarcoma-associated herpesvirus infection | 交叉调控 | 69.6% | 11.7 | 8.3 | 12.7 | 44.7 |
| 12 | Colorectal cancer | 信号主干通路 | 8.9% | 12.0 | 6.4 | 10.2 | 114.4 |
| 12 | Rap1 signaling pathway | 信号主干通路 | 25.9% | 12.0 | 5.4 | 15.3 | 31.7 |
| 12 | icosanoid metabolic process (GO:0006690) | 其他 | 14.3% | 12.2 | 7.2 | 2.6 | 49.0 |
| 12 | prostaglandin biosynthetic process (GO:0001516) | 其他 | 5.4% | 12.3 | 1.5 | 2.4 | 135.4 |

### 6.1 通路稳健性筛选

通路稳健性采用与靶点层同构的二维筛选逻辑：跨参数组合出现频率表示稳健性，平均 -log10(adjusted p-value) 表示富集强度。

- 频率中位数: 0.0804
- 平均 -log10(adjusted p-value) 中位数: 7.3744
- 绝对稳健通路: 1 条
- Q1 候选稳健通路: 71 条

| 象限 | 说明 | 通路数 | 平均频率 | 平均-log10P |
|------|------|-------:|---------:|-----------:|
| Q1: 候选稳健通路 | 高频率 + 高富集强度 | 71 | 0.371 | 10.84 |
| Q2: 低频高显著 | 低频率 + 高富集强度 | 36 | 0.031 | 10.80 |
| Q3: 稳健但效应弱 | 高频率 + 低富集强度 | 43 | 0.149 | 4.40 |
| Q4: 可忽略 | 低频率 + 低富集强度 | 64 | 0.036 | 3.63 |

### 6.2 绝对稳健通路

| 通路 | 来源 | 类别 | 频率 | 平均-log10P | 平均排名 |
|------|------|------|------|-------------|---------|
| Pathways in cancer | KEGG | 交叉调控 | 100.0% | 16.37 | 3.6 |

### 6.3 Q1 候选稳健通路 Top 20

| 排名 | 通路 | 来源 | 类别 | 频率 | 平均-log10P | 平均排名 | 平均EF |
|------|------|------|------|------|-------------|---------|-------|
| 1 | Pathways in cancer | KEGG | 交叉调控 | 100.0% | 16.37 | 3.6 | 23.5 |
| 2 | AGE-RAGE signaling pathway in diabetic complications | KEGG | 交叉调控 | 94.6% | 10.73 | 10.6 | 66.7 |
| 3 | Proteoglycans in cancer | KEGG | 交叉调控 | 87.5% | 16.30 | 5.4 | 51.8 |
| 4 | Hepatitis B | KEGG | 交叉调控 | 80.4% | 13.50 | 8.3 | 54.3 |
| 5 | transmembrane receptor protein tyrosine kinase signalin | GO_BP | 信号主干通路 | 80.4% | 11.54 | 13.0 | 25.1 |
| 6 | cellular response to cytokine stimulus (GO:0071345) | GO_BP | 信号主干通路 | 79.5% | 10.34 | 16.0 | 21.9 |
| 7 | Lipid and atherosclerosis | KEGG | 交叉调控 | 78.6% | 13.35 | 7.9 | 44.9 |
| 8 | cytokine-mediated signaling pathway (GO:0019221) | GO_BP | 信号主干通路 | 78.6% | 9.69 | 17.7 | 18.0 |
| 9 | Chemical carcinogenesis | KEGG | 交叉调控 | 74.1% | 11.67 | 12.9 | 38.5 |
| 10 | Fluid shear stress and atherosclerosis | KEGG | 交叉调控 | 70.5% | 10.02 | 15.6 | 51.1 |
| 11 | Kaposi sarcoma-associated herpesvirus infection | KEGG | 交叉调控 | 69.6% | 12.71 | 11.7 | 44.7 |
| 12 | MicroRNAs in cancer | KEGG | 交叉调控 | 69.6% | 9.76 | 16.1 | 24.7 |
| 13 | negative regulation of apoptotic process (GO:0043066) | GO_BP | 信号主干通路 | 68.8% | 9.77 | 25.0 | 20.1 |
| 14 | positive regulation of protein kinase B signaling (GO:0 | GO_BP | 其他 | 66.1% | 9.14 | 26.7 | 44.1 |
| 15 | positive regulation of intracellular signal transductio | GO_BP | 其他 | 64.3% | 9.41 | 25.4 | 17.7 |
| 16 | regulation of phosphatidylinositol 3-kinase signaling ( | GO_BP | 其他 | 60.7% | 9.21 | 29.9 | 56.3 |
| 17 | regulation of protein kinase B signaling (GO:0051896) | GO_BP | 其他 | 58.9% | 9.55 | 26.8 | 35.4 |
| 18 | Estrogen signaling pathway | KEGG | 交叉调控 | 58.0% | 8.54 | 15.9 | 57.3 |
| 19 | positive regulation of phosphorylation (GO:0042327) | GO_BP | 信号主干通路 | 56.2% | 10.22 | 24.6 | 30.0 |
| 20 | Prostate cancer | KEGG | 交叉调控 | 52.7% | 11.41 | 15.3 | 67.5 |

## 7. Pvalue 阈值敏感性对比

| 指标 | P=0.01 | P=0.05 |
|------|--------|--------|
| 显著通路总数 | 40 | 40 |
| 主干道 | 14 | 14 |
| 交叉调控 | 13 | 13 |
| 直接OX | 1 | 1 |
| 其他 | 12 | 12 |
| 平均EF | 33.92 | 33.92 |
| 共有通路 | 40 | 40 |
| 分类分布一致 | 是（Pvalue 稳健） | 是（Pvalue 稳健） |

## 8. 最终判定

| 检查项 | 结果 | 判定 |
|--------|------|------|
| Q1候选稳健靶点数>=5 | 50 | PASS |
| 关键氧化应激靶点进入Q1 | 8个(TNF,SRC,EGFR,STAT3,BCL2,MMP9,HIF1A,PTGS2) | PASS |
| 排名稳定通路数>=3 | 4 | PASS |
| 信号网络主干道排名稳定 | 0个 | FAIL |
| Q1候选稳健通路数>=5 | 71 | PASS |

**通过: 4/5**

> ✅ **判定: 敏感性分析通过。结果支持故事成立。**