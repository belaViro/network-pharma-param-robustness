# 参数网格敏感性分析报告

**生成时间**: 2026-04-27 10:58:23

## 1. 参数网格设置

| 参数 | 取值 | 组合数 |
|------|------|--------|
| Prob (SwissTarget) | [0, 0.09, 0.1] | 3 |
| STRING (PPI) | [0.3, 0.35, 0.4, 0.45, 0.5, 0.7] | 6 |
| TopK (Hub基因) | [25, 30, 35, 40] | 4 |
| Pvalue 固定 | 0.01 | 1 |
| **总计** | | **72** |

## 2. 数据加载概况

- SwissTarget 并集: 288 靶点
- 疾病基因: 14765 个
- 疾病交集后: 275 靶点
- STRING PPI 边数: 830
- 成功组合: 72 / 72

## 3. 各概率阈值靶点数

| Prob 阈值 | 疾病相关靶点数 |
|-----------|---------------|
| ≥ 0 | 275 |
| ≥ 0.09 | 107 |
| ≥ 0.1 | 69 |

## 4. 每组参数富集分析详情

| 组合 | 靶点数 | PPI节点 | PPI边 | 显著通路 | 主干道 | 交叉调控 | 直接OX | 其他 | 平均EF | 最小P |
|------|--------|--------|-------|---------|-------|---------|--------|-------|--------|-------|
| P0_S0.3_K25 | 112 | 110 | 815 | 40 | 12 | 16 | 0 | 12 | 41.1 | 6.51e-18 |
| P0_S0.3_K30 | 112 | 110 | 815 | 40 | 9 | 17 | 0 | 14 | 47.2 | 3.47e-19 |
| P0_S0.3_K35 | 112 | 110 | 815 | 40 | 10 | 17 | 0 | 13 | 40.9 | 3.84e-19 |
| P0_S0.3_K40 | 112 | 110 | 815 | 40 | 14 | 15 | 0 | 11 | 32.2 | 3.08e-19 |
| P0_S0.35_K25 | 112 | 110 | 815 | 40 | 12 | 16 | 0 | 12 | 41.1 | 6.51e-18 |
| P0_S0.35_K30 | 112 | 110 | 815 | 40 | 9 | 17 | 0 | 14 | 47.2 | 3.47e-19 |
| P0_S0.35_K35 | 112 | 110 | 815 | 40 | 10 | 17 | 0 | 13 | 40.9 | 3.84e-19 |
| P0_S0.35_K40 | 112 | 110 | 815 | 40 | 14 | 15 | 0 | 11 | 32.2 | 3.08e-19 |
| P0_S0.4_K25 | 112 | 110 | 815 | 40 | 12 | 16 | 0 | 12 | 41.1 | 6.51e-18 |
| P0_S0.4_K30 | 112 | 110 | 815 | 40 | 9 | 17 | 0 | 14 | 47.2 | 3.47e-19 |
| P0_S0.4_K35 | 112 | 110 | 815 | 40 | 10 | 17 | 0 | 13 | 40.9 | 3.84e-19 |
| P0_S0.4_K40 | 112 | 110 | 815 | 40 | 14 | 15 | 0 | 11 | 32.2 | 3.08e-19 |
| P0_S0.45_K25 | 112 | 110 | 681 | 40 | 11 | 15 | 0 | 14 | 53.6 | 5.89e-20 |
| P0_S0.45_K30 | 112 | 110 | 681 | 40 | 15 | 14 | 0 | 11 | 43.8 | 3.50e-19 |
| P0_S0.45_K35 | 112 | 110 | 681 | 40 | 14 | 15 | 0 | 11 | 34.5 | 3.86e-19 |
| P0_S0.45_K40 | 112 | 110 | 681 | 40 | 12 | 15 | 0 | 13 | 29.4 | 7.64e-21 |
| P0_S0.5_K25 | 112 | 107 | 573 | 40 | 10 | 17 | 0 | 13 | 57.1 | 1.97e-21 |
| P0_S0.5_K30 | 112 | 107 | 573 | 40 | 14 | 15 | 0 | 11 | 46.5 | 5.67e-21 |
| P0_S0.5_K35 | 112 | 107 | 573 | 40 | 12 | 15 | 0 | 13 | 33.3 | 8.01e-21 |
| P0_S0.5_K40 | 112 | 107 | 573 | 40 | 12 | 15 | 0 | 13 | 30.0 | 7.64e-21 |
| P0_S0.7_K25 | 112 | 98 | 285 | 40 | 15 | 14 | 0 | 11 | 41.0 | 7.70e-20 |
| P0_S0.7_K30 | 112 | 98 | 285 | 40 | 13 | 16 | 0 | 11 | 35.1 | 5.36e-21 |
| P0_S0.7_K35 | 112 | 98 | 285 | 40 | 12 | 15 | 0 | 13 | 77.6 | 3.46e-24 |
| P0_S0.7_K40 | 112 | 98 | 285 | 40 | 10 | 15 | 0 | 15 | 74.6 | 1.88e-22 |
| P0.09_S0.3_K25 | 107 | 106 | 790 | 40 | 12 | 16 | 0 | 12 | 41.1 | 6.51e-18 |
| P0.09_S0.3_K30 | 107 | 106 | 790 | 40 | 9 | 17 | 0 | 14 | 47.2 | 3.47e-19 |
| P0.09_S0.3_K35 | 107 | 106 | 790 | 40 | 10 | 17 | 0 | 13 | 40.9 | 3.84e-19 |
| P0.09_S0.3_K40 | 107 | 106 | 790 | 40 | 14 | 15 | 0 | 11 | 32.2 | 3.08e-19 |
| P0.09_S0.35_K25 | 107 | 106 | 790 | 40 | 12 | 16 | 0 | 12 | 41.1 | 6.51e-18 |
| P0.09_S0.35_K30 | 107 | 106 | 790 | 40 | 9 | 17 | 0 | 14 | 47.2 | 3.47e-19 |
| P0.09_S0.35_K35 | 107 | 106 | 790 | 40 | 10 | 17 | 0 | 13 | 40.9 | 3.84e-19 |
| P0.09_S0.35_K40 | 107 | 106 | 790 | 40 | 14 | 15 | 0 | 11 | 32.2 | 3.08e-19 |
| P0.09_S0.4_K25 | 107 | 106 | 790 | 40 | 12 | 16 | 0 | 12 | 41.1 | 6.51e-18 |
| P0.09_S0.4_K30 | 107 | 106 | 790 | 40 | 9 | 17 | 0 | 14 | 47.2 | 3.47e-19 |
| P0.09_S0.4_K35 | 107 | 106 | 790 | 40 | 10 | 17 | 0 | 13 | 40.9 | 3.84e-19 |
| P0.09_S0.4_K40 | 107 | 106 | 790 | 40 | 14 | 15 | 0 | 11 | 32.2 | 3.08e-19 |
| P0.09_S0.45_K25 | 107 | 106 | 659 | 40 | 11 | 15 | 0 | 14 | 53.6 | 5.89e-20 |
| P0.09_S0.45_K30 | 107 | 106 | 659 | 40 | 15 | 14 | 0 | 11 | 43.8 | 3.50e-19 |
| P0.09_S0.45_K35 | 107 | 106 | 659 | 40 | 14 | 15 | 0 | 11 | 34.5 | 3.86e-19 |
| P0.09_S0.45_K40 | 107 | 106 | 659 | 40 | 12 | 15 | 0 | 13 | 29.4 | 7.64e-21 |
| P0.09_S0.5_K25 | 107 | 103 | 556 | 40 | 10 | 17 | 0 | 13 | 57.1 | 1.97e-21 |
| P0.09_S0.5_K30 | 107 | 103 | 556 | 40 | 14 | 15 | 0 | 11 | 46.5 | 5.67e-21 |
| P0.09_S0.5_K35 | 107 | 103 | 556 | 40 | 12 | 15 | 0 | 13 | 33.3 | 8.01e-21 |
| P0.09_S0.5_K40 | 107 | 103 | 556 | 40 | 12 | 15 | 0 | 13 | 30.0 | 7.64e-21 |
| P0.09_S0.7_K25 | 107 | 94 | 274 | 40 | 15 | 14 | 0 | 11 | 41.0 | 7.70e-20 |
| P0.09_S0.7_K30 | 107 | 94 | 274 | 40 | 13 | 16 | 0 | 11 | 35.1 | 5.36e-21 |
| P0.09_S0.7_K35 | 107 | 94 | 274 | 40 | 12 | 15 | 0 | 13 | 77.6 | 3.46e-24 |
| P0.09_S0.7_K40 | 107 | 94 | 274 | 40 | 10 | 15 | 0 | 15 | 74.6 | 1.88e-22 |
| P0.1_S0.3_K25 | 69 | 68 | 355 | 40 | 15 | 15 | 0 | 10 | 35.4 | 3.88e-14 |
| P0.1_S0.3_K30 | 69 | 68 | 355 | 40 | 14 | 14 | 0 | 12 | 31.9 | 1.90e-13 |
| P0.1_S0.3_K35 | 69 | 68 | 355 | 40 | 13 | 14 | 0 | 13 | 29.0 | 2.87e-12 |
| P0.1_S0.3_K40 | 69 | 68 | 355 | 40 | 13 | 14 | 0 | 13 | 27.5 | 2.73e-11 |
| P0.1_S0.35_K25 | 69 | 68 | 355 | 40 | 15 | 15 | 0 | 10 | 35.4 | 3.88e-14 |
| P0.1_S0.35_K30 | 69 | 68 | 355 | 40 | 14 | 14 | 0 | 12 | 31.9 | 1.90e-13 |
| P0.1_S0.35_K35 | 69 | 68 | 355 | 40 | 13 | 14 | 0 | 13 | 29.0 | 2.87e-12 |
| P0.1_S0.35_K40 | 69 | 68 | 355 | 40 | 13 | 14 | 0 | 13 | 27.5 | 2.73e-11 |
| P0.1_S0.4_K25 | 69 | 68 | 355 | 40 | 15 | 15 | 0 | 10 | 35.4 | 3.88e-14 |
| P0.1_S0.4_K30 | 69 | 68 | 355 | 40 | 14 | 14 | 0 | 12 | 31.9 | 1.90e-13 |
| P0.1_S0.4_K35 | 69 | 68 | 355 | 40 | 13 | 14 | 0 | 13 | 29.0 | 2.87e-12 |
| P0.1_S0.4_K40 | 69 | 68 | 355 | 40 | 13 | 14 | 0 | 13 | 27.5 | 2.73e-11 |
| P0.1_S0.45_K25 | 69 | 68 | 301 | 40 | 14 | 15 | 0 | 11 | 38.9 | 5.99e-15 |
| P0.1_S0.45_K30 | 69 | 68 | 301 | 40 | 14 | 15 | 0 | 11 | 29.8 | 1.90e-13 |
| P0.1_S0.45_K35 | 69 | 68 | 301 | 40 | 14 | 14 | 0 | 12 | 32.4 | 7.82e-14 |
| P0.1_S0.45_K40 | 69 | 68 | 301 | 40 | 13 | 14 | 0 | 13 | 29.9 | 9.25e-13 |
| P0.1_S0.5_K25 | 69 | 66 | 254 | 40 | 14 | 15 | 0 | 11 | 38.9 | 5.98e-15 |
| P0.1_S0.5_K30 | 69 | 66 | 254 | 40 | 14 | 15 | 0 | 11 | 29.1 | 4.04e-15 |
| P0.1_S0.5_K35 | 69 | 66 | 254 | 40 | 13 | 14 | 0 | 13 | 60.3 | 4.28e-17 |
| P0.1_S0.5_K40 | 69 | 66 | 254 | 40 | 13 | 14 | 0 | 13 | 53.2 | 7.98e-16 |
| P0.1_S0.7_K25 | 69 | 60 | 136 | 40 | 13 | 14 | 0 | 13 | 103.5 | 9.16e-23 |
| P0.1_S0.7_K30 | 69 | 60 | 136 | 40 | 14 | 14 | 0 | 12 | 85.2 | 1.65e-20 |
| P0.1_S0.7_K35 | 69 | 60 | 136 | 40 | 15 | 14 | 0 | 11 | 73.8 | 8.23e-19 |
| P0.1_S0.7_K40 | 69 | 60 | 136 | 40 | 15 | 13 | 0 | 12 | 67.2 | 1.93e-17 |

### 4.1 各组合富集通路明细

<details>
<summary><b>P0.09_S0.35_K25</b> — 显著40条, 主干12, 交叉16, 直接OX0, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.51e-18 | 52.0 |
| Pathways in cancer | 交叉调控 | 6.51e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 6.72e-14 | 38.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.86e-13 | 26.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.40e-13 | 18.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.11e-13 | 22.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.62e-12 | 19.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.00e-11 | 45.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.36e-10 | 47.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.80e-10 | 23.2 |
| Bladder cancer | 交叉调控 | 2.84e-10 | 120.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-10 | 35.7 |
| MAPK signaling pathway | 信号主干通路 | 1.28e-09 | 25.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.79e-09 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 2.95e-09 | 30.5 |
| Estrogen signaling pathway | 交叉调控 | 4.76e-09 | 41.9 |
| Chemical carcinogenesis | 交叉调控 | 5.47e-09 | 27.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.49e-09 | 54.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.49e-09 | 54.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.19e-09 | 87.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.35_K30</b> — 显著40条, 主干9, 交叉17, 直接OX0, EF=47.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.47e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.47e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.99e-14 | 35.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.05e-13 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.21e-12 | 36.3 |
| positive regulation of intracellular signal transductio | 其他 | 3.40e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.81e-12 | 42.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.05e-11 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-11 | 18.4 |
| MAPK signaling pathway | 信号主干通路 | 1.36e-11 | 25.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.18e-11 | 21.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.35e-10 | 58.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 2.28e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.72e-10 | 51.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.72e-10 | 51.6 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 6.05e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.35_K35</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=40.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.84e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.46e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.47e-13 | 30.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.64e-12 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.72e-11 | 31.1 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-11 | 36.4 |
| positive regulation of intracellular signal transductio | 其他 | 4.09e-11 | 15.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 6.59e-11 | 46.2 |
| MAPK signaling pathway | 信号主干通路 | 8.10e-11 | 21.9 |
| Rap1 signaling pathway | 信号主干通路 | 8.10e-11 | 27.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-10 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.47e-10 | 13.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.41e-10 | 18.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 5.70e-10 | 49.8 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.97e-10 | 73.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.14e-09 | 44.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.14e-09 | 44.2 |
| Bladder cancer | 交叉调控 | 1.89e-09 | 85.7 |
| Estrogen signaling pathway | 交叉调控 | 1.89e-09 | 34.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.89e-09 | 33.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.35_K40</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=32.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.08e-19 | 19.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.70e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-14 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.26e-12 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.39e-12 | 35.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.33e-11 | 26.8 |
| MAPK signaling pathway | 信号主干通路 | 1.35e-11 | 20.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.77e-11 | 14.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.67e-11 | 17.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.32e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.94e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.57e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.96e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 5.51e-10 | 20.6 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.23e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.10e-10 | 11.6 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.20e-09 | 43.6 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.20e-09 | 18.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.3_K25</b> — 显著40条, 主干12, 交叉16, 直接OX0, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.51e-18 | 52.0 |
| Pathways in cancer | 交叉调控 | 6.51e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 6.72e-14 | 38.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.86e-13 | 26.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.40e-13 | 18.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.11e-13 | 22.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.62e-12 | 19.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.00e-11 | 45.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.36e-10 | 47.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.80e-10 | 23.2 |
| Bladder cancer | 交叉调控 | 2.84e-10 | 120.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-10 | 35.7 |
| MAPK signaling pathway | 信号主干通路 | 1.28e-09 | 25.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.79e-09 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 2.95e-09 | 30.5 |
| Estrogen signaling pathway | 交叉调控 | 4.76e-09 | 41.9 |
| Chemical carcinogenesis | 交叉调控 | 5.47e-09 | 27.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.49e-09 | 54.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.49e-09 | 54.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.19e-09 | 87.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.3_K30</b> — 显著40条, 主干9, 交叉17, 直接OX0, EF=47.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.47e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.47e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.99e-14 | 35.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.05e-13 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.21e-12 | 36.3 |
| positive regulation of intracellular signal transductio | 其他 | 3.40e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.81e-12 | 42.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.05e-11 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-11 | 18.4 |
| MAPK signaling pathway | 信号主干通路 | 1.36e-11 | 25.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.18e-11 | 21.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.35e-10 | 58.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 2.28e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.72e-10 | 51.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.72e-10 | 51.6 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 6.05e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.3_K35</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=40.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.84e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.46e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.47e-13 | 30.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.64e-12 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.72e-11 | 31.1 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-11 | 36.4 |
| positive regulation of intracellular signal transductio | 其他 | 4.09e-11 | 15.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 6.59e-11 | 46.2 |
| MAPK signaling pathway | 信号主干通路 | 8.10e-11 | 21.9 |
| Rap1 signaling pathway | 信号主干通路 | 8.10e-11 | 27.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-10 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.47e-10 | 13.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.41e-10 | 18.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 5.70e-10 | 49.8 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.97e-10 | 73.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.14e-09 | 44.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.14e-09 | 44.2 |
| Bladder cancer | 交叉调控 | 1.89e-09 | 85.7 |
| Estrogen signaling pathway | 交叉调控 | 1.89e-09 | 34.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.89e-09 | 33.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.3_K40</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=32.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.08e-19 | 19.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.70e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-14 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.26e-12 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.39e-12 | 35.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.33e-11 | 26.8 |
| MAPK signaling pathway | 信号主干通路 | 1.35e-11 | 20.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.77e-11 | 14.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.67e-11 | 17.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.32e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.94e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.57e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.96e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 5.51e-10 | 20.6 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.23e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.10e-10 | 11.6 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.20e-09 | 43.6 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.20e-09 | 18.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.45_K25</b> — 显著40条, 主干11, 交叉15, 直接OX0, EF=53.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.89e-20 | 56.0 |
| Pathways in cancer | 交叉调控 | 6.47e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 8.61e-16 | 42.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.33e-15 | 28.4 |
| positive regulation of intracellular signal transductio | 其他 | 1.33e-13 | 21.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.33e-13 | 43.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 6.17e-13 | 50.9 |
| Rap1 signaling pathway | 信号主干通路 | 2.02e-12 | 39.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.62e-12 | 25.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.65e-11 | 17.2 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.28e-11 | 20.4 |
| MAPK signaling pathway | 信号主干通路 | 3.52e-11 | 27.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.81e-11 | 61.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.81e-11 | 61.9 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 5.81e-11 | 102.5 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 9.57e-11 | 57.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.15e-10 | 14.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.17e-10 | 47.2 |
| peptidyl-tyrosine phosphorylation (GO:0018108) | 信号主干通路 | 1.69e-10 | 83.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.69e-10 | 83.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.45_K30</b> — 显著40条, 主干15, 交叉14, 直接OX0, EF=43.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.50e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.53e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.93e-16 | 37.8 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.23e-15 | 25.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.46e-14 | 39.6 |
| positive regulation of intracellular signal transductio | 其他 | 8.04e-14 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 8.04e-14 | 46.7 |
| Rap1 signaling pathway | 信号主干通路 | 3.38e-13 | 35.8 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.14e-12 | 23.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.54e-12 | 58.0 |
| MAPK signaling pathway | 信号主干通路 | 8.23e-12 | 25.6 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.34e-11 | 79.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 5.49e-11 | 24.8 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.48e-10 | 12.5 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 1.58e-10 | 56.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.90e-10 | 14.3 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 1.90e-10 | 85.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.94e-10 | 17.0 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.33e-10 | 51.6 |
| positive regulation of cell migration (GO:0030335) | 其他 | 3.34e-10 | 25.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.45_K35</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=34.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.86e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.47e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 4.65e-15 | 32.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.48e-14 | 21.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-13 | 34.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 7.79e-13 | 40.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.54e-12 | 16.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.13e-12 | 30.7 |
| MAPK signaling pathway | 信号主干通路 | 2.13e-12 | 23.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.55e-11 | 19.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 3.32e-11 | 49.7 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 4.89e-11 | 46.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.89e-11 | 67.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.73e-11 | 15.8 |
| positive regulation of cell migration (GO:0030335) | 其他 | 9.73e-11 | 24.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-10 | 13.2 |
| regulation of cell migration (GO:0030334) | 其他 | 2.54e-10 | 17.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.54e-10 | 21.3 |
| regulation of kinase activity (GO:0043549) | 其他 | 3.61e-10 | 49.8 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 4.36e-10 | 48.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.45_K40</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=29.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.64e-21 | 20.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.59e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.63e-14 | 20.3 |
| Rap1 signaling pathway | 信号主干通路 | 3.22e-13 | 29.3 |
| MAPK signaling pathway | 信号主干通路 | 4.14e-13 | 22.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.23e-12 | 29.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.59e-12 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.34e-12 | 35.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.74e-11 | 14.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.31e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.87e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.92e-10 | 40.5 |
| positive regulation of cell migration (GO:0030335) | 其他 | 5.45e-10 | 21.0 |
| positive regulation of cell population proliferation (G | 其他 | 5.99e-10 | 14.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.72e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.42e-10 | 11.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 8.42e-10 | 10.1 |
| positive regulation of protein modification process (GO | 其他 | 1.03e-09 | 23.9 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.11e-09 | 43.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.4_K25</b> — 显著40条, 主干12, 交叉16, 直接OX0, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.51e-18 | 52.0 |
| Pathways in cancer | 交叉调控 | 6.51e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 6.72e-14 | 38.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.86e-13 | 26.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.40e-13 | 18.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.11e-13 | 22.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.62e-12 | 19.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.00e-11 | 45.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.36e-10 | 47.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.80e-10 | 23.2 |
| Bladder cancer | 交叉调控 | 2.84e-10 | 120.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-10 | 35.7 |
| MAPK signaling pathway | 信号主干通路 | 1.28e-09 | 25.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.79e-09 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 2.95e-09 | 30.5 |
| Estrogen signaling pathway | 交叉调控 | 4.76e-09 | 41.9 |
| Chemical carcinogenesis | 交叉调控 | 5.47e-09 | 27.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.49e-09 | 54.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.49e-09 | 54.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.19e-09 | 87.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.4_K30</b> — 显著40条, 主干9, 交叉17, 直接OX0, EF=47.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.47e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.47e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.99e-14 | 35.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.05e-13 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.21e-12 | 36.3 |
| positive regulation of intracellular signal transductio | 其他 | 3.40e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.81e-12 | 42.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.05e-11 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-11 | 18.4 |
| MAPK signaling pathway | 信号主干通路 | 1.36e-11 | 25.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.18e-11 | 21.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.35e-10 | 58.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 2.28e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.72e-10 | 51.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.72e-10 | 51.6 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 6.05e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.4_K35</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=40.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.84e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.46e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.47e-13 | 30.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.64e-12 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.72e-11 | 31.1 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-11 | 36.4 |
| positive regulation of intracellular signal transductio | 其他 | 4.09e-11 | 15.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 6.59e-11 | 46.2 |
| MAPK signaling pathway | 信号主干通路 | 8.10e-11 | 21.9 |
| Rap1 signaling pathway | 信号主干通路 | 8.10e-11 | 27.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-10 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.47e-10 | 13.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.41e-10 | 18.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 5.70e-10 | 49.8 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.97e-10 | 73.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.14e-09 | 44.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.14e-09 | 44.2 |
| Bladder cancer | 交叉调控 | 1.89e-09 | 85.7 |
| Estrogen signaling pathway | 交叉调控 | 1.89e-09 | 34.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.89e-09 | 33.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.4_K40</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=32.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.08e-19 | 19.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.70e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-14 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.26e-12 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.39e-12 | 35.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.33e-11 | 26.8 |
| MAPK signaling pathway | 信号主干通路 | 1.35e-11 | 20.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.77e-11 | 14.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.67e-11 | 17.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.32e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.94e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.57e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.96e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 5.51e-10 | 20.6 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.23e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.10e-10 | 11.6 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.20e-09 | 43.6 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.20e-09 | 18.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.5_K25</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=57.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.97e-21 | 27.8 |
| Proteoglycans in cancer | 交叉调控 | 2.75e-20 | 56.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 8.12e-16 | 42.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.14e-15 | 28.4 |
| positive regulation of intracellular signal transductio | 其他 | 1.40e-13 | 21.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 5.82e-13 | 18.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.82e-13 | 50.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.02e-13 | 22.1 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 5.30e-12 | 39.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 8.23e-12 | 25.5 |
| MAPK signaling pathway | 信号主干通路 | 4.11e-11 | 27.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.48e-11 | 61.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.48e-11 | 61.9 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 5.48e-11 | 102.5 |
| Rap1 signaling pathway | 信号主干通路 | 6.97e-11 | 35.1 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 9.03e-11 | 57.5 |
| Estrogen signaling pathway | 交叉调控 | 9.38e-11 | 47.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.38e-11 | 47.2 |
| Chemical carcinogenesis | 交叉调控 | 1.40e-10 | 30.9 |
| Bladder cancer | 交叉调控 | 1.46e-10 | 120.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.5_K30</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=46.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.67e-21 | 24.5 |
| Proteoglycans in cancer | 交叉调控 | 9.30e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.79e-16 | 37.8 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.12e-15 | 25.4 |
| positive regulation of intracellular signal transductio | 其他 | 9.55e-14 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 9.55e-14 | 46.7 |
| MAPK signaling pathway | 信号主干通路 | 3.04e-13 | 27.9 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.26e-12 | 36.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.09e-12 | 23.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.22e-12 | 58.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 8.85e-12 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.85e-12 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 8.85e-12 | 18.4 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 8.92e-12 | 79.2 |
| Rap1 signaling pathway | 信号主干通路 | 1.08e-11 | 32.5 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.05e-10 | 58.2 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 1.25e-10 | 56.4 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 1.67e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.22e-10 | 51.6 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 3.82e-10 | 24.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.5_K35</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=33.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 8.01e-21 | 22.1 |
| Proteoglycans in cancer | 交叉调控 | 1.45e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 7.52e-17 | 34.7 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.06e-15 | 23.2 |
| MAPK signaling pathway | 信号主干通路 | 7.10e-14 | 25.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.86e-13 | 21.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 9.93e-13 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.48e-12 | 30.7 |
| positive regulation of intracellular signal transductio | 其他 | 1.83e-12 | 16.1 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 9.72e-12 | 31.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 3.17e-11 | 49.7 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 4.67e-11 | 46.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.67e-11 | 67.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.30e-11 | 15.8 |
| positive regulation of cell migration (GO:0030335) | 其他 | 9.30e-11 | 24.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.00e-10 | 13.2 |
| positive regulation of protein modification process (GO | 其他 | 2.42e-10 | 27.4 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.42e-10 | 21.3 |
| regulation of kinase activity (GO:0043549) | 其他 | 3.45e-10 | 49.8 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 4.16e-10 | 48.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.5_K40</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=30.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.64e-21 | 20.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.57e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.63e-14 | 20.3 |
| Rap1 signaling pathway | 信号主干通路 | 3.22e-13 | 29.3 |
| MAPK signaling pathway | 信号主干通路 | 4.14e-13 | 22.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.23e-12 | 29.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.59e-12 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.33e-12 | 35.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.74e-11 | 14.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.31e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.87e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.91e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.52e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.90e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 4.98e-10 | 20.6 |
| positive regulation of cell population proliferation (G | 其他 | 4.98e-10 | 14.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 5.67e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.28e-10 | 11.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 7.28e-10 | 10.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.7_K25</b> — 显著40条, 主干15, 交叉14, 直接OX0, EF=41.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 7.70e-20 | 48.6 |
| Pathways in cancer | 交叉调控 | 1.70e-19 | 26.3 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.83e-19 | 32.5 |
| Proteoglycans in cancer | 交叉调控 | 3.33e-18 | 52.0 |
| MAPK signaling pathway | 信号主干通路 | 1.81e-14 | 33.5 |
| positive regulation of intracellular signal transductio | 其他 | 1.44e-13 | 21.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 7.48e-13 | 50.9 |
| Rap1 signaling pathway | 信号主干通路 | 1.42e-12 | 39.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.97e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 7.63e-12 | 39.6 |
| positive regulation of protein modification process (GO | 其他 | 8.89e-12 | 38.3 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 4.59e-11 | 24.3 |
| regulation of protein phosphorylation (GO:0001932) | 信号主干通路 | 5.91e-11 | 30.8 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.25e-11 | 61.9 |
| Estrogen signaling pathway | 交叉调控 | 9.44e-11 | 47.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.44e-11 | 47.2 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.02e-10 | 57.5 |
| Chemical carcinogenesis | 交叉调控 | 1.40e-10 | 30.9 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 2.09e-10 | 83.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.55e-10 | 15.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.7_K30</b> — 显著40条, 主干13, 交叉16, 直接OX0, EF=35.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.36e-21 | 24.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.79e-20 | 30.4 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.77e-18 | 40.5 |
| Proteoglycans in cancer | 交叉调控 | 7.47e-17 | 43.3 |
| MAPK signaling pathway | 信号主干通路 | 2.87e-13 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.98e-12 | 23.2 |
| positive regulation of intracellular signal transductio | 其他 | 4.23e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 6.79e-12 | 42.4 |
| Rap1 signaling pathway | 信号主干通路 | 1.02e-11 | 32.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.30e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.30e-11 | 18.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 4.90e-11 | 33.0 |
| positive regulation of protein modification process (GO | 其他 | 5.98e-11 | 31.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 3.38e-10 | 51.6 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 3.71e-10 | 20.3 |
| regulation of protein phosphorylation (GO:0001932) | 信号主干通路 | 3.81e-10 | 25.7 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 4.60e-10 | 48.0 |
| Estrogen signaling pathway | 交叉调控 | 5.06e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.06e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 5.62e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.7_K35</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=77.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.46e-24 | 30.4 |
| Pathways in cancer | 交叉调控 | 7.44e-21 | 22.1 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.48e-19 | 37.0 |
| Proteoglycans in cancer | 交叉调控 | 1.35e-17 | 40.0 |
| MAPK signaling pathway | 信号主干通路 | 6.60e-14 | 25.9 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 3.48e-13 | 227.8 |
| amyloid-beta formation (GO:0034205) | 其他 | 4.80e-13 | 390.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.45e-13 | 21.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.77e-13 | 40.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.18e-12 | 16.1 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-12 | 30.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.73e-12 | 31.1 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 8.53e-12 | 251.0 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.05e-11 | 49.7 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 2.79e-11 | 67.9 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 2.79e-11 | 46.2 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 2.79e-11 | 206.7 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 2.79e-11 | 418.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.44e-11 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.13e-11 | 13.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.09_S0.7_K40</b> — 显著40条, 主干10, 交叉15, 直接OX0, EF=74.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.88e-22 | 26.6 |
| Pathways in cancer | 交叉调控 | 2.85e-19 | 19.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 8.44e-18 | 32.4 |
| Proteoglycans in cancer | 交叉调控 | 1.32e-16 | 35.0 |
| MAPK signaling pathway | 信号主干通路 | 5.13e-13 | 22.7 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 1.01e-12 | 199.3 |
| amyloid-beta formation (GO:0034205) | 其他 | 1.19e-12 | 341.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 3.25e-12 | 35.0 |
| membrane protein ectodomain proteolysis (GO:0006509) | 其他 | 3.87e-12 | 156.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.18e-12 | 18.8 |
| Rap1 signaling pathway | 信号主干通路 | 7.40e-12 | 26.8 |
| positive regulation of intracellular signal transductio | 其他 | 1.17e-11 | 14.1 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 2.12e-11 | 219.6 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.92e-11 | 27.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.13e-11 | 43.5 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 6.13e-11 | 180.9 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 6.13e-11 | 366.1 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 9.44e-11 | 59.4 |
| membrane protein proteolysis (GO:0033619) | 其他 | 9.57e-11 | 92.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 9.57e-11 | 40.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.35_K25</b> — 显著40条, 主干15, 交叉15, 直接OX0, EF=35.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.88e-14 | 21.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.47e-13 | 26.4 |
| Proteoglycans in cancer | 交叉调控 | 1.95e-12 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.95e-12 | 39.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.20e-12 | 35.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.37e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.17e-11 | 39.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.60e-11 | 45.8 |
| MAPK signaling pathway | 信号主干通路 | 1.47e-09 | 25.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.87e-09 | 16.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 5.26e-09 | 12.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.77e-09 | 54.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.77e-09 | 54.2 |
| Ras signaling pathway | 信号主干通路 | 6.49e-09 | 28.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.49e-09 | 41.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 7.41e-09 | 24.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.69e-09 | 14.5 |
| positive regulation of protein modification process (GO | 其他 | 1.28e-08 | 30.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.28e-08 | 17.0 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.89e-08 | 71.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.35_K30</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=31.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.90e-13 | 23.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.32e-13 | 32.4 |
| Pathways in cancer | 交叉调控 | 1.16e-12 | 18.0 |
| Proteoglycans in cancer | 交叉调控 | 1.77e-11 | 33.3 |
| Rap1 signaling pathway | 信号主干通路 | 1.77e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.49e-11 | 21.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.10e-10 | 33.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.08e-10 | 38.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.94e-09 | 22.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.27e-09 | 13.2 |
| MAPK signaling pathway | 信号主干通路 | 1.00e-08 | 20.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.51e-08 | 45.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.51e-08 | 45.1 |
| positive regulation of intracellular signal transductio | 其他 | 2.51e-08 | 13.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 3.16e-08 | 34.4 |
| Ras signaling pathway | 信号主干通路 | 3.16e-08 | 23.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 4.74e-08 | 10.7 |
| positive regulation of protein modification process (GO | 其他 | 7.28e-08 | 25.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 7.28e-08 | 59.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-07 | 14.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.35_K35</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=29.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.87e-12 | 20.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.04e-12 | 27.8 |
| Pathways in cancer | 交叉调控 | 1.70e-11 | 15.4 |
| Proteoglycans in cancer | 交叉调控 | 1.07e-10 | 28.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.07e-10 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.48e-10 | 18.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.70e-10 | 28.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.05e-09 | 13.9 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.62e-09 | 32.7 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.00e-09 | 44.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.06e-09 | 59.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 4.06e-09 | 21.4 |
| regulation of cell migration (GO:0030334) | 其他 | 8.04e-09 | 15.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 8.72e-09 | 19.3 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 2.48e-08 | 10.0 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 2.92e-08 | 42.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.92e-08 | 11.3 |
| MAPK signaling pathway | 信号主干通路 | 4.81e-08 | 17.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 4.94e-08 | 38.7 |
| positive regulation of cell migration (GO:0030335) | 其他 | 5.28e-08 | 19.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.35_K40</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=27.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.73e-11 | 17.8 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.32e-11 | 24.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 8.35e-11 | 27.2 |
| Pathways in cancer | 交叉调控 | 1.56e-10 | 13.5 |
| Proteoglycans in cancer | 交叉调控 | 4.91e-10 | 25.0 |
| Rap1 signaling pathway | 信号主干通路 | 4.91e-10 | 24.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 1.34e-09 | 20.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.91e-09 | 15.9 |
| positive regulation of intracellular signal transductio | 其他 | 5.26e-09 | 12.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.26e-09 | 28.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.72e-09 | 38.7 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.14e-08 | 52.0 |
| Gap junction | 交叉调控 | 1.44e-08 | 40.8 |
| Calcium signaling pathway | 信号主干通路 | 2.55e-08 | 19.2 |
| regulation of cell migration (GO:0030334) | 其他 | 3.92e-08 | 13.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 3.92e-08 | 16.9 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 9.49e-08 | 37.0 |
| MAPK signaling pathway | 信号主干通路 | 1.29e-07 | 15.7 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.50e-07 | 8.7 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.50e-07 | 33.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.3_K25</b> — 显著40条, 主干15, 交叉15, 直接OX0, EF=35.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.88e-14 | 21.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.47e-13 | 26.4 |
| Proteoglycans in cancer | 交叉调控 | 1.95e-12 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.95e-12 | 39.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.20e-12 | 35.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.37e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.17e-11 | 39.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.60e-11 | 45.8 |
| MAPK signaling pathway | 信号主干通路 | 1.47e-09 | 25.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.87e-09 | 16.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 5.26e-09 | 12.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.77e-09 | 54.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.77e-09 | 54.2 |
| Ras signaling pathway | 信号主干通路 | 6.49e-09 | 28.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.49e-09 | 41.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 7.41e-09 | 24.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.69e-09 | 14.5 |
| positive regulation of protein modification process (GO | 其他 | 1.28e-08 | 30.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.28e-08 | 17.0 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.89e-08 | 71.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.3_K30</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=31.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.90e-13 | 23.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.32e-13 | 32.4 |
| Pathways in cancer | 交叉调控 | 1.16e-12 | 18.0 |
| Proteoglycans in cancer | 交叉调控 | 1.77e-11 | 33.3 |
| Rap1 signaling pathway | 信号主干通路 | 1.77e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.49e-11 | 21.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.10e-10 | 33.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.08e-10 | 38.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.94e-09 | 22.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.27e-09 | 13.2 |
| MAPK signaling pathway | 信号主干通路 | 1.00e-08 | 20.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.51e-08 | 45.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.51e-08 | 45.1 |
| positive regulation of intracellular signal transductio | 其他 | 2.51e-08 | 13.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 3.16e-08 | 34.4 |
| Ras signaling pathway | 信号主干通路 | 3.16e-08 | 23.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 4.74e-08 | 10.7 |
| positive regulation of protein modification process (GO | 其他 | 7.28e-08 | 25.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 7.28e-08 | 59.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-07 | 14.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.3_K35</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=29.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.87e-12 | 20.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.04e-12 | 27.8 |
| Pathways in cancer | 交叉调控 | 1.70e-11 | 15.4 |
| Proteoglycans in cancer | 交叉调控 | 1.07e-10 | 28.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.07e-10 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.48e-10 | 18.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.70e-10 | 28.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.05e-09 | 13.9 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.62e-09 | 32.7 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.00e-09 | 44.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.06e-09 | 59.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 4.06e-09 | 21.4 |
| regulation of cell migration (GO:0030334) | 其他 | 8.04e-09 | 15.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 8.72e-09 | 19.3 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 2.48e-08 | 10.0 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 2.92e-08 | 42.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.92e-08 | 11.3 |
| MAPK signaling pathway | 信号主干通路 | 4.81e-08 | 17.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 4.94e-08 | 38.7 |
| positive regulation of cell migration (GO:0030335) | 其他 | 5.28e-08 | 19.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.3_K40</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=27.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.73e-11 | 17.8 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.32e-11 | 24.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 8.35e-11 | 27.2 |
| Pathways in cancer | 交叉调控 | 1.56e-10 | 13.5 |
| Proteoglycans in cancer | 交叉调控 | 4.91e-10 | 25.0 |
| Rap1 signaling pathway | 信号主干通路 | 4.91e-10 | 24.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 1.34e-09 | 20.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.91e-09 | 15.9 |
| positive regulation of intracellular signal transductio | 其他 | 5.26e-09 | 12.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.26e-09 | 28.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.72e-09 | 38.7 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.14e-08 | 52.0 |
| Gap junction | 交叉调控 | 1.44e-08 | 40.8 |
| Calcium signaling pathway | 信号主干通路 | 2.55e-08 | 19.2 |
| regulation of cell migration (GO:0030334) | 其他 | 3.92e-08 | 13.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 3.92e-08 | 16.9 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 9.49e-08 | 37.0 |
| MAPK signaling pathway | 信号主干通路 | 1.29e-07 | 15.7 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.50e-07 | 8.7 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.50e-07 | 33.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.45_K25</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=38.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.99e-15 | 28.4 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.17e-14 | 38.9 |
| Pathways in cancer | 交叉调控 | 3.88e-14 | 21.6 |
| Proteoglycans in cancer | 交叉调控 | 1.95e-12 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.95e-12 | 39.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.37e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.18e-11 | 39.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.66e-11 | 45.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 5.67e-10 | 15.8 |
| MAPK signaling pathway | 信号主干通路 | 1.47e-09 | 25.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.26e-09 | 16.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 4.55e-09 | 12.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.18e-09 | 54.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.18e-09 | 54.2 |
| Ras signaling pathway | 信号主干通路 | 6.49e-09 | 28.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.49e-09 | 41.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 6.73e-09 | 24.4 |
| positive regulation of protein modification process (GO | 其他 | 1.30e-08 | 30.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.30e-08 | 17.0 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.91e-08 | 71.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.45_K30</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=29.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.90e-13 | 23.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.33e-13 | 32.4 |
| Pathways in cancer | 交叉调控 | 1.16e-12 | 18.0 |
| Proteoglycans in cancer | 交叉调控 | 1.76e-11 | 33.3 |
| Rap1 signaling pathway | 信号主干通路 | 1.76e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.45e-11 | 21.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.11e-10 | 33.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-10 | 38.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.94e-09 | 22.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.28e-09 | 13.2 |
| MAPK signaling pathway | 信号主干通路 | 9.97e-09 | 20.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.51e-08 | 45.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.51e-08 | 45.1 |
| positive regulation of intracellular signal transductio | 其他 | 2.51e-08 | 13.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 3.14e-08 | 34.4 |
| Ras signaling pathway | 信号主干通路 | 3.14e-08 | 23.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 4.74e-08 | 10.7 |
| positive regulation of protein modification process (GO | 其他 | 7.30e-08 | 25.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 7.30e-08 | 59.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-07 | 14.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.45_K35</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=32.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.82e-14 | 21.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 4.94e-12 | 27.8 |
| Pathways in cancer | 交叉调控 | 1.68e-11 | 15.4 |
| Proteoglycans in cancer | 交叉调控 | 1.06e-10 | 28.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.06e-10 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.42e-10 | 18.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.56e-10 | 28.3 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.98e-09 | 32.7 |
| positive regulation of intracellular signal transductio | 其他 | 1.42e-08 | 12.9 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.42e-08 | 19.3 |
| MAPK signaling pathway | 信号主干通路 | 4.76e-08 | 17.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 5.31e-08 | 11.3 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 7.53e-08 | 38.7 |
| positive regulation of protein serine/threonine kinase  | 其他 | 7.53e-08 | 38.7 |
| negative regulation of programmed cell death (GO:004306 | 其他 | 7.96e-08 | 15.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 8.30e-08 | 19.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.06e-07 | 29.5 |
| Ras signaling pathway | 信号主干通路 | 1.24e-07 | 20.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 2.00e-07 | 50.9 |
| positive regulation of protein modification process (GO | 其他 | 2.52e-07 | 21.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.45_K40</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=29.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.25e-13 | 19.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.27e-11 | 24.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 8.21e-11 | 27.2 |
| Pathways in cancer | 交叉调控 | 1.55e-10 | 13.5 |
| Proteoglycans in cancer | 交叉调控 | 4.88e-10 | 25.0 |
| Rap1 signaling pathway | 信号主干通路 | 4.88e-10 | 24.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 1.32e-09 | 20.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.89e-09 | 15.9 |
| positive regulation of intracellular signal transductio | 其他 | 5.17e-09 | 12.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.17e-09 | 28.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.63e-09 | 38.7 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.12e-08 | 52.0 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 4.29e-08 | 16.9 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 1.03e-07 | 37.0 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.60e-07 | 8.7 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.60e-07 | 33.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.74e-07 | 9.9 |
| MAPK signaling pathway | 信号主干通路 | 1.80e-07 | 15.7 |
| regulation of protein phosphorylation (GO:0001932) | 信号主干通路 | 1.92e-07 | 17.3 |
| negative regulation of programmed cell death (GO:004306 | 其他 | 2.36e-07 | 13.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K25</b> — 显著40条, 主干15, 交叉15, 直接OX0, EF=35.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.88e-14 | 21.6 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.47e-13 | 26.4 |
| Proteoglycans in cancer | 交叉调控 | 1.95e-12 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.95e-12 | 39.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.20e-12 | 35.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.37e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.17e-11 | 39.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.60e-11 | 45.8 |
| MAPK signaling pathway | 信号主干通路 | 1.47e-09 | 25.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.87e-09 | 16.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 5.26e-09 | 12.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.77e-09 | 54.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.77e-09 | 54.2 |
| Ras signaling pathway | 信号主干通路 | 6.49e-09 | 28.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.49e-09 | 41.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 7.41e-09 | 24.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.69e-09 | 14.5 |
| positive regulation of protein modification process (GO | 其他 | 1.28e-08 | 30.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.28e-08 | 17.0 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.89e-08 | 71.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K30</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=31.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.90e-13 | 23.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.32e-13 | 32.4 |
| Pathways in cancer | 交叉调控 | 1.16e-12 | 18.0 |
| Proteoglycans in cancer | 交叉调控 | 1.77e-11 | 33.3 |
| Rap1 signaling pathway | 信号主干通路 | 1.77e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.49e-11 | 21.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.10e-10 | 33.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.08e-10 | 38.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.94e-09 | 22.6 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.27e-09 | 13.2 |
| MAPK signaling pathway | 信号主干通路 | 1.00e-08 | 20.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.51e-08 | 45.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.51e-08 | 45.1 |
| positive regulation of intracellular signal transductio | 其他 | 2.51e-08 | 13.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 3.16e-08 | 34.4 |
| Ras signaling pathway | 信号主干通路 | 3.16e-08 | 23.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 4.74e-08 | 10.7 |
| positive regulation of protein modification process (GO | 其他 | 7.28e-08 | 25.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 7.28e-08 | 59.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-07 | 14.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K35</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=29.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.87e-12 | 20.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.04e-12 | 27.8 |
| Pathways in cancer | 交叉调控 | 1.70e-11 | 15.4 |
| Proteoglycans in cancer | 交叉调控 | 1.07e-10 | 28.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.07e-10 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.48e-10 | 18.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.70e-10 | 28.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.05e-09 | 13.9 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.62e-09 | 32.7 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.00e-09 | 44.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.06e-09 | 59.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 4.06e-09 | 21.4 |
| regulation of cell migration (GO:0030334) | 其他 | 8.04e-09 | 15.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 8.72e-09 | 19.3 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 2.48e-08 | 10.0 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 2.92e-08 | 42.3 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.92e-08 | 11.3 |
| MAPK signaling pathway | 信号主干通路 | 4.81e-08 | 17.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 4.94e-08 | 38.7 |
| positive regulation of cell migration (GO:0030335) | 其他 | 5.28e-08 | 19.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.4_K40</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=27.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.73e-11 | 17.8 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.32e-11 | 24.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 8.35e-11 | 27.2 |
| Pathways in cancer | 交叉调控 | 1.56e-10 | 13.5 |
| Proteoglycans in cancer | 交叉调控 | 4.91e-10 | 25.0 |
| Rap1 signaling pathway | 信号主干通路 | 4.91e-10 | 24.4 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 1.34e-09 | 20.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.91e-09 | 15.9 |
| positive regulation of intracellular signal transductio | 其他 | 5.26e-09 | 12.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.26e-09 | 28.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.72e-09 | 38.7 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.14e-08 | 52.0 |
| Gap junction | 交叉调控 | 1.44e-08 | 40.8 |
| Calcium signaling pathway | 信号主干通路 | 2.55e-08 | 19.2 |
| regulation of cell migration (GO:0030334) | 其他 | 3.92e-08 | 13.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 3.92e-08 | 16.9 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 9.49e-08 | 37.0 |
| MAPK signaling pathway | 信号主干通路 | 1.29e-07 | 15.7 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.50e-07 | 8.7 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.50e-07 | 33.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K25</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=38.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.98e-15 | 28.4 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.17e-14 | 38.9 |
| Pathways in cancer | 交叉调控 | 3.91e-14 | 21.6 |
| Proteoglycans in cancer | 交叉调控 | 1.96e-12 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.96e-12 | 39.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.41e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.18e-11 | 39.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.65e-11 | 45.8 |
| MAPK signaling pathway | 信号主干通路 | 1.48e-09 | 25.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.90e-09 | 16.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 5.30e-09 | 12.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.82e-09 | 54.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.82e-09 | 54.2 |
| Ras signaling pathway | 信号主干通路 | 6.53e-09 | 28.3 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 6.53e-09 | 41.3 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 7.47e-09 | 24.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.75e-09 | 14.5 |
| positive regulation of protein modification process (GO | 其他 | 1.29e-08 | 30.7 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.29e-08 | 17.0 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-08 | 71.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K30</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=29.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.04e-15 | 25.4 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.34e-13 | 32.4 |
| Pathways in cancer | 交叉调控 | 1.16e-12 | 18.0 |
| Proteoglycans in cancer | 交叉调控 | 1.77e-11 | 33.3 |
| Rap1 signaling pathway | 信号主干通路 | 1.77e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.49e-11 | 21.2 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.11e-10 | 33.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.10e-10 | 38.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.92e-10 | 51.6 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.40e-09 | 69.3 |
| positive regulation of intracellular signal transductio | 其他 | 1.40e-09 | 15.0 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.85e-09 | 22.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 2.85e-09 | 11.6 |
| MAPK signaling pathway | 信号主干通路 | 1.00e-08 | 20.9 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 1.10e-08 | 49.3 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.88e-08 | 45.1 |
| regulation of cell migration (GO:0030334) | 其他 | 2.28e-08 | 16.7 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 3.16e-08 | 34.4 |
| Ras signaling pathway | 信号主干通路 | 3.16e-08 | 23.6 |
| positive regulation of protein modification process (GO | 其他 | 6.32e-08 | 25.5 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K35</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=60.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.28e-17 | 24.6 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 4.93e-12 | 27.8 |
| Pathways in cancer | 交叉调控 | 1.70e-11 | 15.4 |
| Proteoglycans in cancer | 交叉调控 | 1.07e-10 | 28.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.07e-10 | 27.9 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 1.32e-10 | 195.2 |
| amyloid-beta formation (GO:0034205) | 其他 | 3.93e-10 | 325.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.93e-10 | 28.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.48e-10 | 18.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.32e-09 | 32.7 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.68e-09 | 44.2 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 3.49e-09 | 209.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 3.49e-09 | 59.4 |
| positive regulation of intracellular signal transductio | 其他 | 7.32e-09 | 12.9 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 7.32e-09 | 19.3 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 7.32e-09 | 172.3 |
| ephrin receptor signaling pathway (GO:0048013) | 其他 | 8.37e-09 | 50.0 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 1.33e-08 | 10.3 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.78e-08 | 10.0 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 2.03e-08 | 334.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.5_K40</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=53.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 7.98e-16 | 21.6 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.27e-11 | 24.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 8.20e-11 | 27.2 |
| Pathways in cancer | 交叉调控 | 1.55e-10 | 13.5 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 2.44e-10 | 170.8 |
| Proteoglycans in cancer | 交叉调控 | 4.88e-10 | 25.0 |
| Rap1 signaling pathway | 信号主干通路 | 4.88e-10 | 24.4 |
| amyloid-beta formation (GO:0034205) | 其他 | 7.65e-10 | 284.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.89e-09 | 15.9 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.16e-09 | 28.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.62e-09 | 38.7 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 7.55e-09 | 183.0 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 9.97e-09 | 52.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 1.44e-08 | 18.7 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 1.69e-08 | 150.7 |
| ephrin receptor signaling pathway (GO:0048013) | 其他 | 2.59e-08 | 43.7 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.97e-08 | 16.9 |
| positive regulation of intracellular signal transductio | 其他 | 3.74e-08 | 11.3 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 3.94e-08 | 292.9 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 6.03e-08 | 37.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K25</b> — 显著40条, 主干13, 交叉14, 直接OX0, EF=103.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 9.16e-23 | 36.5 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.93e-14 | 38.9 |
| Pathways in cancer | 交叉调控 | 3.71e-14 | 21.6 |
| amyloid-beta formation (GO:0034205) | 其他 | 5.62e-14 | 546.7 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 1.50e-12 | 351.4 |
| Proteoglycans in cancer | 交叉调控 | 2.19e-12 | 40.0 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 4.93e-12 | 289.4 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 6.02e-12 | 273.3 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 6.02e-12 | 585.7 |
| ephrin receptor signaling pathway (GO:0048013) | 其他 | 7.25e-12 | 80.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 8.12e-12 | 25.5 |
| membrane protein ectodomain proteolysis (GO:0006509) | 其他 | 2.09e-11 | 213.9 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 2.09e-11 | 45.8 |
| Notch receptor processing (GO:0007220) | 信号主干通路 | 2.30e-11 | 455.6 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 6.59e-11 | 14.4 |
| Rap1 signaling pathway | 信号主干通路 | 8.59e-11 | 35.1 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.57e-10 | 35.7 |
| membrane protein proteolysis (GO:0033619) | 其他 | 4.56e-10 | 126.2 |
| positive regulation of intracellular signal transductio | 其他 | 1.20e-09 | 16.5 |
| MAPK signaling pathway | 信号主干通路 | 1.40e-09 | 25.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K30</b> — 显著40条, 主干14, 交叉14, 直接OX0, EF=85.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-20 | 30.4 |
| amyloid-beta formation (GO:0034205) | 其他 | 3.10e-13 | 455.6 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 3.38e-13 | 32.4 |
| Pathways in cancer | 交叉调控 | 1.14e-12 | 18.0 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 5.52e-12 | 292.9 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 6.19e-12 | 13.8 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 1.51e-11 | 241.2 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 1.70e-11 | 488.1 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 1.70e-11 | 227.8 |
| Proteoglycans in cancer | 交叉调控 | 2.04e-11 | 33.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.39e-11 | 33.0 |
| ephrin receptor signaling pathway (GO:0048013) | 其他 | 3.39e-11 | 66.7 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 3.43e-11 | 24.8 |
| membrane protein ectodomain proteolysis (GO:0006509) | 其他 | 5.72e-11 | 178.3 |
| Notch receptor processing (GO:0007220) | 信号主干通路 | 5.72e-11 | 379.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.75e-11 | 21.2 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 1.11e-10 | 38.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.87e-10 | 51.6 |
| negative regulation of programmed cell death (GO:004306 | 其他 | 2.84e-10 | 19.7 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 3.04e-10 | 24.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K35</b> — 显著40条, 主干15, 交叉14, 直接OX0, EF=73.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 8.23e-19 | 26.1 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 4.80e-13 | 227.8 |
| amyloid-beta formation (GO:0034205) | 其他 | 5.88e-13 | 390.5 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.41e-12 | 27.8 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 1.25e-11 | 251.0 |
| Pathways in cancer | 交叉调控 | 1.65e-11 | 15.4 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 4.21e-11 | 206.7 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 4.21e-11 | 418.4 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 7.09e-11 | 11.8 |
| Proteoglycans in cancer | 交叉调控 | 1.23e-10 | 28.6 |
| ephrin receptor signaling pathway (GO:0048013) | 其他 | 1.55e-10 | 57.1 |
| Notch receptor processing (GO:0007220) | 信号主干通路 | 1.74e-10 | 325.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.74e-10 | 28.3 |
| membrane protein ectodomain proteolysis (GO:0006509) | 其他 | 1.74e-10 | 152.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.16e-10 | 21.3 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.52e-10 | 32.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 7.09e-10 | 18.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 7.66e-10 | 44.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.25e-09 | 12.3 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.72e-09 | 59.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0.1_S0.7_K40</b> — 显著40条, 主干15, 交叉13, 直接OX0, EF=67.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.93e-17 | 22.8 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 1.36e-12 | 199.3 |
| amyloid-beta formation (GO:0034205) | 其他 | 1.43e-12 | 341.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.56e-11 | 24.3 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 3.04e-11 | 219.6 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 8.91e-11 | 366.1 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 8.91e-11 | 180.9 |
| Pathways in cancer | 交叉调控 | 1.52e-10 | 13.5 |
| Notch receptor processing (GO:0007220) | 信号主干通路 | 4.58e-10 | 284.7 |
| ephrin receptor signaling pathway (GO:0048013) | 其他 | 5.05e-10 | 50.0 |
| membrane protein ectodomain proteolysis (GO:0006509) | 其他 | 5.05e-10 | 133.7 |
| regulation of apoptotic process (GO:0042981) | 信号主干通路 | 5.52e-10 | 10.4 |
| Proteoglycans in cancer | 交叉调控 | 5.63e-10 | 25.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 7.25e-10 | 24.8 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.15e-09 | 18.6 |
| Calcium signaling pathway | 信号主干通路 | 1.79e-09 | 21.4 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 2.12e-09 | 28.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.51e-09 | 38.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.82e-09 | 15.9 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 5.37e-09 | 52.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.35_K25</b> — 显著40条, 主干12, 交叉16, 直接OX0, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.51e-18 | 52.0 |
| Pathways in cancer | 交叉调控 | 6.51e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 6.72e-14 | 38.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.86e-13 | 26.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.40e-13 | 18.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.11e-13 | 22.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.62e-12 | 19.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.00e-11 | 45.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.36e-10 | 47.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.80e-10 | 23.2 |
| Bladder cancer | 交叉调控 | 2.84e-10 | 120.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-10 | 35.7 |
| MAPK signaling pathway | 信号主干通路 | 1.28e-09 | 25.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.79e-09 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 2.95e-09 | 30.5 |
| Estrogen signaling pathway | 交叉调控 | 4.76e-09 | 41.9 |
| Chemical carcinogenesis | 交叉调控 | 5.47e-09 | 27.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.49e-09 | 54.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.49e-09 | 54.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.19e-09 | 87.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.35_K30</b> — 显著40条, 主干9, 交叉17, 直接OX0, EF=47.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.47e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.47e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.99e-14 | 35.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.05e-13 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.21e-12 | 36.3 |
| positive regulation of intracellular signal transductio | 其他 | 3.40e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.81e-12 | 42.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.05e-11 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-11 | 18.4 |
| MAPK signaling pathway | 信号主干通路 | 1.36e-11 | 25.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.18e-11 | 21.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.35e-10 | 58.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 2.28e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.72e-10 | 51.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.72e-10 | 51.6 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 6.05e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.35_K35</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=40.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.84e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.46e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.47e-13 | 30.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.64e-12 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.72e-11 | 31.1 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-11 | 36.4 |
| positive regulation of intracellular signal transductio | 其他 | 4.09e-11 | 15.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 6.59e-11 | 46.2 |
| MAPK signaling pathway | 信号主干通路 | 8.10e-11 | 21.9 |
| Rap1 signaling pathway | 信号主干通路 | 8.10e-11 | 27.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-10 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.47e-10 | 13.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.41e-10 | 18.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 5.70e-10 | 49.8 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.97e-10 | 73.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.14e-09 | 44.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.14e-09 | 44.2 |
| Bladder cancer | 交叉调控 | 1.89e-09 | 85.7 |
| Estrogen signaling pathway | 交叉调控 | 1.89e-09 | 34.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.89e-09 | 33.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.35_K40</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=32.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.08e-19 | 19.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.70e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-14 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.26e-12 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.39e-12 | 35.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.33e-11 | 26.8 |
| MAPK signaling pathway | 信号主干通路 | 1.35e-11 | 20.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.77e-11 | 14.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.67e-11 | 17.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.32e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.94e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.57e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.96e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 5.51e-10 | 20.6 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.23e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.10e-10 | 11.6 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.20e-09 | 43.6 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.20e-09 | 18.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.3_K25</b> — 显著40条, 主干12, 交叉16, 直接OX0, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.51e-18 | 52.0 |
| Pathways in cancer | 交叉调控 | 6.51e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 6.72e-14 | 38.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.86e-13 | 26.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.40e-13 | 18.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.11e-13 | 22.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.62e-12 | 19.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.00e-11 | 45.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.36e-10 | 47.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.80e-10 | 23.2 |
| Bladder cancer | 交叉调控 | 2.84e-10 | 120.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-10 | 35.7 |
| MAPK signaling pathway | 信号主干通路 | 1.28e-09 | 25.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.79e-09 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 2.95e-09 | 30.5 |
| Estrogen signaling pathway | 交叉调控 | 4.76e-09 | 41.9 |
| Chemical carcinogenesis | 交叉调控 | 5.47e-09 | 27.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.49e-09 | 54.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.49e-09 | 54.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.19e-09 | 87.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.3_K30</b> — 显著40条, 主干9, 交叉17, 直接OX0, EF=47.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.47e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.47e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.99e-14 | 35.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.05e-13 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.21e-12 | 36.3 |
| positive regulation of intracellular signal transductio | 其他 | 3.40e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.81e-12 | 42.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.05e-11 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-11 | 18.4 |
| MAPK signaling pathway | 信号主干通路 | 1.36e-11 | 25.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.18e-11 | 21.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.35e-10 | 58.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 2.28e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.72e-10 | 51.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.72e-10 | 51.6 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 6.05e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.3_K35</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=40.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.84e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.46e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.47e-13 | 30.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.64e-12 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.72e-11 | 31.1 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-11 | 36.4 |
| positive regulation of intracellular signal transductio | 其他 | 4.09e-11 | 15.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 6.59e-11 | 46.2 |
| MAPK signaling pathway | 信号主干通路 | 8.10e-11 | 21.9 |
| Rap1 signaling pathway | 信号主干通路 | 8.10e-11 | 27.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-10 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.47e-10 | 13.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.41e-10 | 18.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 5.70e-10 | 49.8 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.97e-10 | 73.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.14e-09 | 44.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.14e-09 | 44.2 |
| Bladder cancer | 交叉调控 | 1.89e-09 | 85.7 |
| Estrogen signaling pathway | 交叉调控 | 1.89e-09 | 34.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.89e-09 | 33.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.3_K40</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=32.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.08e-19 | 19.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.70e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-14 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.26e-12 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.39e-12 | 35.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.33e-11 | 26.8 |
| MAPK signaling pathway | 信号主干通路 | 1.35e-11 | 20.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.77e-11 | 14.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.67e-11 | 17.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.32e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.94e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.57e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.96e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 5.51e-10 | 20.6 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.23e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.10e-10 | 11.6 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.20e-09 | 43.6 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.20e-09 | 18.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.45_K25</b> — 显著40条, 主干11, 交叉15, 直接OX0, EF=53.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 5.89e-20 | 56.0 |
| Pathways in cancer | 交叉调控 | 6.47e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 8.61e-16 | 42.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.33e-15 | 28.4 |
| positive regulation of intracellular signal transductio | 其他 | 1.33e-13 | 21.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.33e-13 | 43.6 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 6.17e-13 | 50.9 |
| Rap1 signaling pathway | 信号主干通路 | 2.02e-12 | 39.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.62e-12 | 25.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.65e-11 | 17.2 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 2.28e-11 | 20.4 |
| MAPK signaling pathway | 信号主干通路 | 3.52e-11 | 27.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.81e-11 | 61.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.81e-11 | 61.9 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 5.81e-11 | 102.5 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 9.57e-11 | 57.5 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.15e-10 | 14.0 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.17e-10 | 47.2 |
| peptidyl-tyrosine phosphorylation (GO:0018108) | 信号主干通路 | 1.69e-10 | 83.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.69e-10 | 83.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.45_K30</b> — 显著40条, 主干15, 交叉14, 直接OX0, EF=43.8</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.50e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.53e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.93e-16 | 37.8 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.23e-15 | 25.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.46e-14 | 39.6 |
| positive regulation of intracellular signal transductio | 其他 | 8.04e-14 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 8.04e-14 | 46.7 |
| Rap1 signaling pathway | 信号主干通路 | 3.38e-13 | 35.8 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.14e-12 | 23.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.54e-12 | 58.0 |
| MAPK signaling pathway | 信号主干通路 | 8.23e-12 | 25.6 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.34e-11 | 79.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 5.49e-11 | 24.8 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 1.48e-10 | 12.5 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 1.58e-10 | 56.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.90e-10 | 14.3 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 1.90e-10 | 85.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.94e-10 | 17.0 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.33e-10 | 51.6 |
| positive regulation of cell migration (GO:0030335) | 其他 | 3.34e-10 | 25.4 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.45_K35</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=34.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.86e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.47e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 4.65e-15 | 32.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 4.48e-14 | 21.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-13 | 34.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 7.79e-13 | 40.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.54e-12 | 16.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.13e-12 | 30.7 |
| MAPK signaling pathway | 信号主干通路 | 2.13e-12 | 23.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.55e-11 | 19.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 3.32e-11 | 49.7 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 4.89e-11 | 46.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.89e-11 | 67.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.73e-11 | 15.8 |
| positive regulation of cell migration (GO:0030335) | 其他 | 9.73e-11 | 24.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-10 | 13.2 |
| regulation of cell migration (GO:0030334) | 其他 | 2.54e-10 | 17.2 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.54e-10 | 21.3 |
| regulation of kinase activity (GO:0043549) | 其他 | 3.61e-10 | 49.8 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 4.36e-10 | 48.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.45_K40</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=29.4</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.64e-21 | 20.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.59e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.63e-14 | 20.3 |
| Rap1 signaling pathway | 信号主干通路 | 3.22e-13 | 29.3 |
| MAPK signaling pathway | 信号主干通路 | 4.14e-13 | 22.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.23e-12 | 29.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.59e-12 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.34e-12 | 35.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.74e-11 | 14.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.31e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.87e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.92e-10 | 40.5 |
| positive regulation of cell migration (GO:0030335) | 其他 | 5.45e-10 | 21.0 |
| positive regulation of cell population proliferation (G | 其他 | 5.99e-10 | 14.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.72e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.42e-10 | 11.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 8.42e-10 | 10.1 |
| positive regulation of protein modification process (GO | 其他 | 1.03e-09 | 23.9 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.11e-09 | 43.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.4_K25</b> — 显著40条, 主干12, 交叉16, 直接OX0, EF=41.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Proteoglycans in cancer | 交叉调控 | 6.51e-18 | 52.0 |
| Pathways in cancer | 交叉调控 | 6.51e-18 | 24.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 6.72e-14 | 38.9 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.86e-13 | 26.4 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.40e-13 | 18.5 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.11e-13 | 22.1 |
| positive regulation of intracellular signal transductio | 其他 | 3.62e-12 | 19.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.00e-11 | 45.8 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 2.36e-10 | 47.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.80e-10 | 23.2 |
| Bladder cancer | 交叉调控 | 2.84e-10 | 120.0 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 3.35e-10 | 35.7 |
| MAPK signaling pathway | 信号主干通路 | 1.28e-09 | 25.1 |
| Rap1 signaling pathway | 信号主干通路 | 2.79e-09 | 31.2 |
| Lipid and atherosclerosis | 交叉调控 | 2.95e-09 | 30.5 |
| Estrogen signaling pathway | 交叉调控 | 4.76e-09 | 41.9 |
| Chemical carcinogenesis | 交叉调控 | 5.47e-09 | 27.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.49e-09 | 54.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.49e-09 | 54.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.19e-09 | 87.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.4_K30</b> — 显著40条, 主干9, 交叉17, 直接OX0, EF=47.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.47e-19 | 23.2 |
| Proteoglycans in cancer | 交叉调控 | 9.47e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.99e-14 | 35.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.05e-13 | 23.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.21e-12 | 36.3 |
| positive regulation of intracellular signal transductio | 其他 | 3.40e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.81e-12 | 42.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.05e-11 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.05e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.05e-11 | 18.4 |
| MAPK signaling pathway | 信号主干通路 | 1.36e-11 | 25.6 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-11 | 32.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 6.18e-11 | 21.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.35e-10 | 58.2 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 2.28e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.72e-10 | 51.6 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.72e-10 | 51.6 |
| Estrogen signaling pathway | 交叉调控 | 5.45e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.45e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 6.05e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.4_K35</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=40.9</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.84e-19 | 21.0 |
| Proteoglycans in cancer | 交叉调控 | 1.46e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.47e-13 | 30.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.64e-12 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.72e-11 | 31.1 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.09e-11 | 36.4 |
| positive regulation of intracellular signal transductio | 其他 | 4.09e-11 | 15.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 6.59e-11 | 46.2 |
| MAPK signaling pathway | 信号主干通路 | 8.10e-11 | 21.9 |
| Rap1 signaling pathway | 信号主干通路 | 8.10e-11 | 27.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.41e-10 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.47e-10 | 13.2 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.41e-10 | 18.2 |
| regulation of kinase activity (GO:0043549) | 其他 | 5.70e-10 | 49.8 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 7.97e-10 | 73.2 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 1.14e-09 | 44.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.14e-09 | 44.2 |
| Bladder cancer | 交叉调控 | 1.89e-09 | 85.7 |
| Estrogen signaling pathway | 交叉调控 | 1.89e-09 | 34.2 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 1.89e-09 | 33.7 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.4_K40</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=32.2</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 3.08e-19 | 19.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.70e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.65e-14 | 20.3 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.26e-12 | 29.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.39e-12 | 35.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.33e-11 | 26.8 |
| MAPK signaling pathway | 信号主干通路 | 1.35e-11 | 20.9 |
| positive regulation of intracellular signal transductio | 其他 | 1.77e-11 | 14.1 |
| PI3K-Akt signaling pathway | 信号主干通路 | 9.67e-11 | 17.4 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.32e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.90e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.94e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.57e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.96e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 5.51e-10 | 20.6 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.23e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.10e-10 | 11.6 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.20e-09 | 43.6 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 1.20e-09 | 18.6 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.5_K25</b> — 显著40条, 主干10, 交叉17, 直接OX0, EF=57.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 1.97e-21 | 27.8 |
| Proteoglycans in cancer | 交叉调控 | 2.75e-20 | 56.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 8.12e-16 | 42.1 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.14e-15 | 28.4 |
| positive regulation of intracellular signal transductio | 其他 | 1.40e-13 | 21.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 5.82e-13 | 18.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.82e-13 | 50.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.02e-13 | 22.1 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 5.30e-12 | 39.6 |
| PI3K-Akt signaling pathway | 信号主干通路 | 8.23e-12 | 25.5 |
| MAPK signaling pathway | 信号主干通路 | 4.11e-11 | 27.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 5.48e-11 | 61.9 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 5.48e-11 | 61.9 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 5.48e-11 | 102.5 |
| Rap1 signaling pathway | 信号主干通路 | 6.97e-11 | 35.1 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 9.03e-11 | 57.5 |
| Estrogen signaling pathway | 交叉调控 | 9.38e-11 | 47.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.38e-11 | 47.2 |
| Chemical carcinogenesis | 交叉调控 | 1.40e-10 | 30.9 |
| Bladder cancer | 交叉调控 | 1.46e-10 | 120.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.5_K30</b> — 显著40条, 主干14, 交叉15, 直接OX0, EF=46.5</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.67e-21 | 24.5 |
| Proteoglycans in cancer | 交叉调控 | 9.30e-19 | 46.7 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 2.79e-16 | 37.8 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 2.12e-15 | 25.4 |
| positive regulation of intracellular signal transductio | 其他 | 9.55e-14 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 9.55e-14 | 46.7 |
| MAPK signaling pathway | 信号主干通路 | 3.04e-13 | 27.9 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 1.26e-12 | 36.3 |
| PI3K-Akt signaling pathway | 信号主干通路 | 2.09e-12 | 23.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.22e-12 | 58.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 8.85e-12 | 53.9 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 8.85e-12 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 8.85e-12 | 18.4 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 8.92e-12 | 79.2 |
| Rap1 signaling pathway | 信号主干通路 | 1.08e-11 | 32.5 |
| regulation of kinase activity (GO:0043549) | 其他 | 1.05e-10 | 58.2 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 1.25e-10 | 56.4 |
| peptidyl-tyrosine modification (GO:0018212) | 其他 | 1.67e-10 | 85.4 |
| regulation of phosphatidylinositol 3-kinase signaling ( | 其他 | 2.22e-10 | 51.6 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 3.82e-10 | 24.9 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.5_K35</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=33.3</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 8.01e-21 | 22.1 |
| Proteoglycans in cancer | 交叉调控 | 1.45e-17 | 40.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 7.52e-17 | 34.7 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.06e-15 | 23.2 |
| MAPK signaling pathway | 信号主干通路 | 7.10e-14 | 25.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.86e-13 | 21.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 9.93e-13 | 40.0 |
| Rap1 signaling pathway | 信号主干通路 | 1.48e-12 | 30.7 |
| positive regulation of intracellular signal transductio | 其他 | 1.83e-12 | 16.1 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 9.72e-12 | 31.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 3.17e-11 | 49.7 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 4.67e-11 | 46.2 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 4.67e-11 | 67.9 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 9.30e-11 | 15.8 |
| positive regulation of cell migration (GO:0030335) | 其他 | 9.30e-11 | 24.0 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.00e-10 | 13.2 |
| positive regulation of protein modification process (GO | 其他 | 2.42e-10 | 27.4 |
| MAPK cascade (GO:0000165) | 信号主干通路 | 2.42e-10 | 21.3 |
| regulation of kinase activity (GO:0043549) | 其他 | 3.45e-10 | 49.8 |
| regulation of MAP kinase activity (GO:0043405) | 信号主干通路 | 4.16e-10 | 48.3 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.5_K40</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=30.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 7.64e-21 | 20.3 |
| Proteoglycans in cancer | 交叉调控 | 1.42e-16 | 35.0 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 9.57e-16 | 30.4 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.63e-14 | 20.3 |
| Rap1 signaling pathway | 信号主干通路 | 3.22e-13 | 29.3 |
| MAPK signaling pathway | 信号主干通路 | 4.14e-13 | 22.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.23e-12 | 29.7 |
| PI3K-Akt signaling pathway | 信号主干通路 | 3.59e-12 | 18.8 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 4.33e-12 | 35.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.74e-11 | 14.1 |
| positive regulation of protein serine/threonine kinase  | 其他 | 1.31e-10 | 43.5 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 1.87e-10 | 59.4 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.91e-10 | 40.5 |
| positive regulation of phospholipase C activity (GO:001 | 其他 | 3.52e-10 | 83.4 |
| positive regulation of cell migration (GO:0030335) | 其他 | 4.90e-10 | 21.0 |
| positive regulation of MAPK cascade (GO:0043410) | 信号主干通路 | 4.98e-10 | 20.6 |
| positive regulation of cell population proliferation (G | 其他 | 4.98e-10 | 14.1 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 5.67e-10 | 13.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.28e-10 | 11.6 |
| regulation of cell population proliferation (GO:0042127 | 其他 | 7.28e-10 | 10.1 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.7_K25</b> — 显著40条, 主干15, 交叉14, 直接OX0, EF=41.0</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 7.70e-20 | 48.6 |
| Pathways in cancer | 交叉调控 | 1.70e-19 | 26.3 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 5.83e-19 | 32.5 |
| Proteoglycans in cancer | 交叉调控 | 3.33e-18 | 52.0 |
| MAPK signaling pathway | 信号主干通路 | 1.81e-14 | 33.5 |
| positive regulation of intracellular signal transductio | 其他 | 1.44e-13 | 21.0 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 7.48e-13 | 50.9 |
| Rap1 signaling pathway | 信号主干通路 | 1.42e-12 | 39.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.97e-12 | 25.5 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 7.63e-12 | 39.6 |
| positive regulation of protein modification process (GO | 其他 | 8.89e-12 | 38.3 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 4.59e-11 | 24.3 |
| regulation of protein phosphorylation (GO:0001932) | 信号主干通路 | 5.91e-11 | 30.8 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.25e-11 | 61.9 |
| Estrogen signaling pathway | 交叉调控 | 9.44e-11 | 47.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 9.44e-11 | 47.2 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 1.02e-10 | 57.5 |
| Chemical carcinogenesis | 交叉调控 | 1.40e-10 | 30.9 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 2.09e-10 | 83.2 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 2.55e-10 | 15.8 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.7_K30</b> — 显著40条, 主干13, 交叉16, 直接OX0, EF=35.1</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| Pathways in cancer | 交叉调控 | 5.36e-21 | 24.5 |
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.79e-20 | 30.4 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 1.77e-18 | 40.5 |
| Proteoglycans in cancer | 交叉调控 | 7.47e-17 | 43.3 |
| MAPK signaling pathway | 信号主干通路 | 2.87e-13 | 27.9 |
| PI3K-Akt signaling pathway | 信号主干通路 | 1.98e-12 | 23.2 |
| positive regulation of intracellular signal transductio | 其他 | 4.23e-12 | 17.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 6.79e-12 | 42.4 |
| Rap1 signaling pathway | 信号主干通路 | 1.02e-11 | 32.5 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 1.30e-11 | 15.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 1.30e-11 | 18.4 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 4.90e-11 | 33.0 |
| positive regulation of protein modification process (GO | 其他 | 5.98e-11 | 31.9 |
| positive regulation of protein serine/threonine kinase  | 其他 | 3.38e-10 | 51.6 |
| positive regulation of protein phosphorylation (GO:0001 | 信号主干通路 | 3.71e-10 | 20.3 |
| regulation of protein phosphorylation (GO:0001932) | 信号主干通路 | 3.81e-10 | 25.7 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 4.60e-10 | 48.0 |
| Estrogen signaling pathway | 交叉调控 | 5.06e-10 | 39.9 |
| Fluid shear stress and atherosclerosis | 交叉调控 | 5.06e-10 | 39.3 |
| Bladder cancer | 交叉调控 | 5.62e-10 | 100.0 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.7_K35</b> — 显著40条, 主干12, 交叉15, 直接OX0, EF=77.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 3.46e-24 | 30.4 |
| Pathways in cancer | 交叉调控 | 7.44e-21 | 22.1 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 5.48e-19 | 37.0 |
| Proteoglycans in cancer | 交叉调控 | 1.35e-17 | 40.0 |
| MAPK signaling pathway | 信号主干通路 | 6.60e-14 | 25.9 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 3.48e-13 | 227.8 |
| amyloid-beta formation (GO:0034205) | 其他 | 4.80e-13 | 390.5 |
| PI3K-Akt signaling pathway | 信号主干通路 | 5.45e-13 | 21.5 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 5.77e-13 | 40.0 |
| positive regulation of intracellular signal transductio | 其他 | 1.18e-12 | 16.1 |
| Rap1 signaling pathway | 信号主干通路 | 1.37e-12 | 30.7 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 6.73e-12 | 31.1 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 8.53e-12 | 251.0 |
| positive regulation of protein serine/threonine kinase  | 其他 | 2.05e-11 | 49.7 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 2.79e-11 | 67.9 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 2.79e-11 | 46.2 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 2.79e-11 | 206.7 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 2.79e-11 | 418.4 |
| cellular response to cytokine stimulus (GO:0071345) | 信号主干通路 | 6.44e-11 | 15.8 |
| cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 7.13e-11 | 13.2 |
| ... 还有 20 条通路 ... | | | |

</details>

<details>
<summary><b>P0_S0.7_K40</b> — 显著40条, 主干10, 交叉15, 直接OX0, EF=74.6</summary>

| 通路名 | 类别 | P-value | EF |
|--------|------|--------|----|
| transmembrane receptor protein tyrosine kinase signalin | 信号主干通路 | 1.88e-22 | 26.6 |
| Pathways in cancer | 交叉调控 | 2.85e-19 | 19.3 |
| positive regulation of phosphorylation (GO:0042327) | 信号主干通路 | 8.44e-18 | 32.4 |
| Proteoglycans in cancer | 交叉调控 | 1.32e-16 | 35.0 |
| MAPK signaling pathway | 信号主干通路 | 5.13e-13 | 22.7 |
| amyloid precursor protein metabolic process (GO:0042982 | 其他 | 1.01e-12 | 199.3 |
| amyloid-beta formation (GO:0034205) | 其他 | 1.19e-12 | 341.7 |
| positive regulation of protein kinase B signaling (GO:0 | 其他 | 3.25e-12 | 35.0 |
| membrane protein ectodomain proteolysis (GO:0006509) | 其他 | 3.87e-12 | 156.0 |
| PI3K-Akt signaling pathway | 信号主干通路 | 4.18e-12 | 18.8 |
| Rap1 signaling pathway | 信号主干通路 | 7.40e-12 | 26.8 |
| positive regulation of intracellular signal transductio | 其他 | 1.17e-11 | 14.1 |
| amyloid precursor protein catabolic process (GO:0042987 | 其他 | 2.12e-11 | 219.6 |
| regulation of protein kinase B signaling (GO:0051896) | 其他 | 2.92e-11 | 27.2 |
| positive regulation of protein serine/threonine kinase  | 其他 | 6.13e-11 | 43.5 |
| amyloid-beta metabolic process (GO:0050435) | 其他 | 6.13e-11 | 180.9 |
| Notch receptor processing, ligand-dependent (GO:0035333 | 信号主干通路 | 6.13e-11 | 366.1 |
| positive regulation of MAP kinase activity (GO:0043406) | 信号主干通路 | 9.44e-11 | 59.4 |
| membrane protein proteolysis (GO:0033619) | 其他 | 9.57e-11 | 92.0 |
| positive regulation of kinase activity (GO:0033674) | 其他 | 9.57e-11 | 40.5 |
| ... 还有 20 条通路 ... | | | |

</details>

## 5. 靶点稳健性分析

| 象限 | 说明 | 靶点数 |
|------|------|--------|
| Q1: 核心稳健 | 高频率 + 高度 | 33 |
| Q2: 低频高重要 | 低频率 + 高度 | 1 |
| Q3: 稳健但边缘 | 高频率 + 低度 | 1 |
| Q4: 可忽略 | 低频率 + 低度 | 32 |

### 5.1 Q1 核心稳健靶点列表

| 排名 | 基因 | 频率 | 出现次数 | 平均Degree(原始) | 平均Degree(标准化) | 平均排名 | 平均中介中心性 |
|------|------|------|---------|----------------|-------------------|---------|--------------|
| 1 | KDR | 100.0% | 72 | 25.7 | 0.599 | 10.9 | 0.0263 |
| 2 | BCL2 | 100.0% | 72 | 34.3 | 0.802 | 4.5 | 0.0513 |
| 3 | PTGS2 | 100.0% | 72 | 29.6 | 0.683 | 8.1 | 0.0444 |
| 4 | EGFR | 100.0% | 72 | 36.9 | 0.883 | 2.7 | 0.0742 |
| 5 | SRC | 100.0% | 72 | 40.0 | 0.934 | 2.2 | 0.2010 |
| 6 | HSP90AA1 | 100.0% | 72 | 25.7 | 0.603 | 10.8 | 0.0221 |
| 7 | STAT3 | 100.0% | 72 | 33.1 | 0.784 | 5.2 | 0.0532 |
| 8 | KIT | 100.0% | 72 | 18.4 | 0.435 | 16.7 | 0.0223 |
| 9 | FGF2 | 100.0% | 72 | 24.2 | 0.585 | 12.4 | 0.0100 |
| 10 | MET | 100.0% | 72 | 18.6 | 0.449 | 15.7 | 0.0031 |
| 11 | MMP9 | 100.0% | 72 | 32.9 | 0.777 | 5.6 | 0.0454 |
| 12 | TNF | 100.0% | 72 | 41.3 | 0.956 | 1.6 | 0.1165 |
| 13 | MMP2 | 100.0% | 72 | 23.2 | 0.552 | 12.4 | 0.0080 |
| 14 | IL2 | 97.2% | 70 | 18.0 | 0.414 | 17.7 | 0.0056 |
| 15 | LGALS3 | 94.4% | 68 | 16.9 | 0.398 | 20.7 | 0.0218 |
| 16 | FGFR1 | 86.1% | 62 | 15.2 | 0.364 | 22.8 | 0.0090 |
| 17 | BACE1 | 72.2% | 52 | 12.1 | 0.319 | 26.4 | 0.0142 |
| 18 | FGF1 | 69.4% | 50 | 11.8 | 0.328 | 24.4 | 0.0010 |
| 19 | MAPK14 | 66.7% | 48 | 20.0 | 0.407 | 18.8 | 0.0106 |
| 20 | ERBB2 | 66.7% | 48 | 31.0 | 0.640 | 9.6 | 0.0237 |
| 21 | ESR1 | 66.7% | 48 | 35.0 | 0.703 | 8.6 | 0.0519 |
| 22 | STAT1 | 66.7% | 48 | 19.7 | 0.407 | 18.9 | 0.0040 |
| 23 | DRD2 | 66.7% | 48 | 16.4 | 0.365 | 24.3 | 0.0455 |
| 24 | HIF1A | 66.7% | 48 | 34.0 | 0.695 | 8.2 | 0.0246 |
| 25 | APP | 66.7% | 48 | 27.4 | 0.559 | 13.3 | 0.0536 |
| 26 | FLT1 | 66.7% | 48 | 19.2 | 0.388 | 19.8 | 0.0035 |
| 27 | ABCB1 | 59.7% | 43 | 16.9 | 0.355 | 23.6 | 0.0081 |
| 28 | CA9 | 59.7% | 43 | 14.3 | 0.342 | 24.6 | 0.0324 |
| 29 | MMP1 | 58.3% | 42 | 17.6 | 0.343 | 24.4 | 0.0016 |
| 30 | ESR2 | 55.6% | 40 | 16.4 | 0.341 | 25.6 | 0.0122 |
| 31 | PGF | 52.8% | 38 | 9.7 | 0.322 | 22.8 | 0.0013 |
| 32 | MAOB | 52.8% | 38 | 19.9 | 0.366 | 23.6 | 0.0358 |
| 33 | MAPT | 50.0% | 36 | 16.5 | 0.349 | 26.3 | 0.0207 |

## 6. 通路排名稳定性 (Top 25)

| 排名 | 通路 | 类别 | 频率 | 平均排名 | 排名波动 | 平均logP | 平均EF |
|------|------|------|------|---------|---------|---------|-------|
| 1 | Pathways in cancer | 交叉调控 | 100.0% | 1.9 | 1.4 | 16.3 | 20.7 |
| 2 | transmembrane receptor protein tyrosine kinase sig | 信号主干通路 | 100.0% | 2.8 | 1.4 | 14.3 | 24.2 |
| 2 | positive regulation of phosphorylation (GO:0042327 | 信号主干通路 | 100.0% | 2.8 | 0.7 | 13.8 | 33.8 |
| 3 | Proteoglycans in cancer | 交叉调控 | 100.0% | 3.3 | 2.1 | 14.9 | 39.6 |
| 5 | amyloid-beta formation (GO:0034205) | 其他 | 13.9% | 5.5 | 2.2 | 11.7 | 380.9 |
| 5 | amyloid precursor protein metabolic process (GO:00 | 其他 | 13.9% | 5.6 | 2.2 | 11.5 | 214.8 |
| 7 | regulation of protein kinase B signaling (GO:00518 | 其他 | 100.0% | 7.6 | 3.3 | 10.8 | 33.2 |
| 8 | positive regulation of protein kinase B signaling  | 其他 | 100.0% | 8.0 | 2.3 | 10.7 | 39.7 |
| 8 | positive regulation of intracellular signal transd | 其他 | 97.2% | 8.9 | 3.3 | 10.5 | 16.1 |
| 9 | Rap1 signaling pathway | 信号主干通路 | 100.0% | 9.2 | 4.6 | 10.6 | 30.9 |
| 9 | PI3K-Akt signaling pathway | 信号主干通路 | 100.0% | 9.4 | 2.9 | 10.4 | 20.8 |
| 9 | amyloid precursor protein catabolic process (GO:00 | 其他 | 13.9% | 9.6 | 4.0 | 10.5 | 244.8 |
| 11 | MAPK signaling pathway | 信号主干通路 | 100.0% | 11.3 | 5.5 | 10.0 | 23.2 |
| 12 | amyloid-beta metabolic process (GO:0050435) | 其他 | 13.9% | 12.3 | 4.7 | 10.0 | 201.6 |
| 13 | regulation of apoptotic process (GO:0042981) | 信号主干通路 | 8.3% | 13.7 | 5.9 | 9.3 | 11.6 |
| 13 | positive regulation of protein serine/threonine ki | 其他 | 100.0% | 13.8 | 3.2 | 9.3 | 49.0 |
| 13 | Notch receptor processing, ligand-dependent (GO:00 | 信号主干通路 | 13.9% | 13.8 | 5.3 | 9.9 | 405.5 |
| 14 | positive regulation of kinase activity (GO:0033674 | 其他 | 66.7% | 14.0 | 4.7 | 9.9 | 48.1 |
| 14 | cytokine-mediated signaling pathway (GO:0019221) | 信号主干通路 | 94.4% | 14.1 | 5.2 | 9.5 | 13.9 |
| 15 | cellular response to cytokine stimulus (GO:0071345 | 信号主干通路 | 83.3% | 15.1 | 5.9 | 9.6 | 16.8 |
| 15 | membrane protein ectodomain proteolysis (GO:000650 | 其他 | 13.9% | 15.7 | 6.1 | 9.7 | 153.5 |
| 16 | ephrin receptor signaling pathway (GO:0048013) | 其他 | 13.9% | 16.8 | 6.0 | 9.6 | 56.2 |
| 17 | positive regulation of MAP kinase activity (GO:004 | 信号主干通路 | 98.6% | 17.3 | 4.6 | 8.9 | 65.8 |
| 17 | Notch receptor processing (GO:0007220) | 信号主干通路 | 11.1% | 17.4 | 5.2 | 9.9 | 333.2 |
| 17 | MAPK cascade (GO:0000165) | 信号主干通路 | 66.7% | 17.8 | 5.2 | 8.6 | 20.6 |

## 7. Pvalue 阈值敏感性对比

| 指标 | P=0.01 | P=0.05 |
|------|--------|--------|
| 显著通路总数 | 40 | 40 |
| 主干道 | 14 | 14 |
| 交叉调控 | 15 | 15 |
| 直接OX | 0 | 0 |
| 其他 | 11 | 11 |
| 平均EF | 35.41 | 35.41 |
| 共有通路 | 40 | 40 |
| 分类分布一致 | 是（Pvalue 稳健） | 是（Pvalue 稳健） |

## 8. 最终判定

| 检查项 | 结果 | 判定 |
|--------|------|------|
| Q1核心稳健靶点数>=5 | 33 | PASS |
| 关键氧化应激靶点进入Q1 | 8个(TNF,SRC,EGFR,STAT3,BCL2,MMP9,HIF1A,PTGS2) | PASS |
| 排名稳定通路数>=3 | 9 | PASS |
| 信号网络主干道排名稳定 | 4个 | PASS |

**通过: 4/4**

> ✅ **判定: 敏感性分析通过。结果支持故事成立。**