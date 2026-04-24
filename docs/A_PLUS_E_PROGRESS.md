# A+E 项目进度评估报告

> 评估时间：2026-04-20
> 评估依据：代码审查 + `src/outputs/` 实际数据状态

---

## 一、一句话结论

**传统方法流水线已完成（~95%），A+E 代码框架刚完成（~90%），但 A+E 从未跑通出结果。当前总体进度约 40%，最危险的瓶颈是：A+E 代码尚未验证能否在真实数据上跑完。**

---

## 二、三大板块拆解

### 板块 1：传统网络药理学流水线（Stage 1-6 + 分支A）

| 模块 | 完成度 | 证据 | 剩余工作 |
|:---|:---|:---|:---|
| **Stage 1: PubChem 查 SMILES** | ✅ 100% | 5 个分子的 SMILES 已成功获取 | 无需改动 |
| **Stage 2: SwissTarget 预测** | ✅ 100% | `swiss_targets/` 下有 5 个完整 CSV | 无需改动 |
| **Stage 3: 合并去重靶点** | ✅ 100% | `merged_targets.csv` 存在 (3.6KB) | 无需改动 |
| **Stage 4: 疾病基因合并** | ✅ 100% | `merged_disease_genes.csv` 存在 (1.1MB)，含 GeneCards + OMIM + TTD | 无需改动 |
| **Stage 5: 交集分析** | ✅ 100% | `intersected_targets.csv` 存在 (4.6KB) | 无需改动 |
| **Stage 6: PPI 网络构建** | ✅ 100% | STRING 边文件 14.8KB，748 节点，hub_nodes/degrees 齐全 | 无需改动 |
| **传统富集分析 (GO/KEGG)** | ✅ 100% | `traditional/enrichment/` 下 4 个 CSV 齐全 | 无需改动 |
| **传统可视化** | ✅ 100% | PPI 网络图、韦恩图、富集气泡图、柱状图均已生成 | 无需改动 |
| **传统方法评估报告** | ✅ 100% | 多个 `.md` 报告存在 | 无需改动 |

**板块 1 结论：这是你最坚实的资产。5 个草莓叶活性成分（Sucrose, Trehalose, Catechin, Procyanidin B1, Scutellarin）的完整传统分析已经跑通，数据、图表、报告齐全。可直接作为论文 baseline。**

---

### 板块 2：GNN / Hybrid / Bio-GNN（已官方弃用）

| 模块 | 完成度 | 状态说明 |
|:---|:---|:---|
| `run_gnn_pipeline.py` (Hybrid) | ❌ 已弃用 | 有历史输出但不再使用 |
| `run_bio_gnn_pipeline.py` | ❌ 已弃用 | 有历史输出但不再使用 |
| `src/gnn/` 全部代码 | ❌ 已弃用 | 保留代码但不进入论文 |
| GNN vs Degree 对比 | ❌ 已弃用 | 历史报告在 `gnn/evaluation/` 下，不采用 |

**板块 2 结论：这些代码和输出是"沉没成本"。建议论文中完全不提 GNN，避免审稿人追问"为什么做了又不用的方法"。如果非要提，放在 Introduction 里作为"我们前期尝试发现 GNN 在此场景无增益"一句话带过即可。**

---

### 板块 3：A+E 联合稳健性分析（当前核心创新点）

| 模块 | 完成度 | 证据 | 剩余工作 |
|:---|:---|:---|:---|
| **A 部分代码（参数网格扫描）** | ✅ 90% | 原 `_param_sensitivity_analysis.py` 有框架 | 刚重写完成，待验证 |
| **E 部分代码（Bootstrap + CI）** | ✅ 90% | 刚补完的版本含 `N_BOOTSTRAP=100`、95% CI 计算 | **尚未在真实数据上跑通** |
| **A+E 联合输出（热图/表格）** | ❌ 0% | `param_sensitivity_results.csv` 为空文件 (0 字节) | **必须跑通并产出结果** |
| **稳健靶点识别（含 CI）** | ❌ 0% | 无输出文件 | 依赖 A+E 运行结果 |
| **核心通路识别（含 CI）** | ❌ 0% | 无输出文件 | 依赖 A+E 运行结果 |
| **带 CI 的可视化** | ❌ 0% | 无 `param_sensitivity_A_plus_E.png/pdf` | 依赖 A+E 运行结果 |
| **可靠性评分模型** | ✅ 90% | 代码已写 `reliability_score = ox_mean × 0.6 + (1-CI_width) × 0.4` | 需跑完后调参验证 |

**板块 3 结论：这是整个项目最大的不确定性。代码看起来完整，但从未在真实数据上成功执行。`src/outputs/gnn/param_sensitivity_results.csv` 是空文件，说明之前的版本跑失败了，或者根本没跑完。**

---

## 三、数据资产清单

### 已有且可用的数据

```
src/outputs/shared/
├── swiss_targets/           # 5 个分子的 SwissTarget 预测结果 ✅
├── merged_targets.csv       # 去重合并靶点 ✅
├── disease_genes/           # GeneCards + OMIM + TTD 疾病基因 ✅
├── intersected_targets.csv  # 靶点-疾病基因交集 ✅
├── ppi_network/             # STRING PPI 网络（边+节点+度数+hub）✅
└── figures/                 # PPI 网络图 + 韦恩图 ✅

src/outputs/traditional/
├── enrichment/              # GO_BP, GO_CC, GO_MF, KEGG 富集结果 ✅
├── figures/                 # 富集可视化图 ✅
└── evaluation/              # 传统方法评估报告 ✅
```

### 缺失且必须补的数据

```
src/outputs/gnn/
├── param_sensitivity_A_plus_E_results.csv      ❌ 主结果表（含点估计+CI）
├── robust_targets_with_ci.csv                  ❌ 稳健靶点（含 Bootstrap CI）
├── core_pathways_with_ci.csv                   ❌ 核心通路（含 Bootstrap CI）
├── param_sensitivity_A_plus_E.png              ❌ 9 子图综合可视化
└── param_sensitivity_A_plus_E.pdf              ❌ 论文插图矢量版
```

---

## 四、风险点与瓶颈分析

### 🔴 高风险：A+E 代码可能跑不通

**症状**：`param_sensitivity_results.csv` 是空文件。

**可能原因**：
1. 原代码有 bug，在第一个参数组合就崩溃
2. `run_enrichment_analysis()` 在循环中被反复调用，某次因网络/API 问题失败导致整个脚本中断
3. 参数组合过多（25 组 × 100 Bootstrap = 2500 次富集），运行时间太长被手动中断

**验证方法**：
```bash
cd D:\科研\网络药理学ppi筛选改进\code
python _param_sensitivity_analysis.py
```
建议先用 `N_BOOTSTRAP=5` 跑一个小规模测试，确认能输出非空 CSV 后再调回 100。

### 🟡 中风险：Bootstrap 计算成本

**症状**：如果每次 `run_enrichment_analysis()` 调用 gseapy/Enrichr API 需要 2-5 秒，2500 次需要 1.5-3.5 小时。

**缓解方案**：
- 方案 A：降低 Bootstrap 到 50-100（论文常用 100，够用）
- 方案 B：加缓存机制（同参数组合+同基因集的结果缓存到本地 SQLite/CSV）
- 方案 C：并行化（多进程跑不同参数组合的 Bootstrap）

### 🟢 低风险：传统方法结果不够漂亮

**症状**：如果 A+E 最优参数组合的结果和传统固定参数差异不大，故事性会变弱。

**缓解方案**：A+E 的核心卖点不是"结果更好"，而是"结果更可信"。即使氧化应激通路数相同，带 95% CI 的结论比单点估计更有说服力。论文定位是"方法论可靠性增强"而非"发现更多通路"。

---

## 五、通往投稿的最短路径

### Phase 1：验证 A+E 代码（1-2 天）

| 任务 | 预计时间 | 验收标准 |
|:---|:---|:---|
| 小规模测试跑通（N_BOOTSTRAP=5, 精简参数网格）| 2-4 小时 | 输出非空 `param_sensitivity_results.csv` |
| 全参数网格正式运行（N_BOOTSTRAP=100）| 2-6 小时 | 产出完整的 3 个 CSV + PNG + PDF |
| Debug 和修 bug | 0.5-1 天 | 所有参数组合都有 Bootstrap 均值和 CI |

### Phase 2：对比分析（2-3 天）

| 任务 | 预计时间 | 验收标准 |
|:---|:---|:---|
| 固定参数 vs A+E 最优参数的 Top-30 靶点对比 | 0.5 天 | 有重叠率、独有靶点列表 |
| 固定参数 vs A+E 最优参数的氧化应激通路对比 | 0.5 天 | 有通路数、EF、CI 对比表 |
| 数据库验证：稳健靶点 vs OMIM/DisGeNET 交叉验证 | 1-2 天 | 超几何检验 p < 0.05 |
| 消融实验：A 单独 vs E 单独 vs A+E（证明 1+1>2）| 1 天 | 有对比表格 |

### Phase 3：论文写作（2-3 周）

| 任务 | 预计时间 |
|:---|:---|
| Introduction（背景 + GNN 尝试简述 + A+E 提出动机）| 2-3 天 |
| Methods（传统流水线 + A+E 框架数学描述）| 3-5 天 |
| Results（传统结果 + A+E 结果 + 对比 + 验证）| 5-7 天 |
| Discussion（A+E 意义 + 局限 + 未来方向）| 2-3 天 |
| 图表精修 + 补充材料整理 | 2-3 天 |

### Phase 4：投稿（1 周）

| 任务 | 预计时间 |
|:---|:---|
| 选刊 + 格式调整 | 2-3 天 |
| Cover letter + 亮点提炼 | 1-2 天 |
| 代码整理 + GitHub 开源 | 1-2 天 |

---

## 六、总体进度量化

| 维度 | 权重 | 完成度 | 加权得分 |
|:---|:---|:---|:---|
| 传统流水线数据产出 | 20% | 95% | 19.0 |
| 传统方法分析/图表/报告 | 15% | 95% | 14.25 |
| A+E 方法学代码 | 15% | 90% | 13.5 |
| A+E 结果产出 | 15% | 0% | 0 |
| 对比分析与验证 | 15% | 0% | 0 |
| 论文写作 | 15% | 0% | 0 |
| 代码开源与可复现性 | 5% | 0% | 0 |
| **总计** | **100%** | — | **46.75%** |

**修正后总体进度：约 40-45%**

> 注：之前我按"代码写完就算 25%"估计偏低了。实际上传统方法的数据资产非常完整，占了将近 33% 的权重。真正的卡点只有 A+E 的运行和后续分析写作。

---

## 七、本周优先行动清单

1. **【今天】** 用 `N_BOOTSTRAP=5` 跑 `_param_sensitivity_analysis.py`，确认不报错、输出非空
2. **【明天】** 调回 `N_BOOTSTRAP=100` 跑完整版，产出所有 CSV 和图表
3. **【后天】** 做固定参数 vs A+E 最优参数的对比表（靶点 + 通路）
4. **【本周内】** 做 OMIM/DisGeNET 数据库交叉验证（超几何检验）

---

## 八、给导师汇报可用的版本

如果导师问"做到哪了"，用这个版本：

> "传统网络药理学流水线已经完整跑通，5 个草莓叶活性成分的靶点预测、PPI 网络构建、富集分析和可视化都已经完成。前期尝试的 GNN 改进方案经评估后发现无显著增益，已转向更核心的方法论创新——A+E 联合稳健性分析框架（参数敏感性 + Bootstrap 不确定性量化）。A+E 的代码已经实现，正在验证运行中，预计本周内出完整结果。下一步是对比固定参数法和 A+E 优化参数法的差异，并用疾病数据库验证稳健靶点的可靠性，随后进入论文写作阶段。"
