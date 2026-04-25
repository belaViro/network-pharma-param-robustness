# ML 元模型可行性分析：从化合物描述符预测网药最优参数

## 一、问题重述

**目标**: 建立机器学习模型，输入化合物的理化性质，输出该化合物在网药分析中的最优参数配置。

```
输入: MW, HBD, HBA, LogP, TPSA, NumRotatableBonds, ...
        ↓
模型: ML 回归/分类
        ↓
输出: prob_th, string_th, pval_th, top_k (或某个综合评分)
```

**学术价值**: 如果能实现，将解决网药分析中参数选择主观性强的问题，实现"自适应参数配置"。

---

## 二、可行性评估

### 2.1 数据需求

| 需求项 | 最低要求 | 理想状态 | 当前状态 |
|--------|---------|---------|---------|
| 化合物数量 | 30-50 | 100+ | **5** |
| 每个化合物的参数扫描 | 54 组合 | 200 组合 | 已可实现 |
| 描述符维度 | 5-10 个 | 20+ (Morgan指纹等) | 可从 PubChem 获取 |
| 标注数据 (最优参数) | 30-50 条 | 100+ 条 | **0** |

**结论**: 当前 5 个化合物远远不够。需要至少再收集 25-45 个化合物。

### 2.2 标注数据的获取成本

对每个化合物，需要：

1. **获取 SMILES** (PubChem API, ~0.5s)
2. **SwissTarget 预测** (浏览器自动化, ~1min/化合物)
3. **跑 54 种参数组合** 的分析流程 (~10s/组合 = 9min/化合物)
4. **定义"最优参数"** (需要人工或自动化的评价标准)

**50 个化合物的总成本**:
- SwissTarget: 50 × 1min = ~1 小时
- 参数扫描: 50 × 9min = ~7.5 小时
- **总计**: 约 1 天可完成数据收集

### 2.3 核心难点：如何定义"最优参数"？

这是整个想法的最大问题。

**场景 A：有明确的疾病方向**
- 假设研究的是"抗氧化应激"
- "最优"= 氧化应激通路数最多 + CI 最窄
- 问题：这个定义只对"抗氧化"有效，换"抗炎"就完全不同
- 模型只能对特定疾病类型预测

**场景 B：无明确疾病方向**
- "最优"= 显著通路数最多？
- 问题：通路数多不代表有意义，可能只是噪声
- "最优"= Enrichment Factor 最高？
- 问题：EF 高可能是特定通路的偶然发现

**场景 C：多目标优化**
- 定义一个综合评分（如当前的 composite_score）
- 问题：评分权重是主观的，不同研究者可能用不同权重

**根本问题**: 网药分析的"最优参数"取决于**研究目的**（发现什么类型的通路），而非化合物本身。

### 2.4 更深层的问题：理化性质 → 最优参数，逻辑链条是否成立？

假设模型学到了：

> "分子量 MW > 400 的化合物，最优 STRING 阈值是 0.3；MW < 400 的，最优 STRING 阈值是 0.5"

**这个规律的物理解释是什么？**

- MW 大的化合物 → 靶点预测可能更少/更不可靠 → 需要更宽松的 STRING 阈值来保留足够网络边？
- 这个因果关系非常间接，中间经过了很多环节

**对比：为什么这个方向难发好论文**

| 任务 | 输入-输出关系 | 物理基础 | 可行性 |
|------|-------------|---------|--------|
| 从 MW/LogP 预测溶解度 | 强相关 | 热力学 | 高 |
| 从结构预测毒性 | 中等相关 | ADMET | 中等 |
| 从描述符预测网药最优参数 | **弱相关** | **无直接物理基础** | **低** |

---

## 三、更合理的替代方案

### 方案 1：统计分析而非 ML

**不训练预测模型，而是做相关性分析**：

```python
# 问题：化合物描述符与最优参数之间是否存在统计关联？

# 方法：
# 1. 收集 30-50 个化合物
# 2. 对每个化合物跑参数扫描
# 3. 找出每个化合物的"最优参数"
# 4. 做 Pearson/Spearman 相关：描述符 vs 最优参数
# 5. 报告相关系数和显著性

# 例如：
# MW 与 optimal_string_th 的相关性: r = 0.45, p < 0.05
# LogP 与 optimal_top_k 的相关性: r = -0.32, p = 0.08
```

**优点**:
- 不需要大量样本训练 ML 模型
- 结果可直接解释为统计关联，不 claiming 预测能力
- 论文中可以诚实写"我们发现 XX 描述符与最优参数存在正相关"

**缺点**:
- 只能发现线性/单调关系
- 不能给出"新化合物的最优参数"

### 方案 2：规则引擎而非 ML

基于领域知识建立启发式规则：

```python
def suggest_parameters(compound_descriptors):
    mw = compound_descriptors['MW']
    n_targets = compound_descriptors['predicted_target_count']
    
    # 规则 1: 靶点少的化合物 → 降低 STRING 阈值以保留网络连通性
    if n_targets < 50:
        string_th = 0.3
    else:
        string_th = 0.5
    
    # 规则 2: 小分子 → 降低 TopK（靶点少，Hub 也少）
    if mw < 300:
        top_k = 20
    else:
        top_k = 30
    
    return {'string_th': string_th, 'top_k': top_k}
```

**优点**:
- 可解释性强
- 不需要训练数据
- 可直接写入 Methods

**缺点**:
- 规则是专家经验，需要文献支持
- 不够"智能"

### 方案 3：聚焦真正有意义的问题

与其预测"最优参数"，不如回答：

> **"化合物的哪些特征使其网药分析结果对参数选择更敏感/更稳健？"**

这是一个更有学术价值的问题：

- 输入：化合物描述符
- 输出：参数敏感性指标（如不同参数组合下氧化应激通路数的标准差）
- 发现：MW 大的化合物，参数敏感性更高（结果更不稳定）

**物理解释**: MW 大的化合物靶点预测可能更分散/不可靠，导致下游分析对参数更敏感。

**论文定位**: "Identifying Compound Characteristics Associated with Parameter Sensitivity in Network Pharmacology"

---

## 四、如果坚持要做 ML 预测模型

### 实施路径

**Phase 1: 数据收集（1-2 周）**

```python
# 1. 从文献/数据库收集 30-50 个已知抗氧化/抗炎活性化合物
#    来源：DrugBank, ChEMBL, 文献 SI

# 2. 从 PubChem 批量获取描述符
import pubchempy as pcp

def get_descriptors(compound_name):
    compound = pcp.get_compounds(compound_name, 'name')[0]
    return {
        'MW': compound.molecular_weight,
        'HBD': compound.h_bond_donor_count,
        'HBA': compound.h_bond_acceptor_count,
        'LogP': compound.xlogp,
        'TPSA': compound.tpsa,
        'NumRotatableBonds': compound.rotatable_bond_count,
    }

# 3. 跑 SwissTarget + 参数扫描（54 组合 × 30 Bootstrap）
# 4. 对每个化合物定义最优参数（如 reliability_score 最高）
```

**Phase 2: 模型训练（2-3 天）**

```python
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

# 特征: 化合物描述符
X = df[['MW', 'HBD', 'HBA', 'LogP', 'TPSA', 'NumRotatableBonds']]

# 目标: 最优参数（4个目标变量，可分别训练4个模型）
y_prob = df['optimal_prob_th']
y_string = df['optimal_string_th']
y_pval = df['optimal_pval_th']
y_topk = df['optimal_top_k']

# 模型
model = RandomForestRegressor(n_estimators=100, random_state=42)
scores = cross_val_score(model, X, y_string, cv=5, scoring='neg_mean_squared_error')
```

**Phase 3: 验证（1 周）**

- 留 5-10 个化合物做测试集
- 看模型预测的最优参数 vs 真实最优参数的误差
- 如果 MSE 很小 → 模型有效
- 如果 MSE 很大 → 想法不成立

### 预期结果（悲观）

| 指标 | 预期值 |
|------|--------|
| 交叉验证 R² | 0.1-0.3 |
| 测试集 MSE | 较大 |
| 结论 | 化合物描述符与最优参数之间**无显著预测关系** |

**原因**: 网药分析的最优参数更多取决于**研究目的**和**疾病类型**，而非化合物本身的理化性质。

---

## 五、总结建议

### 我的诚实建议

**不要做 ML 预测模型。**

原因：
1. 需要 30-50+ 个化合物，数据收集成本高
2. "最优参数"的定义与研究目的强相关，不是一个化合物固有的属性
3. 理化性质 → 最优参数之间缺乏直接的物理因果链条
4. 即使做出来，预测精度大概率很低，审稿人会质疑

**建议改为**:

> **"统计分析：化合物特征与参数敏感性的关联"**

聚焦问题：
- MW / LogP 大的化合物，网药分析结果是否对参数选择更敏感？
- 哪类化合物需要更严格的过滤阈值？

方法：
- 收集 20-30 个化合物（比 ML 需求少）
- 跑参数扫描
- 计算每个化合物的"参数敏感性指数"（不同参数组合下结果的标准差）
- 做相关性分析

论文定位：
> "我们首次系统评估了化合物理化性质对网药分析参数敏感性的影响，发现 XX 类化合物需要更谨慎的参数选择。"

这比"ML 预测最优参数"更扎实、更可信、更容易发。
