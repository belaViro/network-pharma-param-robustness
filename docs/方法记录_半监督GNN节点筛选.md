# 半监督GNN节点筛选方法

## 1. 背景与目标

### 1.1 研究问题
传统网络药理学使用 Degree Centrality 筛选PPI网络核心靶点，本研究引入图神经网络(GNN)进行节点重要性排序，与传统方法做消融实验对比。

### 1.2 核心假设
- 节点出现在关键通路中的频率反映其生物学重要性
- GNN能通过邻域聚合学习到比单纯degree更丰富的结构信息
- 边的权重(score)蕴含相互作用置信度信息

### 1.3 数据规模
| 数据 | 数量 |
|------|------|
| 输入靶点数 | 112 |
| PPI边数 | ~500-830 (无阈值/0.4阈值) |
| Top Hub数 | 30 |

---

## 2. 数据描述

### 2.1 PPI网络
- 来源：STRING数据库
- 节点：蛋白质/基因
- 边：蛋白质相互作用
- 边权重：combined_score (0-1)

### 2.2 节点特征
| 特征类型 | 来源 | 维度 |
|----------|------|------|
| 度中心性 | 计算 | 1 |
| 介数中心性 | 计算 | 1 |
| 接近中心性 | 计算 | 1 |
| PageRank | 计算 | 1 |
| 蛋白嵌入 | ESM2/ProtTrans (可选) | 128/512 |

### 2.3 伪标签设计
用富集通路信息构建节点重要性分数：

```
对于节点v，其重要性分数：
score(v) = Σ_{pathway in Pathways} (1 / rank_in_pathway) × pathway_weight

其中：
- pathway_weight = -log10(P-value)
- rank_in_pathway = 节点在该通路基因列表中的排名
```

---

## 3. 模型框架

### 3.1 GraphSAGE + Attention

```python
class SemiSupervisedGNN(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        # SAGEConv层
        self.conv1 = SAGEConv(input_dim, hidden_dim)
        self.conv2 = SAGEConv(hidden_dim, hidden_dim)
        # 注意力聚合
        self.attention = AttentionPooling(hidden_dim)
        # 重要性预测头
        self.importance_head = nn.Linear(hidden_dim, 1)

    def forward(self, x, edge_index, edge_weight=None):
        # 边权重用于注意力
        h = F.relu(self.conv1(x, edge_index))
        h = F.relu(self.conv2(h, edge_index))
        # 注意力加权聚合
        h = self.attention(h, edge_index, edge_weight)
        # 预测重要性分数
        importance = self.importance_head(h)
        return importance
```

### 3.2 损失函数

```python
def semi_supervised_loss(pred_importance, true_importance, mask):
    """
    pred_importance: 模型预测的节点重要性
    true_importance: 基于通路富集的伪标签
    mask: 有标签节点 (True/False)
    """
    # 有监督部分
    supervised_loss = F.mse_loss(pred_importance[mask], true_importance[mask])

    # 无监督部分 - 一致性正则
    # 同一连通分量的节点应该有相似的重要性
    unsupervised_loss = graph_regularization(pred_importance, edge_index)

    return supervised_loss + λ * unsupervised_loss
```

---

## 4. 训练策略

### 4.1 半监督训练流程

```
Epoch 1-N:
    1. 前向传播 → 预测重要性分数
    2. 计算半监督损失
    3. 反向传播更新GNN参数
    4. 用更新后的模型预测 → 重新计算伪标签 (EMA)
    5. 重复
```

### 4.2 参数设置

| 参数 | 值 | 说明 |
|------|-----|------|
| hidden_dim | 64 | 隐层维度 |
| learning_rate | 0.001 | 学习率 |
| weight_decay | 5e-4 | L2正则 |
| λ (正则权重) | 0.1 | 无监督正则强度 |
| epochs | 200 | 训练轮数 |
| patience | 30 | 早停轮数 |

---

## 5. 对比实验设计

### 5.1 实验组 vs 对照组

| 组别 | 方法 | 说明 |
|------|------|------|
| **Control** | Degree Centrality | 传统方法，取degree最大的30个节点 |
| **GNN-Baseline** | GraphSAGE (无注意力) | 不用attention，用mean pooling |
| **GNN-Attention** | GraphSAGE + Attention | 用边权重做attention聚合 |
| **GNN-Semi** | 半监督GraphSAGE | 用通路伪标签监督 |

### 5.2 消融实验

| 消融项 | 变体 | 预期结果 |
|--------|------|----------|
| 边权重 | 有/无 attention | 边权重attention应该提升 |
| 伪标签 | 有/无 半监督 | 半监督应该提升 |
| 节点特征 | 仅degree / 完整特征 | 丰富特征应该提升 |

---

## 6. 评估指标

### 6.1 节点排序评估

| 指标 | 计算方法 | 说明 |
|------|----------|------|
| Top-K 重合率 | \|GNN_Top30 ∩ Degree_Top30\| / 30 | 两方法选出一致的节点比例 |
| Kendall's τ | 排序相关性 | 两方法排序的等级相关 |
| Spearman ρ | 排序相关性 | 更关注相对顺序 |

### 6.2 富集分析评估

| 指标 | 计算方法 | 说明 |
|------|----------|------|
| 通路覆盖数 | Top30节点富集到的通路数 | 越多越好 |
| 通路重合率 | 与论文通路重合比例 | 与原方法对比 |
| P-value 平均值 | 富集分析平均显著水平 | 越小越显著 |

### 6.3 可视化

- PPI网络图：内圈GNN-Hub，外圈传统Hub
- 富集气泡图：对比两组通路
- 节点重要性散点图：GNN score vs Degree

---

## 7. 预期结论

### 7.1 可能的结果

**结果A：GNN显著优于Degree**
- GNN选出的节点富集到更多通路
- 与原论文方法更一致
- 说明邻域信息有价值

**结果B：GNN与Degree无显著差异**
- 两者选出的Hub高度重合
- 说明在小规模网络上，degree已经是充分统计量

**结果C：Degree反而更好**
- 可能因为数据量太小GNN过拟合
- 或者边权重信息噪声太大

### 7.2 论文贡献点

1. 提出半监督GNN框架用于网络药理学节点筛选
2. 证明邻域聚合信息对节点重要性评估的价值
3. 提供可复现的消融实验流程

---

## 8. 实现计划

- [ ] 1. 数据准备：构建PyG格式的图数据
- [ ] 2. 特征工程：计算节点中心性特征
- [ ] 3. Baseline实现：Degree Centrality排序
- [ ] 4. GNN实现：GraphSAGE + Attention
- [ ] 5. 半监督训练：通路伪标签
- [ ] 6. 对比实验：4组方法
- [ ] 7. 富集分析：KEGG/GO分析
- [ ] 8. 可视化：网络图、散点图、气泡图
