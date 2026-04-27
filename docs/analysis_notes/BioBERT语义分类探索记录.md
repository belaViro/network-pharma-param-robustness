# BioBERT 语义通路分类探索记录

**日期**: 2026-04-26  
**背景**: 参数敏感性分析中，现有通路分类体系（GO 层级追溯 + 关键词匹配）存在严重假阴性。尝试引入 BioBERT 语义理解作为补充验证方法。

---

## 一、现有分类体系的问题

当前 `src/sensitivity/classification.py` 采用三级优先级分类：

| 层级 | 方法 | 问题 |
|------|------|------|
| `direct_ox` | GO ID 层级追溯（祖先链匹配） | 只能捕获氧化应激响应本体子树内的通路，对信号转导类通路无能为力 |
| `cross_talk` | 10 条人工精选 KEGG 精确匹配 | 覆盖太窄，KEGG 有 300+ 通路 |
| `signaling_backbone` | 12 个关键词正则匹配 | 假阴性极高，如 `Pathways in cancer` 语义相关但无关键词 |

**核心矛盾**: 用户的 31 个 Q1 稳健靶点（EGFR、MAPK14、STAT3 等激酶/受体类）富集出的通路主要是**信号转导和磷酸化调控**，这些通路在 GO DAG 中与氧化应激分支平行，层级追溯无法归类；关键词匹配又漏掉了大量语义相关通路。

---

## 二、探索方案与测试结果

### 方案 1：BioBERT [CLS] 向量 + 余弦相似度

**实现**: 使用 `dmis-lab/biobert-base-cased-v1.1`，取 `[CLS]` token 的 hidden state 作为句子向量，计算与 `"oxidative stress"` 的余弦相似度。

**测试结果**:

| 通路 | 分数 | 是否合理 |
|------|------|---------|
| oxidative stress | 1.000 | ✅ |
| Pathways in cancer | 0.871 | 相关，分数合理 |
| MAPK signaling pathway | 0.885 | 相关，分数合理 |
| ABC transporters | 0.891 | ❌ 不相关，分数过高 |
| Ribosome | 0.864 | ❌ 不相关，分数过高 |
| RNA polymerase | 0.842 | ❌ 不相关，分数过高 |

**结论**: 分数分布过于集中（0.84–0.96），**无法区分相关与不相关通路**。BioBERT 基础编码器未针对句子相似度任务训练，[CLS] 向量不适合直接做零样本判别。

---

### 方案 2：BioBERT Mean Pooling + 余弦相似度

**改进**: 将所有有效 token 的向量取平均（排除 padding），替代 [CLS]。

**测试结果**:

| 通路 | 分数 |
|------|------|
| oxidative stress | 1.000 |
| response to oxidative stress | 0.961 |
| Pathways in cancer | 0.814 |
| MAPK signaling pathway | 0.770 |
| Ferroptosis | 0.783 |
| Ribosome | 0.778 |
| RNA polymerase | 0.748 |

**结论**: 分数范围扩大到 0.74–0.96，但**仍然无法有效区分**。相关通路与不相关通路分数重叠严重。

---

### 方案 3：通用 Sentence-BERT (all-MiniLM-L6-v2) + 句子模板

**改进**: 使用专门为语义相似度训练的 `sentence-transformers/all-MiniLM-L6-v2`（384 维，~80MB），并将输入包装为完整句子：
- 参考: `"This biological pathway is involved in oxidative stress"`
- 查询: `"This biological pathway is involved in {pathway_name}"`

**测试结果**:

| 通路 | 分数 | 判定 (threshold=0.55) |
|------|------|----------------------|
| oxidative stress | 1.000 | ✅ OX |
| Pathways in cancer | 0.599 | ✅ OX |
| Fluid shear stress and atherosclerosis | 0.612 | ✅ OX |
| positive regulation of phosphorylation | 0.574 | ✅ OX |
| cytokine-mediated signaling | 0.579 | ✅ OX |
| Ferroptosis | 0.489 | ❌ 漏掉（假阴性）|
| MAPK signaling | 0.438 | ❌ 漏掉 |
| Biotin metabolism | 0.596 | ❌ 误中（假阳性）|
| Ribosome | 0.474 | ✅ 正确排除 |

**结论**: 句子模板显著改善了区分度，但仍存在**关键假阴性**（`Ferroptosis`，铁死亡——氧化应激核心通路）和**假阳性**（`Biotin metabolism`）。通用模型缺乏生物医学领域的精确语义理解。

---

### 方案 4：生物医学 Sentence-BERT (S-PubMedBert-MS-MARCO) + 句子模板

**改进**: 使用生物医学专用模型 `pritamdeka/S-PubMedBert-MS-MARCO`（768 维，~400MB，基于 PubMedBERT + MS MARCO 训练）。

**测试结果**:

| 通路 | 分数 | 判定 (threshold=0.55) |
|------|------|----------------------|
| oxidative stress | 1.000 | ✅ OX |
| Ferroptosis | 0.928 | ✅ OX |
| Pathways in cancer | 0.931 | ✅ OX |
| MAPK signaling | 0.911 | ✅ OX |
| Ribosome | 0.914 | ❌ **误中** |
| Spliceosome | 0.928 | ❌ **误中** |
| RNA polymerase | 0.917 | ❌ **误中** |
| Biotin metabolism | 0.912 | ❌ **误中** |

**结论**: **灾难性失效**。该模型在 MS MARCO（检索排序）任务上训练，倾向于给所有查询高相似度分数。**所有通路（包括明显不相关的）都 ≥ 0.89**，完全丧失区分能力。

---

## 三、根因分析

为什么四种方案全部失败？

### 1. 输入文本过短
通路名称仅 2–5 个词（如 `"Ferroptosis"`、`"MAPK signaling"`），缺乏足够上下文让模型理解其生物学功能。

### 2. 零样本相似度任务的固有局限
预训练语言模型学习的是"词共现统计"，而非"生物学因果关系"。模型无法判断：
- `Ferroptosis` 与氧化应激是**因果关联**
- `Ribosome` 与氧化应激只是**间接共现**（因为所有细胞都有核糖体）

### 3. 缺乏判别性训练
所有测试模型均未接受过"通路-疾病相关性判别"任务的训练。句子相似度 ≠ 生物学相关性。

### 4. 领域鸿沟
即使是生物医学预训练模型（BioBERT、PubMedBERT），其训练目标（Masked LM、NSP）也不包含"判断两条通路是否功能相关"这一任务。

---

## 四、代码实现记录

尽管语义分类效果不达标，相关代码已保留在项目中，可作为后续改进的基础：

| 文件 | 说明 |
|------|------|
| `src/sensitivity/biobert_classifier.py` | BioBERT / Sentence-BERT 语义分类器封装，含懒加载、批量编码、缓存机制 |
| `src/sensitivity/classification.py`（已修改）| 新增 `classify_pathway_with_biobert()` 和 `classify_pathways_batch_with_biobert()`，支持规则 + BioBERT 融合分类 |
| `src/sensitivity/__init__.py`（已修改）| 导出新增函数 |

**使用方式**（如后续找到更合适的模型，可直接替换）：
```python
from src.sensitivity.classification import classify_pathway_with_biobert

result = classify_pathway_with_biobert("Pathways in cancer")
# {
#   "keyword_category": "other",
#   "biobert_score": 0.82,
#   "biobert_category": "biobert_ox",
#   "final_category": "biobert_ox"
# }
```

---

## 五、结论与转向建议

**结论**: 纯零样本语义相似度（无论是 BioBERT [CLS]、Mean Pooling、通用 Sentence-BERT 还是生物医学 Sentence-BERT）**均无法可靠地完成通路-疾病相关性判别**。该方案在当前技术条件下不可行。

**推荐转向**: **外部数据库验证**（OMIM / DisGeNET / CTD）

| 对比维度 | BioBERT 语义匹配 | 外部数据库验证 |
|---------|-----------------|---------------|
| 可靠性 | 低（模型幻觉、分数不稳定） | **高（基于 curated 生物医学数据）** |
| 审稿人认可度 | 低（方法新颖但验证标准不透明） | **高（独立第三方数据库，无法质疑）** |
| 实现难度 | 中（模型下载、推理、阈值调优） | 低（API 查询或本地数据文件 + 超几何检验）|
| 论文卖点 | "NLP 创新"（但效果不佳） | "多重独立验证"（方法学严谨） |
| 工作量 | 0.5–1 天（已尝试） | 1–2 天 |

**下一步行动**:
1. 下载/获取 OMIM、DisGeNET、CTD 的氧化应激相关基因列表
2. 对 31 个 Q1 稳健靶点做超几何检验富集分析
3. 产出验证表格 + 森林图，作为论文 Methods 中的"External Validation"章节

---

## 六、参考文档

- `docs/analysis_notes/审稿人质疑与应对策略.md` — 预判了审稿人对验证标准的质疑
- `docs/analysis_notes/NLP通路分类解释.md` — BioBERT 原理说明
- `docs/analysis_notes/网药ML创新点子.md` — BioBERT 被列为推荐创新点之一

---

> **一句话总结**: BioBERT 语义分类在理论上美好，实践中因输入过短、缺乏判别训练、领域鸿沟而失效。转向外部数据库验证是更务实、更可靠、审稿人更认可的选择。
