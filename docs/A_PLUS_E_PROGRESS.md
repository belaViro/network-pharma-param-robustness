# A+E 探索记录归档

> 最后更新：2026-05-04  
> 当前状态：历史探索，不进入当前主线。

## 1. 归档原因

早期 A+E 方案尝试将参数敏感性分析与 Bootstrap 不确定性量化结合，用于改进传统网络药理学流程。后续项目主线已经调整为更清晰、证据链更完整的：

> 参数扰动 + 成分扰动的 dual-perturbation robustness framework。

因此，本文件不再作为当前进度或论文口径来源。当前主结果以 `result_lock_v1.md`、`result_lock_robustness_validation_addendum.md`、`研究方法与过程.md` 第十九至二十一章为准。

## 2. 当前主线替代方案

| 旧 A+E 设想 | 当前主线 |
|:---|:---|
| 参数敏感性 + Bootstrap | 112 组参数网格 + 18 个化合物 LOO |
| 固定参数 vs A+E 最优参数 | 参数 Top30 消融：固定参数、单一最佳参数、frequency Top30、Q1 robust Top30 |
| 数据库 overlap 验证 | PubMed 共现、超几何/FDR、degree-matched permutation |
| A+E 消融 | 参数稳健消融 + 成分扰动消融分开做 |

## 3. 当前锁定数字

| 项目 | 当前值 |
|:---|---:|
| 参数组合 | 112 |
| canonical 化合物 | 18 |
| 成分 LOO 组合 | 2016 |
| Q1 candidate robust targets | 50 |
| dual-robust targets | 29 |
| dual-robust functional terms | 58 |
| Tier 1 dual-robust target | 1 |

## 4. 不再使用的旧口径

以下内容仅作为历史探索，不进入当前论文主文、补充表或汇报：

- A+E 作为主方法名；
- 早期固定小规模 Bootstrap 计划；
- 旧外部数据库计划作为 P0；
- 旧“5 个草莓叶成分”的进度描述；
- 任何未进入 112 组合 + 成分 LOO 主线的结果。

## 5. 当前汇报口径

如果需要解释 A+E 探索，可以写：

> 早期曾尝试参数敏感性与 Bootstrap 不确定性量化结合，但后续评估认为该路线会增加方法复杂度且不如成分扰动直观。因此项目主线已转向双扰动稳健性框架：通过 112 组参数网格评估方法参数稳定性，再通过 18 个化合物的 leave-one-compound-out 分析评估成分集合稳定性，最终筛选 29 个 dual-robust targets 和 58 个 dual-robust functional terms。
