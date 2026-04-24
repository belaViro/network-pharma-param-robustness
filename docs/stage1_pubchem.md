# Stage 1: PubChem获取SMILES

## 概述

从PubChem数据库获取草莓叶活性成分的SMILES结构信息。

## 活性成分

| 成分 | 英文名 | PubChem CID |
|------|--------|-------------|
| 儿茶素 | Catechin | 73160 |
| 原花青素B1 | Procyanidin B1 | 11250133 |
| 野黄芩苷 | Scutellarin | 185617 |
| 蔗糖 | Sucrose | 5988 |
| 海藻糖 | Trehalose | 7425 |

## 输入

- 成分名称列表

## 输出

- SMILES 结果直接通过函数返回，供 Stage 2 使用；如需持久化可自行保存到 `src/outputs/shared/`

## 输出格式

```csv
compound_name,smiles,pubchem_cid
Catechin,O[C@H]1[C@H](O)[C@H](Oc2... ,73160
Procyanidin_B1,O[C@H]1[C@H](O)[C@H](Oc2...,11250133
...
```

## 运行命令

```bash
python -m src.pipeline.stage1_fetch
```

## 状态

✅ 已完成
