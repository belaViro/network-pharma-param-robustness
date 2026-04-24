# Stage 2: SwissTargetPrediction靶点预测

## 概述

使用SwissTargetPrediction网站预测各活性成分的作用靶点。

## 方法

- 使用Playwright自动化浏览器
- 模拟真实用户操作
- 物种选择: Homo sapiens (人类)
- 预测阈值: Probability > 0

## 输入

- Stage 1 获取的SMILES文件

## 输出

每个成分生成单独的CSV文件:
- `src/outputs/shared/swiss_targets/Catechin_targets.csv`
- `src/outputs/shared/swiss_targets/Procyanidin_B1_targets.csv`
- `src/outputs/shared/swiss_targets/Scutellarin_targets.csv`
- `src/outputs/shared/swiss_targets/Sucrose_targets.csv`
- `src/outputs/shared/swiss_targets/Trehalose_targets.csv`

## 输出格式

```csv
Target,Common name,Uniprot ID,ChEMBL ID,Target Class,Probability*,Known actives (3D/2D)
Carbonic anhydrase 1,CA1,P00915,CHEMBL261,Lyase,1.0,8 / 13
...
```

## 运行命令

```bash
python -m src.pipeline.stage2_predict
```

## 状态

✅ 已完成
