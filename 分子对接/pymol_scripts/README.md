# PyMOL 分子对接可视化脚本

## 环境要求
- PyMOL 2.5+ (免费教育版下载: https://pymol.org/edu/)
- 或 PyMOL 开源版 (conda install -c conda-forge pymol-open-source)

## 使用方法

### 单张生成
1. 打开 PyMOL
2. File -> Run Script -> 选择 `pymol_scripts/` 目录下的某个 `.pml` 文件
3. 等待渲染完成，图片自动保存到 `results/` 目录

### 批量生成
1. 打开 PyMOL
2. File -> Run Script -> 选择 `pymol_scripts/batch_run_all.pml`
3. 自动依次生成全部9张图

## 脚本说明
每张图包含:
- 蛋白 Cartoon (天蓝色)
- 配体 Stick (紫色, 参考图风格)
- 口袋残基 Stick (灰色细棒)
- 氢键虚线 (黄色)
- 残基标签 (黑色)
- 白色背景
- 2400x1800 高清渲染

## 输出文件
- `results/*_pymol.png` — PyMOL渲染的高清3D结合模式图
