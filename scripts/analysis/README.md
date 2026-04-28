# Sensitivity Analysis Scripts

Analysis scripts for the parameter robustness sensitivity analysis.

- `run_sensitivity_analysis.py` — Main sensitivity analysis: 72 parameter combos, target robustness, pathway stability
- `compare_robust_vs_parametric.py` — Compare Q1/Tier1 robust targets vs each parameter combo's enrichment
- `augment_pathway_robustness.py` — Add pathway quadrant labels from the current pathway stability table
- `plot_sensitivity_figures.py` — Generate Fig.1-Fig.7
- `plot_pathway_quadrant.py` / `plot_pathway_quadrant_by_source.py` — Generate pathway quadrant figures

## 用法

```bash
cd code
python scripts/analysis/run_sensitivity_analysis.py  # 运行主分析
python scripts/analysis/compare_robust_vs_parametric.py  # 对比分析
python scripts/analysis/augment_pathway_robustness.py  # 生成通路象限增强表
python scripts/analysis/plot_sensitivity_figures.py  # 生成主要图件
```
