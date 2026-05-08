# Analysis Scripts

Analysis scripts for the parameter and compound dual-perturbation robustness workflow.

- `run_sensitivity_analysis.py` — Main sensitivity analysis: 112 parameter combos, target robustness, pathway stability
- `run_compound_sensitivity_analysis.py` — Leave-one-compound-out compound sensitivity with resumable checkpoints
- `compare_robust_vs_parametric.py` — Compare Q1/Tier1 robust targets vs each parameter combo's enrichment
- `augment_pathway_robustness.py` — Add pathway quadrant labels from the current pathway stability table
- `robustness_validation_bundle.py` — External evidence, Top30 ablation, and degree-matched validation bundle
- `generate_paper_figures.py` — Generate manuscript-facing figures
- `export_supplementary_tables.py` — Export Tables S1-S6
- `plot_sensitivity_figures.py` — Generate core sensitivity figures
- `plot_pathway_quadrant.py` / `plot_pathway_quadrant_by_source.py` — Generate pathway quadrant figures

For the manuscript-level end-to-end workflow, use the repository root entry point:

```bash
python reproduce.py          # full workflow from packaged inputs
python reproduce.py --core   # parameter-robustness workflow only
```

## 用法

```bash
cd code
python scripts/analysis/run_sensitivity_analysis.py  # 运行主分析
python scripts/analysis/run_compound_sensitivity_analysis.py  # 运行成分敏感性分析
python scripts/analysis/compare_robust_vs_parametric.py  # 对比分析
python scripts/analysis/augment_pathway_robustness.py  # 生成通路象限增强表
python scripts/analysis/plot_sensitivity_figures.py  # 生成主要图件
```
