#!/usr/bin/env bash
# ============================================================================
# 网络药理学参数敏感性分析 — 复现脚本
# ============================================================================
# 用法: bash reproduce.sh
# 要求: Python >= 3.10, pip
#
# 这将按顺序执行:
#   1. 参数网格分析（72个组合，靶点稳健性，通路稳定性）
#   2. 稳健靶点 vs 各参数组合的富集效果对比
#   3. 通路象限增强与论文图重新生成
#
# 所有输出保存在 src/outputs/sensitivity/
# ============================================================================
set -euo pipefail

echo "=============================================="
echo "  网络药理学参数敏感性分析 — 复现"
echo "=============================================="

# 切换到项目根目录
cd "$(dirname "$0")"

# 检查 Python 版本
PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")' 2>/dev/null || echo "0.0")
REQUIRED="3.10"
if [ "$(printf '%s\n' "$REQUIRED" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED" ]; then
    echo "错误: 需要 Python >= 3.10 (当前: $PYTHON_VERSION)"
    exit 1
fi

# 1. 安装依赖
echo ""
echo "[1/5] 安装依赖..."
pip install -e . --quiet

# 2. 运行主分析
echo ""
echo "[2/5] 运行参数网格分析..."
python scripts/analysis/run_sensitivity_analysis.py

# 3. 运行对比分析
echo ""
echo "[3/5] 稳健靶点 vs 参数组合对比..."
python scripts/analysis/compare_robust_vs_parametric.py

# 4. 追加通路稳健性象限
echo ""
echo "[4/5] 通路稳健性象限增强..."
python scripts/analysis/augment_pathway_robustness.py

# 5. 重新生成主要图件
echo ""
echo "[5/5] 重新生成图件..."
python scripts/analysis/plot_sensitivity_figures.py
python scripts/analysis/plot_pathway_quadrant.py
python scripts/analysis/plot_pathway_quadrant_by_source.py
python scripts/analysis/plot_param_space_fig10.py
python scripts/analysis/plot_upset_fig8.py
python scripts/analysis/plot_network_fig9.py

# 完成
echo ""
echo "=============================================="
echo "  复现完成!"
echo "  输出目录: src/outputs/sensitivity/"
echo "=============================================="
