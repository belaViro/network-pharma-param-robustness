"""
Generate PyMOL scripts (.pml) and a batch runner for the 9 core docking pairs.
"""
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
RECEPTORS_DIR = BASE_DIR / "receptors"
LIGANDS_DIR = BASE_DIR / "ligands"
PYMOL_DIR = BASE_DIR / "pymol_scripts"
PYMOL_DIR.mkdir(exist_ok=True)

combinations = [
    ('STAT3', 'Catechin', -6.265),
    ('STAT3', 'Miquelianin', -7.234),
    ('STAT3', 'Procyanidin_B1', -7.024),
    ('EGFR', 'Catechin', -8.267),
    ('EGFR', 'Miquelianin', -9.499),
    ('EGFR', 'Procyanidin_B1', -7.795),
    ('SRC', 'Catechin', -9.160),
    ('SRC', 'Miquelianin', -9.764),
    ('SRC', 'Procyanidin_B1', -8.310),
]


def generate_pml(receptor, ligand, affinity):
    """Generate a PyMOL script for one combination."""
    
    rec_pdbqt = (RECEPTORS_DIR / f"{receptor}.pdbqt").as_posix()
    lig_pdbqt = (RESULTS_DIR / f"{receptor}_{ligand}_docked.pdbqt").as_posix()
    output = (RESULTS_DIR / f"{receptor}_{ligand}_pymol.png").as_posix()
    
    pml = f'''# PyMOL Script: {receptor} x {ligand}
# Affinity: {affinity:.2f} kcal/mol
# Run in PyMOL: File -> Run Script -> select this .pml file

# --- Load structures ---
cmd.load("{rec_pdbqt}", "protein")
cmd.load("{lig_pdbqt}", "ligand")

# --- Remove water and hetero atoms ---
cmd.remove("solvent")
cmd.remove("resn HOH")

# --- Protein style: Cartoon ---
cmd.hide("everything", "protein")
cmd.show("cartoon", "protein")
cmd.color("skyblue", "protein")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_highlight_color", "marine")

# --- Ligand style: Purple sticks (like reference image) ---
cmd.show("sticks", "ligand")
cmd.color("magenta", "ligand")
cmd.set("stick_radius", 0.15)
cmd.set("stick_ball", 0)  # No ball ends = pure stick model
cmd.set("stick_ball_ratio", 1.5)

# --- Select pocket residues (within 5A of ligand) ---
cmd.select("pocket", "protein within 5 of ligand")
cmd.select("pocket_residues", "byres pocket")

# Show pocket residues as thin gray sticks
cmd.show("sticks", "pocket_residues")
cmd.color("gray60", "pocket_residues")
cmd.set("stick_radius", 0.10, "pocket_residues")

# --- Select and highlight key residues (closest to ligand) ---
# These will be labeled and shown with thicker sticks
cmd.select("key_res", "pocket_residues")
cmd.set("stick_radius", 0.15, "key_res")

# --- Labels for key residues ---
cmd.label("name CA and pocket_residues", "'%s-%s' % (resn, resi)")
cmd.set("label_size", 20)
cmd.set("label_color", "black")
cmd.set("label_font_id", 7)  # Arial

# --- Hydrogen bonds (yellow dashed lines) ---
cmd.distance("hbonds", "ligand and (elem N+O)", "pocket_residues and (elem N+O)", 3.5)
cmd.set("dash_color", "yellow", "hbonds")
cmd.set("dash_width", 3.0)
cmd.set("dash_radius", 0.05)
cmd.hide("labels", "hbonds")  # Hide distance labels, just show lines

# --- Camera and view ---
cmd.zoom("ligand", 10)
cmd.center("ligand")

# --- White background ---
cmd.bg_color("white")
cmd.set("ray_shadows", 0)

# --- High quality rendering ---
cmd.set("ray_trace_mode", 1)  # High quality
cmd.set("antialias", 2)
cmd.set("ray_texture", 1)
cmd.set("opaque_background", 1)

# --- Save image ---
cmd.ray(2400, 1800)
cmd.png("{output}", dpi=300)

print("Image saved: {output}")
'''
    return pml


# Generate individual .pml files
print("Generating PyMOL scripts...")
for receptor, ligand, affinity in combinations:
    pml_content = generate_pml(receptor, ligand, affinity)
    pml_file = PYMOL_DIR / f"{receptor}_{ligand}.pml"
    with pml_file.open('w', encoding='utf-8') as f:
        f.write(pml_content)
    print(f"  {pml_file.name}")

# Generate batch runner (runs all 9)
batch_pml = '''# PyMOL Batch Script: Run all 9 docking visualizations
# Usage: File -> Run Script -> select batch_run_all.pml

'''
for receptor, ligand, affinity in combinations:
    batch_pml += f'run "{receptor}_{ligand}.pml"\n'

batch_file = PYMOL_DIR / "batch_run_all.pml"
with batch_file.open('w', encoding='utf-8') as f:
    f.write(batch_pml)
print(f"\n  {batch_file.name} (runs all 9)")

# Generate README
readme = '''# PyMOL 分子对接可视化脚本

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
'''

readme_file = PYMOL_DIR / "README.md"
with readme_file.open('w', encoding='utf-8') as f:
    f.write(readme)

print(f"\n{'='*50}")
print(f"PyMOL脚本已生成到: {PYMOL_DIR}")
print(f"{'='*50}")
print(f"共 {len(combinations)} 个脚本 + 1 个批量运行脚本")
print(f"\n使用步骤:")
print(f"1. 下载 PyMOL 教育版: https://pymol.org/edu/")
print(f"2. 打开 PyMOL → File → Run Script → 选择 batch_run_all.pml")
print(f"3. 图片自动保存到 results/ 目录")
