import os
import sys
from rdkit import Chem

# 添加site-packages路径
sys.path.insert(0, r"C:\Users\liuzhen\AppData\Local\Programs\Python\Python311\Lib\site-packages")

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
except ImportError as e:
    print(f"meeko导入失败: {e}")
    sys.exit(1)

compounds = ['Catechin', 'Miquelianin', 'Procyanidin_B1']

for name in compounds:
    sdf_path = f"{name}.sdf"
    if not os.path.exists(sdf_path):
        print(f"{name}: SDF文件不存在")
        continue
    
    mols = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = mols[0]
    if mol is None:
        print(f"{name}: 无法读取SDF")
        continue
    
    try:
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        writer = PDBQTWriterLegacy()
        pdbqt_string, _, _ = writer.write_string(preparator.setup)
        
        pdbqt_path = f"{name}.pdbqt"
        with open(pdbqt_path, "w") as f:
            f.write(pdbqt_string)
        print(f"{name}: PDBQT已生成 ({os.path.getsize(pdbqt_path)} bytes)")
    except Exception as e:
        print(f"{name}: 转换失败 - {e}")
