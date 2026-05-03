import os
import re
import urllib.request

def clean_pdb(input_pdb, output_pdb):
    """去除水分子、配体、杂原子，保留ATOM记录"""
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    cleaned = []
    for line in lines:
        if line.startswith('ATOM  ') or line.startswith('TER   ') or line.startswith('END'):
            cleaned.append(line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(cleaned)
    return output_pdb

def get_ligand_center(pdb_file, ligand_resname):
    """从PDB中提取指定配体的几何中心"""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if (line.startswith('HETATM') or line.startswith('ATOM  ')) and ligand_resname in line[17:20]:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except:
                    pass
    if not coords:
        return None
    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    return (cx, cy, cz)

def guess_ligand_resname(pdb_file):
    """猜测配体残基名（最常见的药物分子）"""
    hetatms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                resname = line[17:20].strip()
                if resname not in ['HOH', 'DOD', 'WAT', 'H2O', ' CL', ' NA', ' MG', ' CA', ' ZN', ' K ']:
                    hetatms.append(resname)
    from collections import Counter
    if hetatms:
        return Counter(hetatms).most_common(1)[0][0]
    return None

# 配置
proteins = {
    'STAT3': {'pdb': 'STAT3_6NJS.pdb', 'ligand': None},  # 自动检测
    'EGFR': {'pdb': 'EGFR_1M17.pdb', 'ligand': None},
    'SRC': {'pdb': 'SRC_2H8H.pdb', 'ligand': None},
}

os.chdir(r'D:\科研\网络药理学参数敏感性分析\分子对接\receptors')

for name, info in proteins.items():
    pdb_file = info['pdb']
    if not os.path.exists(pdb_file):
        print(f"{name}: PDB文件不存在")
        continue
    
    # 检测配体
    ligand = guess_ligand_resname(pdb_file)
    print(f"{name}: 检测到配体 = {ligand}")
    
    # 计算口袋中心
    center = get_ligand_center(pdb_file, ligand) if ligand else None
    if center:
        print(f"{name}: 口袋中心 = ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    else:
        print(f"{name}: 无法确定口袋中心")
    
    # 清理PDB
    clean_file = f"{name}_clean.pdb"
    clean_pdb(pdb_file, clean_file)
    
    # 转换为PDBQT（使用Open Babel）
    pdbqt_file = f"{name}.pdbqt"
    obabel = r'C:\Users\liuzhen\AppData\Local\Programs\Python\Python311\Scripts\obabel.exe'
    cmd = f'"{obabel}" {clean_file} -O {pdbqt_file} -xr -p 7.4'
    os.system(cmd)
    
    if os.path.exists(pdbqt_file):
        print(f"{name}: 受体PDBQT已生成 ({os.path.getsize(pdbqt_file)} bytes)")
    else:
        print(f"{name}: 受体PDBQT生成失败")
