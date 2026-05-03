import os
import subprocess
import json
import re

VINA_EXE = r'D:\科研\网络药理学参数敏感性分析\分子对接\vina_1.2.7_win.exe'
BASE_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接'
RECEPTORS_DIR = os.path.join(BASE_DIR, 'receptors')
LIGANDS_DIR = os.path.join(BASE_DIR, 'ligands')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

receptors = {
    'STAT3': {'pdbqt': os.path.join(RECEPTORS_DIR, 'STAT3.pdbqt'), 'center': (13.028, 55.637, 0.250)},
    'EGFR':  {'pdbqt': os.path.join(RECEPTORS_DIR, 'EGFR.pdbqt'),  'center': (22.014, 0.253, 52.794)},
    'SRC':   {'pdbqt': os.path.join(RECEPTORS_DIR, 'SRC.pdbqt'),   'center': (20.255, 20.924, 58.253)},
}
ligands = {
    'Catechin': os.path.join(LIGANDS_DIR, 'Catechin.pdbqt'),
    'Miquelianin': os.path.join(LIGANDS_DIR, 'Miquelianin.pdbqt'),
    'Procyanidin_B1': os.path.join(LIGANDS_DIR, 'Procyanidin_B1.pdbqt'),
}

results_summary = []

for rname, rinfo in receptors.items():
    for lname, lpath in ligands.items():
        out_pdbqt = os.path.join(RESULTS_DIR, f'{rname}_{lname}_docked.pdbqt')
        cmd = [
            VINA_EXE,
            '--receptor', rinfo['pdbqt'],
            '--ligand', lpath,
            '--center_x', str(rinfo['center'][0]),
            '--center_y', str(rinfo['center'][1]),
            '--center_z', str(rinfo['center'][2]),
            '--size_x', '20', '--size_y', '20', '--size_z', '20',
            '--out', out_pdbqt,
            '--cpu', '4',
            '--num_modes', '9',
            '--exhaustiveness', '16',
        ]
        print(f"\n对接: {rname} × {lname}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        stdout = result.stdout

        best_affinity = None
        for line in stdout.split('\n'):
            m = re.match(r'^\s*1\s+(-?\d+\.\d+)', line.strip())
            if m:
                best_affinity = float(m.group(1))
                break

        results_summary.append({
            'receptor': rname, 'ligand': lname,
            'affinity_kcal_mol': best_affinity,
            'output': out_pdbqt if os.path.exists(out_pdbqt) else None,
        })
        aff_str = f"{best_affinity:.2f}" if best_affinity is not None else "N/A"
        print(f"  → {aff_str} kcal/mol")

# 保存汇总
with open(os.path.join(RESULTS_DIR, 'docking_summary.json'), 'w') as f:
    json.dump(results_summary, f, indent=2)

print(f"\n{'='*55}")
print("分子对接结果汇总")
print(f"{'='*55}")
print(f"{'受体':<10} {'配体':<18} {'结合能(kcal/mol)':>16}")
print("-" * 55)
for r in results_summary:
    aff = f"{r['affinity_kcal_mol']:.2f}" if r['affinity_kcal_mol'] is not None else "N/A"
    print(f"{r['receptor']:<10} {r['ligand']:<18} {aff:>16}")
