import json
import re
import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
VINA_EXE = BASE_DIR / 'vina_1.2.7_win.exe'
RECEPTORS_DIR = BASE_DIR / 'receptors'
LIGANDS_DIR = BASE_DIR / 'ligands'
RESULTS_DIR = BASE_DIR / 'results'
RESULTS_DIR.mkdir(exist_ok=True)

EXHAUSTIVENESS = 16
CPU = 4
NUM_MODES = 9
BOX_SIZE = (20, 20, 20)

receptors = {
    'STAT3': {'pdbqt': RECEPTORS_DIR / 'STAT3.pdbqt', 'center': (13.028, 55.637, 0.250), 'pdb': '6NJS', 'domain': 'SH2'},
    'EGFR':  {'pdbqt': RECEPTORS_DIR / 'EGFR.pdbqt',  'center': (22.014, 0.253, 52.794), 'pdb': '1M17', 'domain': 'kinase'},
    'SRC':   {'pdbqt': RECEPTORS_DIR / 'SRC.pdbqt',   'center': (20.255, 20.924, 58.253), 'pdb': '2H8H', 'domain': 'kinase'},
}
ligands = {
    'Catechin': LIGANDS_DIR / 'Catechin.pdbqt',
    'Miquelianin': LIGANDS_DIR / 'Miquelianin.pdbqt',
    'Procyanidin_B1': LIGANDS_DIR / 'Procyanidin_B1.pdbqt',
}

results_summary = []

for rname, rinfo in receptors.items():
    for lname, lpath in ligands.items():
        out_pdbqt = RESULTS_DIR / f'{rname}_{lname}_docked.pdbqt'
        cmd = [
            str(VINA_EXE),
            '--receptor', str(rinfo['pdbqt']),
            '--ligand', str(lpath),
            '--center_x', str(rinfo['center'][0]),
            '--center_y', str(rinfo['center'][1]),
            '--center_z', str(rinfo['center'][2]),
            '--size_x', str(BOX_SIZE[0]), '--size_y', str(BOX_SIZE[1]), '--size_z', str(BOX_SIZE[2]),
            '--out', str(out_pdbqt),
            '--cpu', str(CPU),
            '--num_modes', str(NUM_MODES),
            '--exhaustiveness', str(EXHAUSTIVENESS),
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
            'output': str(out_pdbqt.relative_to(BASE_DIR)) if out_pdbqt.exists() else None,
        })
        aff_str = f"{best_affinity:.2f}" if best_affinity is not None else "N/A"
        print(f"  → {aff_str} kcal/mol")

summary = {
    'tool': 'AutoDock Vina v1.2.7',
    'settings': {
        'exhaustiveness': EXHAUSTIVENESS,
        'cpu': CPU,
        'num_modes': NUM_MODES,
        'grid_size': f'{BOX_SIZE[0]}x{BOX_SIZE[1]}x{BOX_SIZE[2]} Å',
        'spacing': 0.375,
    },
    'receptors': {
        name: {
            'pdb': info['pdb'],
            'domain': info['domain'],
            'center': info['center'],
        }
        for name, info in receptors.items()
    },
    'results': results_summary,
}

with (RESULTS_DIR / 'docking_summary.json').open('w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(f"\n{'='*55}")
print("分子对接结果汇总")
print(f"{'='*55}")
print(f"{'受体':<10} {'配体':<18} {'结合能(kcal/mol)':>16}")
print("-" * 55)
for r in results_summary:
    aff = f"{r['affinity_kcal_mol']:.2f}" if r['affinity_kcal_mol'] is not None else "N/A"
    print(f"{r['receptor']:<10} {r['ligand']:<18} {aff:>16}")
