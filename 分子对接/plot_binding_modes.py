"""
Generate 3D binding mode visualization from Vina docking results.
No PyMOL required — pure Python with matplotlib 3D.
"""
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from collections import defaultdict

RESULTS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\results'
RECEPTORS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\receptors'
LIGANDS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\ligands'
OUTPUT_DIR = RESULTS_DIR
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Color palette for elements
ELEM_COLOR = {
    'C': '#444444', 'N': '#3050F8', 'O': '#FF0D0D', 'S': '#FFFF30',
    'H': '#EEEEEE', 'P': '#FF8000', 'F': '#90E050', 'CL': '#1FF01F',
    'BR': '#A62929', 'I': '#940094', 'B': '#FFB5B5'
}
VDW_RADII = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2, 'P': 1.8}


def parse_pdbqt_atoms(filepath):
    """Parse PDBQT file, return list of atom dicts."""
    atoms = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': line[21].strip(),
                    'resseq': int(line[22:26].strip()) if line[22:26].strip() else 0,
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[77:79].strip().upper() if len(line) > 78 else line[12:16].strip()[0].upper(),
                    'charge': float(line[70:76].strip()) if len(line) > 75 and line[70:76].strip() else 0.0,
                    'atype': line[77:79].strip() if len(line) > 78 else '',
                }
                # Fix element from atom name
                if atom['element'] not in ELEM_COLOR:
                    atom['element'] = atom['name'][0].upper()
                atoms.append(atom)
    return atoms


def parse_ligand_poses(filepath):
    """Parse multi-model PDBQT ligand, return list of poses (each is atom list)."""
    poses = []
    current = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('MODEL'):
                current = []
            elif line.startswith('ENDMDL'):
                if current:
                    poses.append(current)
                current = []
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[77:79].strip().upper() if len(line) > 78 else line[12:16].strip()[0].upper(),
                    'charge': float(line[70:76].strip()) if len(line) > 75 and line[70:76].strip() else 0.0,
                    'atype': line[77:79].strip() if len(line) > 78 else '',
                }
                if atom['element'] not in ELEM_COLOR:
                    atom['element'] = atom['name'][0].upper()
                current.append(atom)
        if current and not poses:
            poses.append(current)
    return poses


def extract_pocket_residues(receptor_atoms, ligand_atoms, cutoff=5.0):
    """Extract receptor residues within cutoff of any ligand atom."""
    ligand_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    pocket_atoms = []
    pocket_residues = defaultdict(list)
    
    for ratom in receptor_atoms:
        rcoord = np.array([ratom['x'], ratom['y'], ratom['z']])
        dists = np.linalg.norm(ligand_coords - rcoord, axis=1)
        if np.min(dists) <= cutoff:
            pocket_atoms.append(ratom)
            key = (ratom['resname'], ratom['resseq'], ratom['chain'])
            pocket_residues[key].append(ratom)
    
    return pocket_atoms, pocket_residues


def find_hbonds(ligand_atoms, receptor_atoms, max_dist=3.2):
    """Find potential H-bonds between ligand and receptor."""
    # Simple H-bond detection: H-bond donor (N-H, O-H) to acceptor (N, O)
    # More accurate: check atom types in PDBQT (A=acceptor, HD=donor H, NA=aromatic N, etc.)
    hbonds = []
    
    # Define acceptors and donors by element
    acceptors = ['O', 'N']
    donors_H = ['H']
    
    for latom in ligand_atoms:
        for ratom in receptor_atoms:
            if ratom['element'] == 'H':
                continue
            dist = np.linalg.norm(
                np.array([latom['x'], latom['y'], latom['z']]) -
                np.array([ratom['x'], ratom['y'], ratom['z']])
            )
            if dist <= max_dist:
                # Check if one is donor-H and other is acceptor
                lelem = latom['element']
                relem = ratom['element']
                
                # Simplified: O/N to O/N contacts within H-bond distance
                if lelem in acceptors and relem in acceptors:
                    hbonds.append({
                        'latom': latom, 'ratom': ratom,
                        'dist': dist, 'type': 'HB'
                    })
    
    # Sort by distance, take closest unique pairs
    hbonds.sort(key=lambda x: x['dist'])
    seen = set()
    filtered = []
    for hb in hbonds:
        key = (hb['latom']['serial'], hb['ratom']['serial'])
        if key not in seen:
            seen.add(key)
            filtered.append(hb)
    return filtered[:8]  # Limit to top 8


def draw_binding_mode(receptor_file, ligand_file, receptor_name, ligand_name, 
                      affinity, output_file):
    """Draw a 3D binding mode figure."""
    
    receptor_atoms = parse_pdbqt_atoms(receptor_file)
    poses = parse_ligand_poses(ligand_file)
    if not poses:
        print(f"No poses found for {receptor_name}x{ligand_name}")
        return
    
    ligand_atoms = poses[0]  # Best pose (first model)
    
    # Extract pocket
    pocket_atoms, pocket_residues = extract_pocket_residues(
        receptor_atoms, ligand_atoms, cutoff=5.5
    )
    
    # Find H-bonds
    hbonds = find_hbonds(ligand_atoms, pocket_atoms)
    
    # Create figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # --- Draw receptor pocket (wireframe style) ---
    rec_coords = np.array([[a['x'], a['y'], a['z']] for a in pocket_atoms])
    if len(rec_coords) > 0:
        # Draw C-alpha trace
        ca_atoms = [a for a in pocket_atoms if a['name'] == 'CA']
        if len(ca_atoms) >= 2:
            ca_coords = np.array([[a['x'], a['y'], a['z']] for a in ca_atoms])
            # Sort by residue sequence
            ca_sorted = sorted(ca_atoms, key=lambda a: a['resseq'])
            ca_coords_sorted = np.array([[a['x'], a['y'], a['z']] for a in ca_sorted])
            ax.plot(ca_coords_sorted[:, 0], ca_coords_sorted[:, 1], ca_coords_sorted[:, 2],
                   'k-', alpha=0.3, linewidth=0.8)
        
        # Draw pocket atoms as small dots
        for atom in pocket_atoms:
            if atom['element'] == 'C':
                ax.scatter(atom['x'], atom['y'], atom['z'],
                          c='#AAAAAA', s=8, alpha=0.4, edgecolors='none')
    
    # --- Draw ligand (ball-and-stick) ---
    lig_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    lig_center = lig_coords.mean(axis=0)
    
    # Shift view to ligand center
    for atom in ligand_atoms:
        color = ELEM_COLOR.get(atom['element'], '#888888')
        size = VDW_RADII.get(atom['element'], 1.5) * 80
        ax.scatter(atom['x'], atom['y'], atom['z'],
                  c=color, s=size, edgecolors='black', linewidths=0.5, alpha=0.9, zorder=10)
    
    # Draw ligand bonds
    # Simple distance-based bond detection
    lig_bonds = []
    for i, a1 in enumerate(ligand_atoms):
        for j, a2 in enumerate(ligand_atoms[i+1:], i+1):
            d = np.linalg.norm(
                np.array([a1['x'], a1['y'], a1['z']]) -
                np.array([a2['x'], a2['y'], a2['z']])
            )
            # Typical bond distances
            max_bond = 1.8
            if a1['element'] == 'H' or a2['element'] == 'H':
                max_bond = 1.3
            if d <= max_bond:
                lig_bonds.append((a1, a2))
                ax.plot([a1['x'], a2['x']], [a1['y'], a2['y']], [a1['z'], a2['z']],
                       'k-', linewidth=1.5, alpha=0.8, zorder=9)
    
    # --- Draw H-bonds (dashed yellow lines) ---
    for hb in hbonds:
        a1, a2 = hb['latom'], hb['ratom']
        ax.plot([a1['x'], a2['x']], [a1['y'], a2['y']], [a1['z'], a2['z']],
               'y--', linewidth=1.5, alpha=0.8, zorder=15)
    
    # --- Label key residues ---
    # Pick representative atom per residue (CA or first atom)
    labeled = []
    for (resname, resseq, chain), atoms in sorted(pocket_residues.items(), key=lambda x: x[0][1]):
        if len(atoms) == 0:
            continue
        # Find atom closest to ligand center
        closest = min(atoms, key=lambda a: np.linalg.norm(
            np.array([a['x'], a['y'], a['z']]) - lig_center
        ))
        # Only label if within reasonable distance and not too many labels
        dist_to_lig = np.linalg.norm(
            np.array([closest['x'], closest['y'], closest['z']]) - lig_center
        )
        if dist_to_lig < 6.0:
            label_text = f"{resname}{resseq}"
            # Offset label slightly
            offset = 0.8
            ax.text(closest['x'] + offset, closest['y'] + offset, closest['z'] + offset,
                   label_text, fontsize=7, color='#444444', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.15', facecolor='white', 
                            edgecolor='gray', alpha=0.7))
            labeled.append(label_text)
    
    # Set view to ligand center
    margin = 8
    ax.set_xlim(lig_center[0] - margin, lig_center[0] + margin)
    ax.set_ylim(lig_center[1] - margin, lig_center[1] + margin)
    ax.set_zlim(lig_center[2] - margin, lig_center[2] + margin)
    
    # Remove axis
    ax.set_axis_off()
    
    # Title
    title = f"{receptor_name} × {ligand_name}\nBinding Affinity: {affinity:.2f} kcal/mol"
    fig.suptitle(title, fontsize=13, fontweight='bold', y=0.95)
    
    # Legend / info box
    info_text = f"H-bonds: {len(hbonds)} | Pocket residues: {len(pocket_residues)}"
    fig.text(0.5, 0.02, info_text, ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.92])
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {os.path.basename(output_file)}")


# === MAIN ===
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

print("Generating binding mode visualizations...")
for receptor, ligand, affinity in combinations:
    rec_file = os.path.join(RECEPTORS_DIR, f"{receptor}.pdbqt")
    lig_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_docked.pdbqt")
    out_file = os.path.join(OUTPUT_DIR, f"{receptor}_{ligand}_binding_mode.png")
    
    if not os.path.exists(lig_file):
        print(f"  SKIP: {ligand_file} not found")
        continue
    
    print(f"\nProcessing: {receptor} × {ligand} ({affinity:.2f})")
    draw_binding_mode(rec_file, lig_file, receptor, ligand, affinity, out_file)

print(f"\n✓ All figures saved to: {OUTPUT_DIR}")
