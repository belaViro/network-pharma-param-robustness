"""
Generate 2D interaction diagrams (LigPlot-style) using RDKit + PIL.
Shows ligand scaffold with surrounding residue labels and interaction dashed lines.
"""
import os
import re
import math
import io
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

RESULTS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\results'
RECEPTORS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\receptors'
LIGANDS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\ligands'

# Interaction type colors
HBOND_COLOR = (34, 139, 34)       # Green for H-bonds
HYDROPHOBIC_COLOR = (128, 0, 128) # Purple for hydrophobic
PIPI_COLOR = (0, 128, 128)        # Teal for pi-pi
OTHER_COLOR = (100, 100, 100)     # Gray for others

ATOM_COLORS_2D = {
    'C': (50, 50, 50),
    'N': (30, 80, 200),
    'O': (220, 30, 30),
    'H': (150, 150, 150),
    'S': (200, 180, 50),
}


def parse_pdbqt_atoms(filepath):
    atoms = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'resseq': int(line[22:26].strip()) if line[22:26].strip() else 0,
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[77:79].strip().upper() if len(line) > 78 else line[12:16].strip()[0].upper(),
                }
                if atom['element'] not in ['C', 'N', 'O', 'H', 'S', 'P']:
                    atom['element'] = atom['name'][0].upper()
                atoms.append(atom)
    return atoms


def parse_ligand_poses(filepath):
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
                }
                if atom['element'] not in ['C', 'N', 'O', 'H', 'S', 'P']:
                    atom['element'] = atom['name'][0].upper()
                current.append(atom)
        if current and not poses:
            poses.append(current)
    return poses


def get_ligand_smiles(name):
    """Get SMILES for the 3 ligands."""
    smiles_map = {
        'Catechin': 'Oc1cc(O)c2c(c1)O[C@H](c1ccc(O)c(O)c1)[C@@H](O)C2',
        'Miquelianin': 'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)cc(O)c12',
        'Procyanidin_B1': 'Oc1cc(O)c2c(c1)O[C@H](c1ccc(O)c(O)c1)[C@@H](O)C2c1c(O)cc(O)c2c1O[C@H](c1ccc(O)c(O)c1)[C@@H](O)C2',
    }
    return smiles_map.get(name, '')


def find_interactions(ligand_atoms, receptor_atoms):
    """Find interactions between ligand and pocket residues."""
    lig_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    
    # Group receptor atoms by residue
    res_atoms = {}
    for ratom in receptor_atoms:
        key = (ratom['resname'], ratom['resseq'])
        if key not in res_atoms:
            res_atoms[key] = []
        res_atoms[key].append(ratom)
    
    interactions = []
    
    for res_key, ratoms in res_atoms.items():
        resname, resseq = res_key
        if resname in ['', 'UNK', 'UNL']:
            continue
        
        # Find closest atom pair between ligand and this residue
        closest_dist = float('inf')
        closest_latom = None
        closest_ratom = None
        
        for latom in ligand_atoms:
            lcoord = np.array([latom['x'], latom['y'], latom['z']])
            for ratom in ratoms:
                rcoord = np.array([ratom['x'], ratom['y'], ratom['z']])
                d = np.linalg.norm(lcoord - rcoord)
                if d < closest_dist:
                    closest_dist = d
                    closest_latom = latom
                    closest_ratom = ratom
        
        if closest_dist <= 5.0:
            # Classify interaction type
            itype = 'hydrophobic'
            color = HYDROPHOBIC_COLOR
            
            if closest_dist <= 3.5:
                lelem = closest_latom['element']
                relem = closest_ratom['element']
                if lelem in ['O', 'N'] and relem in ['O', 'N']:
                    itype = 'H-bond'
                    color = HBOND_COLOR
                elif lelem == 'S' or relem == 'S':
                    itype = 'H-bond'
                    color = HBOND_COLOR
            
            # Aromatic check
            aromatic_res = ['PHE', 'TYR', 'TRP', 'HIS']
            if resname in aromatic_res:
                # Check if ligand atom is also aromatic (C in ring)
                if closest_latom['element'] == 'C':
                    itype = 'pi-pi'
                    color = PIPI_COLOR
            
            interactions.append({
                'resname': resname,
                'resseq': resseq,
                'distance': closest_dist,
                'type': itype,
                'color': color,
                'latom': closest_latom,
                'ratom': closest_ratom,
            })
    
    # Sort by distance, take top 8
    interactions.sort(key=lambda x: x['distance'])
    return interactions[:10]


def draw_dashed_line(draw, p1, p2, color=(100, 100, 100), dash=5):
    """Draw a dashed line."""
    x1, y1 = int(p1[0]), int(p1[1])
    x2, y2 = int(p2[0]), int(p2[1])
    dist = math.hypot(x2 - x1, y2 - y1)
    if dist < 1:
        return
    n_dashes = max(2, int(dist / dash))
    for i in range(0, n_dashes, 2):
        t1 = i / n_dashes
        t2 = min((i + 1) / n_dashes, 1.0)
        x1d = int(x1 + (x2 - x1) * t1)
        y1d = int(y1 + (y2 - y1) * t1)
        x2d = int(x1 + (x2 - x1) * t2)
        y2d = int(y1 + (y2 - y1) * t2)
        draw.line([(x1d, y1d), (x2d, y2d)], fill=color, width=2)


def draw_rounded_rect(draw, xy, fill, radius=12):
    """Draw a rounded rectangle."""
    x1, y1, x2, y2 = xy
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=(80, 80, 80), width=1)


def draw_2d_interaction(ligand_name, receptor_name, affinity, 
                        ligand_file, receptor_file, output_file,
                        img_size=(900, 750)):
    """Draw a 2D interaction diagram."""
    
    # Parse structures
    receptor_atoms = parse_pdbqt_atoms(receptor_file)
    poses = parse_ligand_poses(ligand_file)
    if not poses:
        return
    ligand_atoms = poses[0]
    
    # Find pocket and interactions
    lig_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    pocket_atoms = []
    for ratom in receptor_atoms:
        rcoord = np.array([ratom['x'], ratom['y'], ratom['z']])
        dists = np.linalg.norm(lig_coords - rcoord, axis=1)
        if np.min(dists) <= 5.5:
            pocket_atoms.append(ratom)
    
    interactions = find_interactions(ligand_atoms, pocket_atoms)
    
    # Create RDKit mol from SMILES for 2D drawing
    smiles = get_ligand_smiles(ligand_name)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Failed to parse SMILES for {ligand_name}")
        return
    
    # Generate 2D coordinates
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    
    # Get 2D coordinates
    conf = mol.GetConformer()
    atom_2d = {}
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        atom_2d[i] = (pos.x, pos.y)
    
    # Draw molecule with RDKit
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 400)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.backgroundColour = (1, 1, 1)
    opts.bondLineWidth = 2
    opts.atomLabelFontSize = 14
    opts.bondColor = (0.15, 0.15, 0.15)
    
    # Color atoms
    atom_colors = {}
    highlight_atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        elem = atom.GetSymbol()
        if elem == 'O':
            atom_colors[i] = (0.9, 0.1, 0.1)
        elif elem == 'N':
            atom_colors[i] = (0.15, 0.35, 0.9)
        elif elem == 'S':
            atom_colors[i] = (0.8, 0.7, 0.2)
    
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightAtomColors=atom_colors)
    drawer.FinishDrawing()
    
    # Get PNG bytes and convert to PIL Image
    png = drawer.GetDrawingText()
    mol_img = Image.open(io.BytesIO(png))
    
    # Create main canvas
    img = Image.new('RGB', img_size, (255, 255, 255))
    draw = ImageDraw.Draw(img)
    
    # Place molecule in center
    mol_w, mol_h = mol_img.size
    mol_x = (img_size[0] - mol_w) // 2
    mol_y = (img_size[1] - mol_h) // 2 + 20
    img.paste(mol_img, (mol_x, mol_y))
    
    # Map molecule atoms to canvas coordinates
    # RDKit draws in a specific range, we need to map to actual pixels
    # Estimate: RDKit canvas is approximately 500x400 with some padding
    rdkit_scale_x = mol_w * 0.75  # approximate effective drawing area
    rdkit_scale_y = mol_h * 0.75
    rdkit_offset_x = mol_x + mol_w * 0.125
    rdkit_offset_y = mol_y + mol_h * 0.125
    
    # Get min/max of 2D coords for normalization
    if atom_2d:
        xs = [p[0] for p in atom_2d.values()]
        ys = [p[1] for p in atom_2d.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        range_x = max(max_x - min_x, 1)
        range_y = max(max_y - min_y, 1)
    else:
        return
    
    def map_2d_to_canvas(ax, ay):
        cx = rdkit_offset_x + (ax - min_x) / range_x * rdkit_scale_x
        cy = rdkit_offset_y + (ay - min_y) / range_y * rdkit_scale_y
        # RDKit y-axis is flipped relative to PIL
        cy = mol_y + mol_h - (cy - mol_y)
        return (cx, cy)
    
    # Place interaction labels around the molecule
    # Use polar placement based on angle from center
    center_x = img_size[0] // 2
    center_y = img_size[1] // 2 + 20
    
    n_interactions = len(interactions)
    if n_interactions == 0:
        # Draw title and save
        draw.text((20, 15), f"{receptor_name} x {ligand_name}", fill=(20, 20, 20))
        draw.text((20, 35), f"Affinity: {affinity:.2f} kcal/mol", fill=(80, 80, 80))
        img.save(output_file)
        return
    
    # Place labels in a fan arrangement
    # Sort by angle to avoid overlap
    labeled_positions = []
    
    for i, inter in enumerate(interactions):
        label = f"{inter['resname']}{inter['resseq']}"
        
        # Determine position angle
        angle = -90 + (180 / max(n_interactions, 1)) * i
        if n_interactions > 6:
            angle = -100 + (200 / n_interactions) * i
        
        rad = math.radians(angle)
        radius = min(img_size[0], img_size[1]) * 0.38
        
        lx = center_x + radius * math.cos(rad)
        ly = center_y + radius * math.sin(rad)
        
        # Adjust if too close to other labels
        for _ in range(10):
            too_close = False
            for px, py in labeled_positions:
                if abs(lx - px) < 80 and abs(ly - py) < 30:
                    too_close = True
                    ly += 35
                    break
            if not too_close:
                break
        
        labeled_positions.append((lx, ly))
        
        # Draw label box
        text_w = len(label) * 10
        box_x1 = int(lx - text_w / 2) - 6
        box_y1 = int(ly) - 10
        box_x2 = int(lx + text_w / 2) + 6
        box_y2 = int(ly) + 10
        
        # Color code the box based on interaction type
        box_color = inter['color']
        # Lighten for background
        bg_color = tuple(min(255, c + 180) for c in box_color)
        
        draw_rounded_rect(draw, (box_x1, box_y1, box_x2, box_y2), bg_color, radius=8)
        draw.text((int(lx) - len(label)*4, int(ly) - 6), label, fill=(20, 20, 20))
        
        # Draw dashed line to closest ligand atom
        # Find closest atom in 2D
        latom = inter['latom']
        
        # Find which RDKit atom index corresponds to this ligand atom
        # Use element match as simple heuristic
        best_idx = None
        best_dist = float('inf')
        for idx, (ax, ay) in atom_2d.items():
            cx, cy = map_2d_to_canvas(ax, ay)
            # We don't have exact mapping, so connect to center of molecule
            d = abs(cx - center_x) + abs(cy - center_y)
            if d < best_dist:
                best_dist = d
                best_idx = idx
        
        if best_idx is not None:
            ax, ay = atom_2d[best_idx]
            cx, cy = map_2d_to_canvas(ax, ay)
            
            # Draw curved or straight dashed line
            draw_dashed_line(draw, (cx, cy), (lx, ly - 10), color=box_color)
    
    # Title
    draw.text((20, 15), f"{receptor_name} x {ligand_name}", fill=(20, 20, 20))
    draw.text((20, 35), f"Binding Affinity: {affinity:.2f} kcal/mol", fill=(80, 80, 80))
    
    # Legend
    legend_y = img_size[1] - 60
    legend_items = [
        ('H-bond', HBOND_COLOR),
        ('Hydrophobic', HYDROPHOBIC_COLOR),
        ('Pi-Pi', PIPI_COLOR),
    ]
    lx = 30
    for text, color in legend_items:
        draw_rounded_rect(draw, (lx, legend_y, lx + 80, legend_y + 20), 
                         tuple(min(255, c + 180) for c in color), radius=6)
        draw.text((lx + 10, legend_y + 4), text, fill=(20, 20, 20))
        lx += 100
    
    img.save(output_file, 'PNG')
    print(f"  Saved: {os.path.basename(output_file)}")


# === MAIN ===
import io

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

print("Generating 2D interaction diagrams...")
for receptor, ligand, affinity in combinations:
    rec_file = os.path.join(RECEPTORS_DIR, f"{receptor}.pdbqt")
    lig_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_docked.pdbqt")
    out_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_2d_interaction.png")
    
    if not os.path.exists(lig_file):
        print(f"  SKIP: {lig_file} not found")
        continue
    
    print(f"\nProcessing: {receptor} x {ligand} ({affinity:.2f})")
    draw_2d_interaction(ligand, receptor, affinity, lig_file, rec_file, out_file)

print(f"\nAll 2D diagrams saved to: {RESULTS_DIR}")
