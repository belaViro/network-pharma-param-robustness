"""
True 3D binding mode visualization using PIL (no PyMOL needed).
Renders spheres with Phong shading for realistic 3D appearance.
"""
import os
import re
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import json

RESULTS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\results'
RECEPTORS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\receptors'

# Colors: C, N, O, H, S
ATOM_COLORS = {
    'C': (60, 60, 60),
    'N': (48, 80, 248),
    'O': (255, 13, 13),
    'H': (220, 220, 220),
    'S': (255, 255, 48),
    'P': (255, 128, 0),
}
ATOM_RADII = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'H': 1.2, 'S': 1.8, 'P': 1.8}


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
                if atom['element'] not in ATOM_COLORS:
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
                if atom['element'] not in ATOM_COLORS:
                    atom['element'] = atom['name'][0].upper()
                current.append(atom)
        if current and not poses:
            poses.append(current)
    return poses


def draw_sphere(draw, cx, cy, r, color, light_dir=(0.3, -0.5, 1.0)):
    """Draw a 3D sphere with Phong-like shading."""
    r = max(r, 2)
    light_dir = np.array(light_dir)
    light_dir = light_dir / np.linalg.norm(light_dir)
    
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            d2 = dx * dx + dy * dy
            if d2 > r * r:
                continue
            # Normal on sphere surface
            dz = np.sqrt(max(0, r * r - d2))
            normal = np.array([dx, dy, dz])
            normal = normal / np.linalg.norm(normal)
            
            # Diffuse lighting
            diffuse = max(0.0, np.dot(normal, light_dir))
            ambient = 0.35
            intensity = ambient + (1 - ambient) * diffuse
            
            # Specular highlight
            view_dir = np.array([0, 0, 1])
            half = light_dir + view_dir
            half = half / np.linalg.norm(half)
            spec = max(0.0, np.dot(normal, half)) ** 20
            intensity += 0.3 * spec
            intensity = min(1.0, intensity)
            
            px = cx + dx
            py = cy + dy
            col = tuple(min(255, int(c * intensity)) for c in color)
            draw.point((px, py), fill=col)


def draw_line_3d(draw, p1, p2, color=(80, 80, 80), width=2):
    """Draw a thick line."""
    x1, y1 = int(p1[0]), int(p1[1])
    x2, y2 = int(p2[0]), int(p2[1])
    draw.line([(x1, y1), (x2, y2)], fill=color, width=width)


def draw_dashed_line(draw, p1, p2, color=(200, 180, 0), dash=6):
    """Draw a dashed line."""
    x1, y1 = int(p1[0]), int(p1[1])
    x2, y2 = int(p2[0]), int(p2[1])
    dist = np.hypot(x2 - x1, y2 - y1)
    if dist < 1:
        return
    n_dashes = int(dist / dash)
    for i in range(0, n_dashes, 2):
        t1 = i / n_dashes
        t2 = min((i + 1) / n_dashes, 1.0)
        x1d = int(x1 + (x2 - x1) * t1)
        y1d = int(y1 + (y2 - y1) * t1)
        x2d = int(x1 + (x2 - x1) * t2)
        y2d = int(y1 + (y2 - y1) * t2)
        draw.line([(x1d, y1d), (x2d, y2d)], fill=color, width=2)


def project(coords, center, scale, offset_x, offset_y):
    """Simple perspective projection."""
    proj = []
    for c in coords:
        x = (c[0] - center[0]) * scale + offset_x
        y = (c[1] - center[1]) * scale + offset_y
        # Simple depth = z
        z = c[2] - center[2]
        proj.append((x, y, z))
    return proj


def render_binding_mode(receptor_file, ligand_file, receptor_name, ligand_name,
                        affinity, output_file, img_size=(1200, 1000)):
    """Render a true 3D binding mode image."""
    
    receptor_atoms = parse_pdbqt_atoms(receptor_file)
    poses = parse_ligand_poses(ligand_file)
    if not poses:
        print(f"  No poses: {receptor_name}x{ligand_name}")
        return
    ligand_atoms = poses[0]
    
    # Find pocket atoms (within 5.5A of any ligand atom)
    lig_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    pocket_atoms = []
    for ratom in receptor_atoms:
        rcoord = np.array([ratom['x'], ratom['y'], ratom['z']])
        dists = np.linalg.norm(lig_coords - rcoord, axis=1)
        if np.min(dists) <= 5.5 and ratom['element'] != 'H':
            pocket_atoms.append(ratom)
    
    # Center on ligand
    lig_center = lig_coords.mean(axis=0)
    
    # Combine all atoms for rendering
    all_render = []
    
    # Pocket atoms (small, gray, low z-priority)
    for a in pocket_atoms:
        r = ATOM_RADII.get(a['element'], 1.5) * 5
        col = (160, 160, 160)
        all_render.append({
            'x': a['x'], 'y': a['y'], 'z': a['z'],
            'r': r, 'color': col, 'type': 'pocket', 'atom': a
        })
    
    # Ligand bonds (we'll draw after atoms)
    ligand_bonds = []
    for i, a1 in enumerate(ligand_atoms):
        for j, a2 in enumerate(ligand_atoms[i+1:], i+1):
            d = np.linalg.norm(
                np.array([a1['x'], a1['y'], a1['z']]) -
                np.array([a2['x'], a2['y'], a2['z']])
            )
            max_bond = 1.8 if a1['element'] != 'H' and a2['element'] != 'H' else 1.3
            if d <= max_bond:
                ligand_bonds.append((a1, a2))
    
    # Find H-bonds
    hbonds = []
    acceptors = ['O', 'N']
    for latom in ligand_atoms:
        for ratom in pocket_atoms:
            if ratom['element'] == 'H':
                continue
            d = np.linalg.norm(
                np.array([latom['x'], latom['y'], latom['z']]) -
                np.array([ratom['x'], ratom['y'], ratom['z']])
            )
            if d <= 3.2 and latom['element'] in acceptors and ratom['element'] in acceptors:
                hbonds.append((latom, ratom))
    
    # Scale and project
    # Auto-scale to fit image
    all_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms + pocket_atoms])
    if len(all_coords) == 0:
        return
    spread = np.ptp(all_coords, axis=0).max()
    scale = min(img_size[0], img_size[1]) * 0.35 / max(spread, 1)
    offset_x = img_size[0] // 2
    offset_y = img_size[1] // 2 + 30
    
    # Project ligand atoms
    lig_proj = project(lig_coords, lig_center, scale, offset_x, offset_y)
    
    # Add ligand atoms to render list (larger, colored)
    for i, a in enumerate(ligand_atoms):
        r = ATOM_RADII.get(a['element'], 1.5) * 9
        col = ATOM_COLORS.get(a['element'], (128, 128, 128))
        all_render.append({
            'x': a['x'], 'y': a['y'], 'z': a['z'],
            'r': r, 'color': col, 'type': 'ligand', 'atom': a,
            'px': lig_proj[i][0], 'py': lig_proj[i][1], 'pz': lig_proj[i][2]
        })
    
    # Project pocket atoms
    pocket_coords = np.array([[a['x'], a['y'], a['z']] for a in pocket_atoms])
    if len(pocket_coords) > 0:
        pocket_proj = project(pocket_coords, lig_center, scale, offset_x, offset_y)
        for i, a in enumerate(pocket_atoms):
            # Update existing pocket entries
            for item in all_render:
                if item['type'] == 'pocket' and item['atom'] is a:
                    item['px'] = pocket_proj[i][0]
                    item['py'] = pocket_proj[i][1]
                    item['pz'] = pocket_proj[i][2]
                    break
    
    # Sort by z (painter's algorithm: back to front)
    all_render = [r for r in all_render if 'px' in r]
    all_render.sort(key=lambda x: x['pz'])
    
    # Create image
    img = Image.new('RGB', img_size, (255, 255, 255))
    draw = ImageDraw.Draw(img)
    
    # Draw pocket bonds (CA trace)
    ca_atoms = [(a, p) for a, p in zip(pocket_atoms, pocket_proj) if a['name'] == 'CA']
    ca_atoms.sort(key=lambda x: x[0]['resseq'])
    for i in range(len(ca_atoms) - 1):
        a1, p1 = ca_atoms[i]
        a2, p2 = ca_atoms[i + 1]
        # Only connect if close in space
        d = np.linalg.norm(np.array([a1['x'], a1['y'], a1['z']]) - np.array([a2['x'], a2['y'], a2['z']]))
        if d < 5.0:
            draw.line([(int(p1[0]), int(p1[1])), (int(p2[0]), int(p2[1]))],
                     fill=(180, 180, 180), width=2)
    
    # Draw H-bonds
    for latom, ratom in hbonds:
        l_proj = project(np.array([[latom['x'], latom['y'], latom['z']]]), lig_center, scale, offset_x, offset_y)[0]
        r_proj = project(np.array([[ratom['x'], ratom['y'], ratom['z']]]), lig_center, scale, offset_x, offset_y)[0]
        draw_dashed_line(draw, (l_proj[0], l_proj[1]), (r_proj[0], r_proj[1]))
    
    # Draw ligand bonds
    for a1, a2 in ligand_bonds:
        p1 = project(np.array([[a1['x'], a1['y'], a1['z']]]), lig_center, scale, offset_x, offset_y)[0]
        p2 = project(np.array([[a2['x'], a2['y'], a2['z']]]), lig_center, scale, offset_x, offset_y)[0]
        draw.line([(int(p1[0]), int(p1[1])), (int(p2[0]), int(p2[1]))],
                 fill=(40, 40, 40), width=3)
    
    # Draw atoms (sorted by depth)
    for item in all_render:
        draw_sphere(draw, int(item['px']), int(item['py']), int(item['r']),
                   item['color'])
    
    # Draw labels for key pocket residues
    # Pick residues closest to ligand
    pocket_residues = {}
    for a, p in zip(pocket_atoms, pocket_proj):
        key = (a['resname'], a['resseq'])
        if key not in pocket_residues:
            pocket_residues[key] = (a, p, float('inf'))
        dist_to_lig = np.min(np.linalg.norm(lig_coords - np.array([a['x'], a['y'], a['z']]), axis=1))
        if dist_to_lig < pocket_residues[key][2]:
            pocket_residues[key] = (a, p, dist_to_lig)
    
    # Sort by distance, label top 6
    sorted_res = sorted(pocket_residues.items(), key=lambda x: x[1][2])
    labeled = 0
    for (resname, resseq), (a, p, dist) in sorted_res:
        if dist > 6.0 or labeled >= 6:
            continue
        if resname in ['', 'UNK', 'UNL']:
            continue
        label = f"{resname}{resseq}"
        tx, ty = int(p[0]) + 15, int(p[1]) - 10
        # Draw text with outline
        for dx, dy in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            draw.text((tx + dx, ty + dy), label, fill=(255, 255, 255), font=None)
        draw.text((tx, ty), label, fill=(50, 50, 50), font=None)
        labeled += 1
    
    # Title
    title = f"{receptor_name} x {ligand_name}    Affinity: {affinity:.2f} kcal/mol"
    draw.text((20, 15), title, fill=(20, 20, 20), font=None)
    
    # Info
    info = f"H-bonds: {len(hbonds)}  |  Pocket residues: {len(pocket_residues)}"
    draw.text((20, img_size[1] - 30), info, fill=(80, 80, 80), font=None)
    
    # Legend
    legend_y = img_size[1] - 80
    legend_items = [
        ('C', ATOM_COLORS['C']),
        ('N', ATOM_COLORS['N']),
        ('O', ATOM_COLORS['O']),
        ('H', ATOM_COLORS['H']),
    ]
    lx = img_size[0] - 120
    for elem, col in legend_items:
        draw_sphere(draw, lx + 10, legend_y, 8, col)
        draw.text((lx + 25, legend_y - 6), elem, fill=(50, 50, 50), font=None)
        legend_y += 20
    
    img.save(output_file, 'PNG')
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

print("Rendering true 3D binding mode images with Phong shading...")
for receptor, ligand, affinity in combinations:
    rec_file = os.path.join(RECEPTORS_DIR, f"{receptor}.pdbqt")
    lig_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_docked.pdbqt")
    out_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_3d.png")
    
    if not os.path.exists(lig_file):
        print(f"  SKIP: {lig_file} not found")
        continue
    
    print(f"\nProcessing: {receptor} x {ligand} ({affinity:.2f})")
    render_binding_mode(rec_file, lig_file, receptor, ligand, affinity, out_file)

print(f"\nAll 3D images saved to: {RESULTS_DIR}")
