"""
Generate publication-quality 3D protein-ligand binding mode images using VTK.
Renders protein cartoon (ribbon) + ligand stick model + key residue labels.
"""
import os
import re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

RESULTS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\results'
RECEPTORS_DIR = r'D:\科研\网络药理学参数敏感性分析\分子对接\receptors'


def parse_pdbqt_atoms(filepath):
    """Parse PDBQT, return atom list."""
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
                    'occupancy': float(line[54:60].strip()) if len(line) > 59 and line[54:60].strip() else 1.0,
                    'temp_factor': float(line[60:66].strip()) if len(line) > 65 and line[60:66].strip() else 0.0,
                }
                if atom['element'] not in ['C', 'N', 'O', 'H', 'S', 'P']:
                    atom['element'] = atom['name'][0].upper()
                atoms.append(atom)
    return atoms


def parse_ligand_poses(filepath):
    """Parse multi-model ligand PDBQT."""
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


def write_temp_pdb(atoms, filepath):
    """Write atoms to a temporary PDB file for VTK reader."""
    with open(filepath, 'w') as f:
        for i, atom in enumerate(atoms):
            f.write(f"ATOM  {i+1:5d} {atom['name']:4s} {atom['resname']:3s} {atom.get('chain', 'A')}"
                   f"{atom.get('resseq', i+1):4d}    "
                   f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                   f"{atom.get('occupancy', 1.0):6.2f}{atom.get('temp_factor', 0.0):6.2f}"
                   f"          {atom['element']:>2s}\n")
        f.write("END\n")


def create_tube_from_points(points, radius=0.3, color=(0.6, 0.7, 0.8)):
    """Create a tube from a list of points (for CA trace)."""
    vtk_points = vtk.vtkPoints()
    for p in points:
        vtk_points.InsertNextPoint(p[0], p[1], p[2])
    
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(len(points))
    for i in range(len(points)):
        polyline.GetPointIds().SetId(i, i)
    
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)
    
    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(vtk_points)
    poly_data.SetLines(cells)
    
    tube_filter = vtk.vtkTubeFilter()
    tube_filter.SetInputData(poly_data)
    tube_filter.SetRadius(radius)
    tube_filter.SetNumberOfSides(12)
    tube_filter.Update()
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube_filter.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(0.85)
    
    return actor


def create_sphere(center, radius, color=(1, 0, 0)):
    """Create a sphere actor."""
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(center[0], center[1], center[2])
    sphere.SetRadius(radius)
    sphere.SetPhiResolution(16)
    sphere.SetThetaResolution(16)
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    
    return actor


def create_cylinder(p1, p2, radius, color=(0.3, 0.3, 0.3)):
    """Create a cylinder (bond) between two points."""
    p1_arr = np.array(p1)
    p2_arr = np.array(p2)
    center = (p1_arr + p2_arr) / 2
    length = np.linalg.norm(p2_arr - p1_arr)
    
    if length < 0.01:
        return None
    
    direction = (p2_arr - p1_arr) / length
    
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(radius)
    cylinder.SetHeight(length)
    cylinder.SetResolution(12)
    
    # Transform to align with bond direction
    transform = vtk.vtkTransform()
    transform.Translate(center)
    
    # Rotate cylinder (default is along Y axis)
    # Calculate rotation axis and angle
    y_axis = np.array([0, 1, 0])
    rot_axis = np.cross(y_axis, direction)
    rot_angle = np.degrees(np.arccos(np.dot(y_axis, direction)))
    
    if np.linalg.norm(rot_axis) > 0.001:
        transform.RotateWXYZ(rot_angle, rot_axis[0], rot_axis[1], rot_axis[2])
    elif direction[1] < 0:
        transform.RotateWXYZ(180, 1, 0, 0)
    
    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputConnection(cylinder.GetOutputPort())
    transform_filter.SetTransform(transform)
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(transform_filter.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    
    return actor


def create_label(text, position, color=(0, 0, 0), font_size=16):
    """Create a text label at a 3D position."""
    label = vtk.vtkBillboardTextActor3D()
    label.SetInput(text)
    label.SetPosition(position)
    label.GetTextProperty().SetFontSize(font_size)
    label.GetTextProperty().SetColor(*color)
    label.GetTextProperty().SetBold(True)
    label.GetTextProperty().SetBackgroundColor(1, 1, 1)
    label.GetTextProperty().SetBackgroundOpacity(0.8)
    
    return label


def render_binding_mode(receptor_file, ligand_file, receptor_name, ligand_name,
                        affinity, output_file, img_size=(1200, 1000)):
    """Render a 3D binding mode image with VTK."""
    
    # Parse structures
    receptor_atoms = parse_pdbqt_atoms(receptor_file)
    poses = parse_ligand_poses(ligand_file)
    if not poses:
        print(f"  No poses found")
        return
    ligand_atoms = poses[0]
    
    # Find pocket residues
    lig_coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    lig_center = lig_coords.mean(axis=0)
    
    pocket_atoms = []
    for ratom in receptor_atoms:
        rcoord = np.array([ratom['x'], ratom['y'], ratom['z']])
        dists = np.linalg.norm(lig_coords - rcoord, axis=1)
        if np.min(dists) <= 6.0:
            pocket_atoms.append(ratom)
    
    # --- Create renderer ---
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    
    # --- Render protein pocket as CA trace tubes ---
    ca_atoms = [a for a in pocket_atoms if a['name'] == 'CA']
    ca_atoms.sort(key=lambda a: a['resseq'])
    
    # Split into segments (don't connect across large gaps)
    segments = []
    current_seg = []
    for i, ca in enumerate(ca_atoms):
        if i == 0:
            current_seg.append(ca)
        else:
            prev = ca_atoms[i-1]
            dist = np.linalg.norm(
                np.array([ca['x'], ca['y'], ca['z']]) - 
                np.array([prev['x'], prev['y'], prev['z']])
            )
            if dist < 5.0:
                current_seg.append(ca)
            else:
                if len(current_seg) >= 2:
                    segments.append(current_seg)
                current_seg = [ca]
    if len(current_seg) >= 2:
        segments.append(current_seg)
    
    # Color by secondary structure approximation
    for seg in segments:
        points = [(a['x'], a['y'], a['z']) for a in seg]
        # Determine color based on segment length (rough approximation)
        if len(seg) >= 8:
            color = (0.5, 0.35, 0.6)  # Purple - alpha helix-ish
        elif len(seg) >= 4:
            color = (0.8, 0.7, 0.3)   # Yellow - beta sheet-ish
        else:
            color = (0.4, 0.6, 0.6)   # Teal - coil
        
        actor = create_tube_from_points(points, radius=0.35, color=color)
        renderer.AddActor(actor)
    
    # --- Render ligand as spheres and cylinders (stick model) ---
    elem_colors = {
        'C': (0.2, 0.2, 0.2),
        'N': (0.2, 0.3, 0.8),
        'O': (0.9, 0.1, 0.1),
        'H': (0.9, 0.9, 0.9),
        'S': (0.9, 0.8, 0.2),
        'P': (0.9, 0.5, 0.1),
    }
    elem_radii = {'C': 0.25, 'N': 0.23, 'O': 0.22, 'H': 0.15, 'S': 0.28, 'P': 0.30}
    bond_max = {'C-C': 1.8, 'C-N': 1.6, 'C-O': 1.6, 'C-H': 1.2, 'N-H': 1.2, 'O-H': 1.2, 'default': 1.8}
    
    # Draw bonds
    for i, a1 in enumerate(ligand_atoms):
        for j, a2 in enumerate(ligand_atoms[i+1:], i+1):
            p1 = np.array([a1['x'], a1['y'], a1['z']])
            p2 = np.array([a2['x'], a2['y'], a2['z']])
            d = np.linalg.norm(p2 - p1)
            
            e1, e2 = a1['element'], a2['element']
            max_d = bond_max.get(f"{e1}-{e2}", bond_max.get(f"{e2}-{e1}", bond_max['default']))
            if e1 == 'H' or e2 == 'H':
                max_d = 1.3
            
            if d <= max_d:
                color = elem_colors.get(e1, (0.5, 0.5, 0.5))
                actor = create_cylinder(p1, p2, 0.12, color)
                if actor:
                    renderer.AddActor(actor)
    
    # Draw atoms
    for atom in ligand_atoms:
        color = elem_colors.get(atom['element'], (0.5, 0.5, 0.5))
        radius = elem_radii.get(atom['element'], 0.2)
        actor = create_sphere((atom['x'], atom['y'], atom['z']), radius, color)
        renderer.AddActor(actor)
    
    # --- Render key pocket residues as spheres ---
    pocket_residues = {}
    for a in pocket_atoms:
        key = (a['resname'], a['resseq'])
        if key not in pocket_residues:
            pocket_residues[key] = []
        pocket_residues[key].append(a)
    
    # Find closest residues to ligand
    residue_distances = []
    for key, atoms in pocket_residues.items():
        resname, resseq = key
        if resname in ['', 'UNK', 'UNL']:
            continue
        # Use CA or first atom as representative
        rep = next((a for a in atoms if a['name'] == 'CA'), atoms[0])
        rcoord = np.array([rep['x'], rep['y'], rep['z']])
        min_dist = np.min(np.linalg.norm(lig_coords - rcoord, axis=1))
        residue_distances.append((key, rep, min_dist))
    
    residue_distances.sort(key=lambda x: x[2])
    
    # Show top 6 residues
    for (resname, resseq), rep, dist in residue_distances[:6]:
        # Draw residue side chain as thin tube
        res_atoms = pocket_residues[(resname, resseq)]
        if len(res_atoms) >= 2:
            points = [(a['x'], a['y'], a['z']) for a in res_atoms]
            actor = create_tube_from_points(points, radius=0.15, color=(0.3, 0.5, 0.5))
            renderer.AddActor(actor)
        
        # Draw CA as slightly larger sphere
        ca = next((a for a in res_atoms if a['name'] == 'CA'), res_atoms[0])
        actor = create_sphere((ca['x'], ca['y'], ca['z']), 0.4, (0.3, 0.5, 0.5))
        renderer.AddActor(actor)
        
        # Label
        label_text = f"{resname}{resseq}"
        label = create_label(label_text, (ca['x'], ca['y'], ca['z'] + 1.5), 
                            color=(0.1, 0.1, 0.1), font_size=14)
        renderer.AddActor(label)
    
    # --- Setup camera ---
    camera = renderer.GetActiveCamera()
    camera.SetFocalPoint(lig_center[0], lig_center[1], lig_center[2])
    camera.SetPosition(lig_center[0] + 15, lig_center[1] + 15, lig_center[2] + 20)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    
    # Zoom in
    renderer.ResetCameraClippingRange()
    
    # --- Setup lighting ---
    light1 = vtk.vtkLight()
    light1.SetPosition(lig_center[0] + 20, lig_center[1] + 20, lig_center[2] + 30)
    light1.SetFocalPoint(lig_center[0], lig_center[1], lig_center[2])
    light1.SetColor(1, 1, 1)
    light1.SetIntensity(1.0)
    renderer.AddLight(light1)
    
    light2 = vtk.vtkLight()
    light2.SetPosition(lig_center[0] - 20, lig_center[1] - 10, lig_center[2] + 10)
    light2.SetFocalPoint(lig_center[0], lig_center[1], lig_center[2])
    light2.SetColor(0.8, 0.8, 0.9)
    light2.SetIntensity(0.5)
    renderer.AddLight(light2)
    
    # --- Render offscreen ---
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(img_size[0], img_size[1])
    render_window.SetOffScreenRendering(1)
    render_window.Render()
    
    # Capture to image
    window_to_image = vtk.vtkWindowToImageFilter()
    window_to_image.SetInput(render_window)
    window_to_image.SetInputBufferTypeToRGB()
    window_to_image.ReadFrontBufferOff()
    window_to_image.Update()
    
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output_file)
    writer.SetInputConnection(window_to_image.GetOutputPort())
    writer.Write()
    
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

print("Generating VTK 3D cartoon images...")
for receptor, ligand, affinity in combinations:
    rec_file = os.path.join(RECEPTORS_DIR, f"{receptor}.pdbqt")
    lig_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_docked.pdbqt")
    out_file = os.path.join(RESULTS_DIR, f"{receptor}_{ligand}_cartoon.png")
    
    if not os.path.exists(lig_file):
        print(f"  SKIP: {lig_file} not found")
        continue
    
    print(f"\nProcessing: {receptor} x {ligand} ({affinity:.2f})")
    try:
        render_binding_mode(rec_file, lig_file, receptor, ligand, affinity, out_file)
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

print(f"\nAll VTK images saved to: {RESULTS_DIR}")
