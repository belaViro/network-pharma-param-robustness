# PyMOL Script: EGFR x Catechin
# Affinity: -8.27 kcal/mol
# Run in PyMOL: File -> Run Script -> select this .pml file

# --- Load structures ---
cmd.load("D:/科研/网络药理学参数敏感性分析/分子对接/receptors/EGFR.pdbqt", "protein")
cmd.load("D:/科研/网络药理学参数敏感性分析/分子对接/results/EGFR_Catechin_docked.pdbqt", "ligand")

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
cmd.png("D:/科研/网络药理学参数敏感性分析/分子对接/results/EGFR_Catechin_pymol.png", dpi=300)

print("Image saved: D:/科研/网络药理学参数敏感性分析/分子对接/results/EGFR_Catechin_pymol.png")
