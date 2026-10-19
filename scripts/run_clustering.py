import os
import sys
import subprocess
import numpy as np

os.makedirs("results/clustering", exist_ok=True)

raw_pdb = "data/raw/5MF6.pdb"
dry_pdb = "results/mcpb/SIRT6_Zn_dry.pdb"

excluded_het = {"HOH", "WAT", "SO4", "EDO", "PGE", "ZN", "ZN1", "AR6", "GOL", "CL", "NA"}
het_atoms_by_res = {}

with open(raw_pdb, "r") as f:
    for line in f:
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            chain = line[21].strip()
            if chain == "A" and resname not in excluded_het:
                resseq = line[22:26].strip()
                key = f"{resname}_{resseq}"
                coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                het_atoms_by_res.setdefault(key, []).append(coords)

target_lig_key = max(het_atoms_by_res, key=lambda k: len(het_atoms_by_res[k]))
lig_coords = np.array(het_atoms_by_res[target_lig_key])
print(f"[*] Ligando alostérico de referencia en 5MF6: {target_lig_key} ({len(lig_coords)} átomos)")

pocket_resids = set()
with open(dry_pdb, "r") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            resid = int(line[22:26].strip())
            if resid <= 282:
                atom_coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if np.min(np.linalg.norm(lig_coords - atom_coord, axis=1)) <= 8.0:
                    pocket_resids.add(resid)

sorted_resids = sorted(pocket_resids)
pocket_mask = ":" + ",".join(map(str, sorted_resids))

pocket_center = np.mean(lig_coords, axis=0)
with open("results/clustering/pocket_box.txt", "w") as f:
    f.write(f"center_x = {pocket_center[0]:.3f}\n")
    f.write(f"center_y = {pocket_center[1]:.3f}\n")
    f.write(f"center_z = {pocket_center[2]:.3f}\n")
    f.write(f"pocket_mask = {pocket_mask}\n")

# Cargamos la topología seca como referencia [crystal] para anclar el origen cartesiano al cristal 5MF6
cpptraj_in = f"""parm results/mcpb/SIRT6_Zn_solv.prmtop [solv]
parm results/mcpb/SIRT6_Zn_dry.prmtop [dry]
reference results/mcpb/SIRT6_Zn_dry.inpcrd parm [dry] [crystal]

trajin results/mcpb/SIRT6_Zn_prod_50ns.dcd 201 1000 parm [solv]
autoimage anchor :1-283
strip :WAT,Na+,Cl-

# 1. Auditar tetraedro cuántico Zn-S4
rms ZN_TETRAHEDRON :CM1,CM2,CM3,CM4,ZN1&!@H= first out results/clustering/rmsd_clean.dat

# 2. Medir RMSD del esqueleto global contra el cristal
rms BB_RMSD :1-282@CA,C,N ref [crystal] :1-282@CA,C,N out results/clustering/rmsd_clean.dat

# 3. Anclar las coordenadas del bolsillo alostérico exactamente sobre el marco cristalográfico [crystal]
rms POCKET_RMSD {pocket_mask}@CA,C,N ref [crystal] {pocket_mask}@CA,C,N out results/clustering/rmsd_clean.dat

# 4. K-means (k=5) sobre todos los átomos pesados del bolsillo ya co-centrado con el cristal
cluster C_POCKET \\
    kmeans clusters 5 randompoint maxit 500 kseed 42 \\
    rms {pocket_mask}&!@H= nofit \\
    out results/clustering/cnumvtime.dat \\
    summary results/clustering/cluster_summary.dat \\
    info results/clustering/cluster_info.dat \\
    repout results/clustering/centroid repfmt pdb

go
quit
"""

in_file = "results/clustering/cluster_pocket.in"
with open(in_file, "w") as f:
    f.write(cpptraj_in)

print("[*] Re-alineando trayectoria contra referencia cristalográfica y extrayendo centroides...")
res = subprocess.run(["cpptraj", "-i", in_file], capture_output=True, text=True)
if res.returncode != 0:
    print("[!] Error en cpptraj:\n", res.stderr)
    sys.exit(1)

print("[✓] Centroides anclados al marco cristalográfico de 5MF6.")
