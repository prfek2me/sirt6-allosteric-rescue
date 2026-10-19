import os
import re
import sys
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

holo_dir = "results/holo_mdl800"
os.makedirs(holo_dir, exist_ok=True)

# 1. Rescatar la pose ganadora Top-1 de MDL800 en c1 y restituir hidrógenos 3D completos
pose_sdf = "results/docking/poses/MDL800_c1.sdf"
mdl_smiles = "COC(=O)c1cc(NS(=O)(=O)c2cc(Cl)cc(Cl)c2)ccc1S(=O)(=O)Nc1cc(Br)c(F)cc1C"

template = Chem.MolFromSmiles(mdl_smiles)
docked_raw = Chem.SDMolSupplier(pose_sdf, removeHs=True)[0]
mol_clean = AllChem.AssignBondOrdersFromTemplate(template, docked_raw)
mol_3d_h = Chem.AddHs(mol_clean, addCoords=True)
mol_3d_h.SetProp("_Name", "MDL800")

lig_sdf_full = os.path.join(holo_dir, "MDL800_c1_fullH.sdf")
with Chem.SDWriter(lig_sdf_full) as w:
    w.write(mol_3d_h)
print(f"[✓] Pose Top-1 de MDL800_c1 reconstruida con valencia 3D completa ({mol_3d_h.GetNumAtoms()} átomos).")

# 2. Parametrizar MDL800 con ANTECHAMBER (GAFF2 + AM1-BCC) y parmchk2
print("[*] Calculando cargas semiempíricas AM1-BCC con ANTECHAMBER/sqm (2-3 min)...")
cmd_ante = [
    "antechamber",
    "-i", "MDL800_c1_fullH.sdf", "-fi", "sdf",
    "-o", "MDL800.mol2", "-fo", "mol2",
    "-c", "bcc", "-nc", "0", "-m", "1",
    "-at", "gaff2", "-rn", "MDL", "-pf", "y"
]
res_ante = subprocess.run(cmd_ante, cwd=holo_dir, capture_output=True, text=True)
if res_ante.returncode != 0:
    print("[!] Error en antechamber:\n", res_ante.stderr, res_ante.stdout)
    sys.exit(1)

cmd_parm = ["parmchk2", "-i", "MDL800.mol2", "-f", "mol2", "-o", "MDL800.frcmod", "-s", "2"]
res_parm = subprocess.run(cmd_parm, cwd=holo_dir, capture_output=True, text=True)
if res_parm.returncode != 0:
    print("[!] Error en parmchk2:\n", res_parm.stderr)
    sys.exit(1)
print("[✓] Parámetros GAFF2 (MDL800.mol2 y MDL800.frcmod) generados con éxito.")

# 3. Transferir coordenadas 3D del centroide c1 al molde topológico exacto de MCPB.py
c1_pdb = "results/clustering/centroid.c1.pdb"
mcpb_pdb = "results/mcpb/SIRT6_Zn_mcpbpy.pdb"
holo_prot_pdb = os.path.join(holo_dir, "SIRT6_c1_mcpbpy.pdb")

# Mapear coordenadas de centroid.c1.pdb por (índice secuencial de residuo 1..283, nombre de átomo)
c1_coords = {}
seq_res = 0
last_res_id = None
with open(c1_pdb, "r") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            res_uid = line[21:27]
            if res_uid != last_res_id:
                seq_res += 1
                last_res_id = res_uid
            aname = line[12:16].strip()
            xyz_str = line[30:54]
            c1_coords[(seq_res, aname)] = xyz_str

seq_res_mcpb = 0
last_res_mcpb = None
updated_atoms = 0
with open(mcpb_pdb, "r") as fin, open(holo_prot_pdb, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            res_uid = line[21:27]
            if res_uid != last_res_id:
                pass
            if res_uid != last_res_mcpb:
                seq_res_mcpb += 1
                last_res_mcpb = res_uid
            aname = line[12:16].strip()
            if (seq_res_mcpb, aname) in c1_coords:
                line = line[:30] + c1_coords[(seq_res_mcpb, aname)] + line[54:]
                updated_atoms += 1
        fout.write(line)
print(f"[✓] Molde cuántico actualizado a las coordenadas del centroide c1 ({updated_atoms} átomos anclados).")

# 4. Copiar artefactos cuánticos de MCPB.py y construir el enrutador tleap Holo
for fname in ["CM1.mol2", "CM2.mol2", "CM3.mol2", "CM4.mol2", "ZN1.mol2", "SIRT6_Zn_mcpbpy.frcmod"]:
    shutil.copy(os.path.join("results/mcpb", fname), os.path.join(holo_dir, fname))

with open("results/mcpb/SIRT6_Zn_tleap.in", "r") as f:
    apo_tleap = f.read()

# Extraer el encabezado cuántico y los enlaces de coordinación exactos hasta antes de saveamberparm
header_and_bonds = []
for line in apo_tleap.splitlines():
    if "saveamberparm" in line or "solvatebox" in line or "addions" in line or line.strip() == "quit":
        break
    if "loadpdb" in line:
        line = "mol = loadpdb SIRT6_c1_mcpbpy.pdb"
    header_and_bonds.append(line)

holo_tleap_str = "\n".join(header_and_bonds) + """
source leaprc.gaff2
loadamberparams MDL800.frcmod
MDL = loadmol2 MDL800.mol2
holo = combine { mol MDL }

savepdb holo SIRT6_Holo_MDL800_dry.pdb
saveamberparm holo SIRT6_Holo_MDL800_dry.prmtop SIRT6_Holo_MDL800_dry.inpcrd

solvatebox holo OPCBOX 10.0
addions holo Na+ 0
addions holo Cl- 0

savepdb holo SIRT6_Holo_MDL800_solv.pdb
saveamberparm holo SIRT6_Holo_MDL800_solv.prmtop SIRT6_Holo_MDL800_solv.inpcrd
quit
"""

tleap_in_path = os.path.join(holo_dir, "SIRT6_Holo_tleap.in")
with open(tleap_in_path, "w") as f:
    f.write(holo_tleap_str)

print("[*] Compilando topología ternaria (SIRT6_c1 + Zn + MDL800) en tleap...")
res_tleap = subprocess.run(["tleap", "-s", "-f", "SIRT6_Holo_tleap.in"], cwd=holo_dir, capture_output=True, text=True)
with open(os.path.join(holo_dir, "tleap_holo.log"), "w") as f:
    f.write(res_tleap.stdout)

solv_prmtop = os.path.join(holo_dir, "SIRT6_Holo_MDL800_solv.prmtop")
if not os.path.exists(solv_prmtop) or os.path.getsize(solv_prmtop) == 0:
    print("[!] Error en tleap. Revisar results/holo_mdl800/tleap_holo.log")
    sys.exit(1)

print(f"[✓] Topología Holo-MD compilada exitosamente: {solv_prmtop}")
