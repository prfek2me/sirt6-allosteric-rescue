import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

os.makedirs("results/docking/receptors", exist_ok=True)
os.makedirs("results/docking/ligands", exist_ok=True)

# 1. Extraer ligando cristalográfico 7M2_410 (UBCS039) de 5MF6.pdb (Cadena A)
raw_pdb = "data/raw/5MF6.pdb"
crystal_block = []
with open(raw_pdb, "r") as f:
    for line in f:
        if line.startswith("HETATM") and line[17:20].strip() == "7M2" and line[21].strip() == "A":
            crystal_block.append(line)

if not crystal_block:
    print("[!] Error: No se encontró 7M2 en la cadena A de 5MF6.pdb")
    sys.exit(1)

crystal_pdb_str = "".join(crystal_block) + "END\n"
with open("results/docking/ligands/UBCS039_crystal_ref.pdb", "w") as f:
    f.write(crystal_pdb_str)

ubcs_smiles = "C1=CC=C2C(=C1)NC(C3=CC=CN32)C4=CN=CC=C4"
template = Chem.MolFromSmiles(ubcs_smiles)
crystal_mol = Chem.MolFromPDBBlock(crystal_pdb_str, removeHs=True)
crystal_mol = AllChem.AssignBondOrdersFromTemplate(template, crystal_mol)
Chem.AssignAtomChiralTagsFromStructure(crystal_mol)
Chem.AssignStereochemistry(crystal_mol, cleanIt=True, force=True)

with Chem.SDWriter("results/docking/ligands/UBCS039_crystal_ref.sdf") as w:
    w.write(crystal_mol)

ubcs_isomeric_smiles = Chem.MolToSmiles(crystal_mol, isomericSmiles=True)

# 2. Generar ligandos 3D minimizados (MMFF94s)
ligands_smiles = {
    "UBCS039": ubcs_isomeric_smiles,
    "MDL800": "COC(=O)c1cc(NS(=O)(=O)c2cc(Cl)cc(Cl)c2)ccc1S(=O)(=O)Nc1cc(Br)c(F)cc1C",
    "NAM": "NC(=O)c1cccnc1"
}

for name, smi in ligands_smiles.items():
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=2000)
    mol.SetProp("_Name", name)
    with Chem.SDWriter(f"results/docking/ligands/{name}.sdf") as w:
        w.write(mol)

# 3. Sanitizar Receptores a formato wwPDB estricto de 80 columnas (columnas 77-78 con elemento)
receptors = {
    "dry": "results/mcpb/SIRT6_Zn_dry.pdb",
    "c0": "results/clustering/centroid.c0.pdb",
    "c1": "results/clustering/centroid.c1.pdb",
    "c2": "results/clustering/centroid.c2.pdb",
    "c3": "results/clustering/centroid.c3.pdb"
}

for label, src_path in receptors.items():
    dst_path = f"results/docking/receptors/SIRT6_{label}.pdb"
    with open(src_path, "r") as fin, open(dst_path, "w") as fout:
        for raw_line in fin:
            if raw_line.startswith(("ATOM", "HETATM")):
                line = raw_line.rstrip("\r\n").ljust(76)
                resname = line[17:20].strip()
                atomname = line[12:16].strip()

                if resname in {"ZN1", "ZN"} or atomname.upper().startswith("ZN"):
                    line = "HETATM" + line[6:12] + " ZN   ZN" + line[20:54] + "  1.00  0.00          "
                    elem = "ZN"
                else:
                    if resname in {"CM1", "CM2", "CM3", "CM4", "CYM"}:
                        line = line[:17] + "CYS" + line[20:]
                    # Extraer elemento químico limpio (C, N, O, S, H) del nombre atómico AMBER
                    clean_name = "".join([c for c in atomname if not c.isdigit()])
                    elem = clean_name[0].upper() if clean_name else "C"

                fout.write(f"{line[:76]}{elem:>2s}  \n")
            elif raw_line.startswith(("TER", "END")):
                fout.write(raw_line)
    print(f"[✓] Receptor wwPDB (80 cols) sanitizado: {dst_path}")
