import os
import sys
import subprocess
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog('rdApp.warning')
os.makedirs("results/docking/poses", exist_ok=True)

RT = 0.59616  # kcal/mol a 300 K
raw_fracs = {}
with open("results/clustering/cluster_summary.dat", "r") as f:
    for line in f:
        if not line.startswith("#") and line.strip():
            parts = line.split()
            cid = int(parts[0])
            if cid in {0, 1, 2, 3}:
                raw_fracs[f"c{cid}"] = float(parts[2])

total_frac = sum(raw_fracs.values())
weights = {k: v / total_frac for k, v in raw_fracs.items()}
p0 = raw_fracs["c0"]
dG_conf = {k: -RT * np.log(v / p0) for k, v in raw_fracs.items()}
dG_conf["dry"] = 0.0
weights["dry"] = 0.0

ref_ubcs_sdf = "results/docking/ligands/UBCS039_crystal_ref.sdf"
ref_ubcs = Chem.SDMolSupplier(ref_ubcs_sdf, removeHs=True)[0]

receptors = ["dry", "c0", "c1", "c2", "c3"]
ligands = ["UBCS039", "MDL800", "NAM"]
records = []

env = os.environ.copy()
env["LD_LIBRARY_PATH"] = f"{sys.prefix}/lib:" + env.get("LD_LIBRARY_PATH", "")

for lig in ligands:
    lig_sdf = f"results/docking/ligands/{lig}.sdf"
    for rec in receptors:
        rec_pdb = f"results/docking/receptors/SIRT6_{rec}.pdb"
        out_sdf = f"results/docking/poses/{lig}_{rec}.sdf"
        log_txt = f"results/docking/poses/{lig}_{rec}.log"

        cmd = [
            "./bin/gnina",
            "-r", rec_pdb,
            "-l", lig_sdf,
            "--autobox_ligand", ref_ubcs_sdf,
            "--autobox_add", "8.0",
            "--cnn_scoring", "rescore",
            "--cpu", "6",
            "--exhaustiveness", "16",
            "--num_modes", "9",
            "--seed", "42",
            "-o", out_sdf,
            "--log", log_txt
        ]

        print(f"[*] Acoplando (6 CPU + GPU CNN-Rescore) {lig} -> SIRT6_{rec} ...")
        res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=180)
        if res.returncode != 0:
            print(f"[!] Error en GNINA ({lig} vs {rec}):\n{res.stderr}")
            sys.exit(1)

        suppl = [m for m in Chem.SDMolSupplier(out_sdf, removeHs=True) if m is not None]
        top_mol = suppl[0]
        vina_scores = [float(m.GetProp("minimizedAffinity")) for m in suppl]
        best_vina = min(vina_scores)
        top1_vina = float(top_mol.GetProp("minimizedAffinity"))
        cnn_score = float(top_mol.GetProp("CNNscore"))
        cnn_aff = float(top_mol.GetProp("CNNaffinity"))

        rmsd_top1 = np.nan
        rmsd_min = np.nan
        if lig == "UBCS039":
            rmsds = [rdMolAlign.CalcRMS(m, ref_ubcs) for m in suppl]
            rmsd_top1 = rmsds[0]
            rmsd_min = min(rmsds)

        records.append({
            "Ligand": lig,
            "Receptor": rec,
            "Pop_Weight": round(weights[rec], 3),
            "dG_conf_kcal": round(dG_conf[rec], 2),
            "Vina_Top1": round(top1_vina, 2),
            "Vina_Best": round(best_vina, 2),
            "E_eff_kcal": round(best_vina + dG_conf[rec], 2),
            "CNNscore": round(cnn_score, 3),
            "CNNaffinity": round(cnn_aff, 2),
            "Crystal_RMSD_Top1": round(rmsd_top1, 2) if not np.isnan(rmsd_top1) else None,
            "Crystal_RMSD_Best": round(rmsd_min, 2) if not np.isnan(rmsd_min) else None
        })

df = pd.DataFrame(records)
df.to_csv("results/docking/ensemble_detailed.csv", index=False)

print("\n=== TABLA DETALLADA DE ACOPLAMIENTO (DRY + ENSEMBLE c0-c3 CO-CENTRADO) ===")
print(df.to_string(index=False))

md_df = df[df["Receptor"] != "dry"].copy()
summary_rows = []
for lig, grp in md_df.groupby("Ligand", sort=False):
    w = grp["Pop_Weight"].values
    boltz_vina = np.sum(w * grp["Vina_Best"].values)
    boltz_cnn_aff = np.sum(w * grp["CNNaffinity"].values)
    best_idx = grp["E_eff_kcal"].idxmin()
    best_row = grp.loc[best_idx]
    summary_rows.append({
        "Ligand": lig,
        "Weighted_Vina_kcal": round(boltz_vina, 2),
        "Best_State": best_row["Receptor"],
        "Best_E_eff_kcal": best_row["E_eff_kcal"],
        "Weighted_CNNaffinity": round(boltz_cnn_aff, 2),
        "Best_State_CNNscore": best_row["CNNscore"]
    })

sum_df = pd.DataFrame(summary_rows)
sum_df.to_csv("results/docking/ensemble_summary.csv", index=False)

print("\n=== RESUMEN TERMODINÁMICO DEL ENSAMBLE (PONDERADO POR POBLACIÓN) ===")
print(sum_df.to_string(index=False))
