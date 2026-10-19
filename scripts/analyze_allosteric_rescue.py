import os
import re
import subprocess
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/analysis", exist_ok=True)

# 1. Mapear índice AMBER (1..282) -> Numeración cristalográfica oficial en 5MF6.pdb
amber_to_pdb = {}
seq_idx = 0
last_uid = None
with open("data/raw/5MF6.pdb") as f:
    for line in f:
        if line.startswith("ATOM") and line[21] == "A":
            uid = line[22:27]
            if uid != last_uid:
                seq_idx += 1
                last_uid = uid
                resname = line[17:20].strip()
                resnum = int(line[22:26].strip())
                amber_to_pdb[seq_idx] = f"{resname}{resnum}"

pocket_mask = ":44,45,47,48,51,52,53,54,55,60,61,62,71,72,73,77,104,105,106,107,108,124,127,144,146,147,148,149,169,170,171"
if os.path.exists("results/clustering/pocket_box.txt"):
    with open("results/clustering/pocket_box.txt") as f:
        for line in f:
            if line.startswith("pocket_mask"):
                pocket_mask = line.split("=")[1].strip()

pocket_resids = [int(x) for x in pocket_mask.replace(":", "").split(",")]

# 2. Procesar trayectorias en cpptraj y exportar DCD seco de la cuenca estacionaria (28-50 ns, frames 561-1000)
apo_in = f"""parm results/mcpb/SIRT6_Zn_solv.prmtop
trajin results/mcpb/SIRT6_Zn_prod_50ns.dcd 201 1000
autoimage anchor :1-283
strip :WAT,Na+,Cl-

rms :25-265@CA,C,N first
average crdset APO_AVG :1-282 start 361 stop 800
run

rms ZN_TETRA :CM1,CM2,CM3,CM4,ZN1&!@H= first out results/analysis/apo_rmsd.dat
rms CORE_RMSD :25-265@CA,C,N ref APO_AVG :25-265@CA,C,N out results/analysis/apo_rmsd.dat
rms POCKET_RMSD {pocket_mask}&!@H= ref APO_AVG {pocket_mask}&!@H= nofit out results/analysis/apo_rmsd.dat
run

clear trajin
trajin results/mcpb/SIRT6_Zn_prod_50ns.dcd 561 1000
autoimage anchor :1-283
strip :WAT,Na+,Cl-
rms :25-265@CA,C,N ref APO_AVG :25-265@CA,C,N
atomicfluct out results/analysis/apo_rmsf.dat :1-282@CA byres
matrix correl name DCCM_APO :1-282@CA out results/analysis/apo_dccm.dat
run
quit
"""

basin_dcd = "results/analysis/holo_basin_dry.dcd"
if os.path.exists(basin_dcd):
    os.remove(basin_dcd)

holo_in = f"""parm results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop
trajin results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd 201 1000
autoimage anchor :1-283
strip :WAT,Na+,Cl-

rms :25-265@CA,C,N first
average crdset HOLO_AVG :1-284 start 361 stop 800
run

rms ZN_TETRA :CM1,CM2,CM3,CM4,ZN1&!@H= first out results/analysis/holo_rmsd.dat
rms CORE_RMSD :25-265@CA,C,N ref HOLO_AVG :25-265@CA,C,N out results/analysis/holo_rmsd.dat
rms POCKET_RMSD {pocket_mask}&!@H= ref HOLO_AVG {pocket_mask}&!@H= nofit out results/analysis/holo_rmsd.dat
rms MDL800_POSE :MDL&!@H= ref HOLO_AVG :MDL&!@H= nofit out results/analysis/holo_rmsd.dat
rms MDL800_INTERNAL :MDL&!@H= ref HOLO_AVG :MDL&!@H= out results/analysis/holo_mdl_internal.dat
run

clear trajin
trajin results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd 561 1000
autoimage anchor :1-283
strip :WAT,Na+,Cl-
rms :25-265@CA,C,N ref HOLO_AVG :25-265@CA,C,N
atomicfluct out results/analysis/holo_rmsf.dat :1-282@CA byres
matrix correl name DCCM_HOLO :1-282@CA out results/analysis/holo_dccm.dat
hbond HB :1-282,MDL avgout results/analysis/holo_hbonds.dat nointramol
trajout {basin_dcd} dcd
run
quit
"""

with open("results/analysis/apo_analysis.in", "w") as f:
    f.write(apo_in)
with open("results/analysis/holo_analysis.in", "w") as f:
    f.write(holo_in)

print("[*] Procesando trayectorias Apo y Holo en la cuenca estacionaria (28-50 ns)...")
subprocess.run(["cpptraj", "-i", "results/analysis/apo_analysis.in"], check=True, stdout=subprocess.DEVNULL)
subprocess.run(["cpptraj", "-i", "results/analysis/holo_analysis.in"], check=True, stdout=subprocess.DEVNULL)

apo_rmsd = pd.read_csv("results/analysis/apo_rmsd.dat", sep=r"\s+")
holo_rmsd = pd.read_csv("results/analysis/holo_rmsd.dat", sep=r"\s+")
mdl_int = pd.read_csv("results/analysis/holo_mdl_internal.dat", sep=r"\s+")

apo_eq = apo_rmsd.iloc[360:].copy()
holo_eq = holo_rmsd.iloc[360:].copy()
mdl_eq = mdl_int.iloc[360:].copy()

# 3. Identificar y exportar el fotograma Medoide exacto de la cuenca estacionaria (561-1000)
combined_dist2 = holo_eq["POCKET_RMSD"]**2 + holo_eq["MDL800_POSE"]**2
best_idx_in_rmsd = int(combined_dist2.idxmin())  # índice 0..799
medoid_dcd_frame = 201 + best_idx_in_rmsd        # fotograma real en el DCD (561..1000)
medoid_time_ns = medoid_dcd_frame * 0.05

medoid_pdb = "results/analysis/SIRT6_Holo_MDL800_basin_rep.pdb"
medoid_in = f"""parm results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop [solv]
parm results/mcpb/SIRT6_Zn_dry.prmtop [dry]
reference results/mcpb/SIRT6_Zn_dry.inpcrd parm [dry] [crystal]
trajin results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd {medoid_dcd_frame} {medoid_dcd_frame} parm [solv]
autoimage anchor :1-283
strip :WAT,Na+,Cl-
rms :25-265@CA,C,N ref [crystal] :25-265@CA,C,N
trajout {medoid_pdb} pdb
run
quit
"""
with open("results/analysis/export_medoid.in", "w") as f:
    f.write(medoid_in)
subprocess.run(["cpptraj", "-i", "results/analysis/export_medoid.in"], check=True, stdout=subprocess.DEVNULL)
print(f"[✓] Confórmero medoide estacionario exportado (Frame {medoid_dcd_frame} = {medoid_time_ns:.2f} ns): {medoid_pdb}")

# 4. Descomposición Energética Exacta del Hamiltoniano AMBER (ΔE_vdW + ΔE_elec por residuo)
def parse_prmtop_sections(prmtop_path):
    sections = {}
    current_flag = None
    with open(prmtop_path, "r") as f:
        for line in f:
            if line.startswith("%FLAG"):
                current_flag = line.split()[1].strip()
                sections[current_flag] = []
            elif line.startswith("%FORMAT") or line.startswith("%VERSION") or current_flag is None:
                continue
            else:
                if current_flag in {"POINTERS", "CHARGE", "ATOM_TYPE_INDEX", "NONBONDED_PARM_INDEX",
                                    "LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF", "RESIDUE_POINTER"}:
                    sections[current_flag].extend(line.split())
    return sections

print("[*] Calculando descomposición energética AMBER (Lennard-Jones 12-6 + Coulomb) por residuo...")
sec = parse_prmtop_sections("results/holo_mdl800/SIRT6_Holo_MDL800_dry.prmtop")
pointers = [int(x) for x in sec["POINTERS"]]
n_types = pointers[1]
charges = np.array([float(x) for x in sec["CHARGE"]])
atom_types = np.array([int(x) - 1 for x in sec["ATOM_TYPE_INDEX"]])
nb_index = np.array([int(x) - 1 for x in sec["NONBONDED_PARM_INDEX"]])
acoef = np.array([float(x) for x in sec["LENNARD_JONES_ACOEF"]])
bcoef = np.array([float(x) for x in sec["LENNARD_JONES_BCOEF"]])
res_ptr = [int(x) - 1 for x in sec["RESIDUE_POINTER"]] + [pointers[0]]

# Índices atómicos: Proteína (residuos 1..282 -> átomos 0..4429), MDL800 (residuo 284 -> átomos 4431..4481)
prot_end = res_ptr[282]
lig_start, lig_end = res_ptr[283], res_ptr[284]

q_prot = charges[:prot_end]
q_lig = charges[lig_start:lig_end]
t_prot = atom_types[:prot_end]
t_lig = atom_types[lig_start:lig_end]

pair_nb_idx = nb_index[t_prot[:, None] * n_types + t_lig[None, :]]
A_mat = acoef[pair_nb_idx]
B_mat = bcoef[pair_nb_idx]
Q_mat = q_prot[:, None] * q_lig[None, :]

u_basin = mda.Universe("results/holo_mdl800/SIRT6_Holo_MDL800_dry.prmtop", basin_dcd)
n_frames = len(u_basin.trajectory)
atom_vdw_sum = np.zeros(prot_end)
atom_elec_sum = np.zeros(prot_end)

for ts in u_basin.trajectory:
    pos = ts.positions
    diff = pos[:prot_end, None, :] - pos[lig_start:lig_end][None, :, :]
    r2 = np.sum(diff * diff, axis=-1)
    inv_r = 1.0 / np.sqrt(r2)
    inv_r6 = (1.0 / r2) ** 3
    inv_r12 = inv_r6 * inv_r6

    e_vdw = A_mat * inv_r12 - B_mat * inv_r6
    e_elec = Q_mat * inv_r

    atom_vdw_sum += np.sum(e_vdw, axis=1)
    atom_elec_sum += np.sum(e_elec, axis=1)

atom_vdw_mean = atom_vdw_sum / n_frames
atom_elec_mean = atom_elec_sum / n_frames

res_energy_rows = []
for r_idx in range(282):
    a0, a1 = res_ptr[r_idx], res_ptr[r_idx + 1]
    evdw = float(np.sum(atom_vdw_mean[a0:a1]))
    eelec = float(np.sum(atom_elec_mean[a0:a1]))
    amber_r = r_idx + 1
    res_energy_rows.append({
        "AMBER_Res": amber_r,
        "PDB_5MF6": amber_to_pdb.get(amber_r, f"RES{amber_r}"),
        "E_vdW_kcal": round(evdw, 3),
        "E_elec_kcal": round(eelec, 3),
        "E_total_kcal": round(evdw + eelec, 3),
        "In_Allosteric_Pocket": amber_r in pocket_resids
    })

ene_df = pd.DataFrame(res_energy_rows).sort_values("E_total_kcal")
ene_df.to_csv("results/analysis/mdl800_per_residue_energy.csv", index=False)

# 5. Compilar tablas de RMSF, DCCM y Resumen Biofísico
apo_rmsf = pd.read_csv("results/analysis/apo_rmsf.dat", sep=r"\s+")
holo_rmsf = pd.read_csv("results/analysis/holo_rmsf.dat", sep=r"\s+")

rmsf_df = pd.DataFrame({
    "AMBER_Res": apo_rmsf.iloc[:, 0].astype(int),
    "PDB_5MF6": [amber_to_pdb.get(int(r), f"RES{int(r)}") for r in apo_rmsf.iloc[:, 0]],
    "RMSF_Apo": apo_rmsf.iloc[:, 1],
    "RMSF_Holo": holo_rmsf.iloc[:, 1]
})
rmsf_df["Delta_RMSF"] = rmsf_df["RMSF_Holo"] - rmsf_df["RMSF_Apo"]
rmsf_df["In_Allosteric_Pocket"] = rmsf_df["AMBER_Res"].isin(pocket_resids)
rmsf_df = rmsf_df.merge(ene_df[["AMBER_Res", "E_vdW_kcal", "E_elec_kcal", "E_total_kcal"]], on="AMBER_Res")
rmsf_df.to_csv("results/analysis/rmsf_comparison.csv", index=False)

dccm_apo = np.loadtxt("results/analysis/apo_dccm.dat")
dccm_holo = np.loadtxt("results/analysis/holo_dccm.dat")
delta_dccm = dccm_holo - dccm_apo

summary_stats = pd.DataFrame([
    {
        "Metric (Stationary Basin 28-50 ns)": "Folded Core Backbone RMSD (:25-265) (Å)",
        "Apo (Mean ± SD)": f"{apo_eq['CORE_RMSD'].mean():.2f} ± {apo_eq['CORE_RMSD'].std():.2f}",
        "Holo-MDL800 (Mean ± SD)": f"{holo_eq['CORE_RMSD'].mean():.2f} ± {holo_eq['CORE_RMSD'].std():.2f}"
    },
    {
        "Metric (Stationary Basin 28-50 ns)": "Zn-S4 Quantum Tetrahedron RMSD (Å)",
        "Apo (Mean ± SD)": f"{apo_eq['ZN_TETRA'].mean():.2f} ± {apo_eq['ZN_TETRA'].std():.2f}",
        "Holo-MDL800 (Mean ± SD)": f"{holo_eq['ZN_TETRA'].mean():.2f} ± {holo_eq['ZN_TETRA'].std():.2f}"
    },
    {
        "Metric (Stationary Basin 28-50 ns)": "Allosteric Pocket Basin RMSD (Å)",
        "Apo (Mean ± SD)": f"{apo_eq['POCKET_RMSD'].mean():.2f} ± {apo_eq['POCKET_RMSD'].std():.2f}",
        "Holo-MDL800 (Mean ± SD)": f"{holo_eq['POCKET_RMSD'].mean():.2f} ± {holo_eq['POCKET_RMSD'].std():.2f}"
    },
    {
        "Metric (Stationary Basin 28-50 ns)": "MDL800 Pose RMSD in Basin (nofit) (Å)",
        "Apo (Mean ± SD)": "N/A",
        "Holo-MDL800 (Mean ± SD)": f"{holo_eq['MDL800_POSE'].mean():.2f} ± {holo_eq['MDL800_POSE'].std():.2f}"
    },
    {
        "Metric (Stationary Basin 28-50 ns)": "MDL800 Internal Torsional RMSD (Å)",
        "Apo (Mean ± SD)": "N/A",
        "Holo-MDL800 (Mean ± SD)": f"{mdl_eq['MDL800_INTERNAL'].mean():.2f} ± {mdl_eq['MDL800_INTERNAL'].std():.2f}"
    },
    {
        "Metric (Stationary Basin 28-50 ns)": "Total Protein-MDL800 Non-Bonded Energy (kcal/mol)",
        "Apo (Mean ± SD)": "N/A",
        "Holo-MDL800 (Mean ± SD)": f"{ene_df['E_total_kcal'].sum():.2f} (vdW: {ene_df['E_vdW_kcal'].sum():.2f}, Elec: {ene_df['E_elec_kcal'].sum():.2f})"
    }
])
summary_stats.to_csv("results/analysis/biophysical_summary.csv", index=False)

print("\n=== RESUMEN BIOFÍSICO Y ENERGÉTICO EN CUENCA ESTACIONARIA (28-50 ns) ===")
print(summary_stats.to_string(index=False))

print("\n=== TOP 12 RESIDUOS DE ANCLAJE ENERGÉTICO DE MDL800 (ΔE_vdW + ΔE_elec) ===")
print(ene_df.head(12).to_string(index=False))

# 6. Generar Dashboard de Publicación de 4 Paneles
time_ns = np.linspace(10.05, 50.0, len(apo_rmsd))
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])

# Panel A: RMSD y Medoide
ax1 = fig.add_subplot(gs[0, 0])
ax1.axvspan(10.0, 28.0, color="lightgray", alpha=0.35, label="Induced-Fit Relaxation (10-28 ns)")
ax1.plot(time_ns, apo_rmsd["POCKET_RMSD"], label="Apo Pocket (vs Basin Avg)", color="slategray", alpha=0.75)
ax1.plot(time_ns, holo_rmsd["POCKET_RMSD"], label="Holo Pocket (vs Basin Avg)", color="navy", linewidth=1.4)
ax1.plot(time_ns, holo_rmsd["MDL800_POSE"], label="MDL800 Pose (vs Basin Avg)", color="crimson", alpha=0.85)
ax1.plot(time_ns, holo_rmsd["ZN_TETRA"], label="Holo Zn-S4 Cluster", color="darkorange", linewidth=1.2)
ax1.axvline(28.0, color="black", linestyle="--", linewidth=1.0)
ax1.scatter([medoid_time_ns], [holo_rmsd.loc[best_idx_in_rmsd, "MDL800_POSE"]], color="gold", edgecolor="black", s=70, zorder=6, label=f"Basin Medoid ({medoid_time_ns:.1f} ns)")
ax1.set_xlabel("Simulation Time (ns)", fontweight="bold")
ax1.set_ylabel("RMSD vs. Stationary Basin Avg (Å)", fontweight="bold")
ax1.set_title("A. Induced-Fit Convergence & Basin Stability (28-50 ns)", fontweight="bold")
ax1.legend(loc="upper right", frameon=True, fontsize=8)
ax1.grid(alpha=0.25)

# Panel B: ΔDCCM
ax2 = fig.add_subplot(gs[0, 1])
im = ax2.imshow(delta_dccm, cmap="seismic", vmin=-0.6, vmax=0.6, origin="lower", extent=[1, 282, 1, 282])
ax2.set_xlabel("AMBER Residue Index (1-282)", fontweight="bold")
ax2.set_ylabel("AMBER Residue Index (1-282)", fontweight="bold")
ax2.set_title("B. Stationary Allosteric Reorganization (ΔDCCM: 28-50 ns)", fontweight="bold")
cbar = plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
cbar.set_label("ΔC_ij (Holo - Apo)", fontweight="bold")

# Panel C: ΔRMSF
ax3 = fig.add_subplot(gs[1, 0])
res_x = rmsf_df["AMBER_Res"].values
d_rmsf = rmsf_df["Delta_RMSF"].values
ax3.plot(res_x, d_rmsf, color="black", linewidth=1.1, label="ΔRMSF (Holo - Apo)")
ax3.fill_between(res_x, d_rmsf, 0, where=(d_rmsf < 0), color="royalblue", alpha=0.4, label="Hinge Rigidity / Damping (ΔRMSF < 0)")
ax3.fill_between(res_x, d_rmsf, 0, where=(d_rmsf >= 0), color="indianred", alpha=0.35, label="Allosteric Plasticity (ΔRMSF > 0)")
ax3.scatter(
    rmsf_df[rmsf_df["In_Allosteric_Pocket"]]["AMBER_Res"],
    rmsf_df[rmsf_df["In_Allosteric_Pocket"]]["Delta_RMSF"],
    color="gold", edgecolor="black", zorder=5, s=30, label="Allosteric Pocket (31 res)"
)
ax3.axhline(0, color="black", linestyle="--", linewidth=0.8)
ax3.set_xlim(20, 270)
ax3.set_ylim(-0.6, 2.0)
ax3.set_xlabel("AMBER Residue Index (20-270)", fontweight="bold")
ax3.set_ylabel("ΔRMSF (Å)", fontweight="bold")
ax3.set_title("C. Stationary Thermal Fluctuation Modulation (28-50 ns)", fontweight="bold")
ax3.legend(loc="upper right", frameon=True, fontsize=8)
ax3.grid(alpha=0.25)

# Panel D: Descomposición Energética por Residuo (Top 10 Anclas de MDL800)
ax4 = fig.add_subplot(gs[1, 1])
top10_ene = ene_df.head(10).iloc[::-1]
y_pos = np.arange(len(top10_ene))
labels = [f"{r['PDB_5MF6']} (:{r['AMBER_Res']})" for _, r in top10_ene.iterrows()]
ax4.barh(y_pos, top10_ene["E_vdW_kcal"], color="darkslateblue", alpha=0.85, label="ΔE_vdW (London Dispersion / Hydrophobic)")
ax4.barh(y_pos, top10_ene["E_elec_kcal"], left=np.minimum(top10_ene["E_vdW_kcal"], 0), color="darkcyan", alpha=0.85, label="ΔE_elec (Coulomb / H-Bonds)")
ax4.set_yticks(y_pos)
ax4.set_yticklabels(labels, fontweight="bold")
ax4.axvline(0, color="black", linewidth=0.8)
ax4.set_xlabel("Mean Interaction Energy with MDL800 (kcal/mol)", fontweight="bold")
ax4.set_title("D. Top 10 Energetic Anchor Residues for MDL800 (AMBER)", fontweight="bold")
ax4.legend(loc="lower left", frameon=True, fontsize=8.5)
ax4.grid(alpha=0.25, axis="x")

plt.tight_layout()
out_fig = "results/analysis/allosteric_rescue_dashboard.png"
plt.savefig(out_fig, dpi=300)
print(f"\n[✓] Dashboard biofísico y energético de 4 paneles generado en: {out_fig}")
