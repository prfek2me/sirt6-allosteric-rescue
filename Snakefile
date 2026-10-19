import yaml

configfile: "config/config.yaml"
PDB_ID = config["system"]["pdb_id"]

rule all:
    input:
        "results/docking/ensemble_detailed.csv",
        "results/docking/ensemble_summary.csv",
        "results/holo_mdl800/SIRT6_Holo_prod_final.xml",
        "results/analysis/biophysical_summary.csv",
        "results/analysis/mdl800_per_residue_energy.csv",
        "results/analysis/SIRT6_Holo_MDL800_basin_rep.pdb",
        "results/analysis/allosteric_rescue_dashboard.png"

rule prep_metal_center:
    input:
        f"data/raw/{PDB_ID}.pdb"
    output:
        apo_pdb="results/SIRT6_apo_cym.pdb",
        opt_inp="results/mcpb/SIRT6_Zn_small_opt.inp",
        mk_inp="results/mcpb/SIRT6_Zn_large_mk.inp"
    shell:
        "python scripts/prep_mcpb_fase1.py --input-pdb {input} --out-pdb {output.apo_pdb} --mcpb-dir results/mcpb"

rule gamess_small_opt:
    input:
        inp="results/mcpb/SIRT6_Zn_small_opt.inp"
    output:
        log="results/mcpb/SIRT6_Zn_small_opt.log",
        dat="results/mcpb/SIRT6_Zn_small_opt.dat"
    shell:
        "python scripts/run_gamess_fase2.py --inp {input.inp} --log {output.log} --dat {output.dat} --gamess-dir bin/gamess"

rule extract_hessian:
    input:
        opt_log="results/mcpb/SIRT6_Zn_small_opt.log"
    output:
        fc_inp="results/mcpb/SIRT6_Zn_small_fc.inp"
    shell:
        "python scripts/extract_hessian_inp.py {input.opt_log} {output.fc_inp}"

rule gamess_small_fc:
    input:
        fc_inp="results/mcpb/SIRT6_Zn_small_fc.inp"
    output:
        log="results/mcpb/SIRT6_Zn_small_fc.log",
        dat="results/mcpb/SIRT6_Zn_small_fc.dat"
    shell:
        "python scripts/run_gamess_fase2.py --inp {input.fc_inp} --log {output.log} --dat {output.dat} --gamess-dir bin/gamess"

rule gamess_large_mk:
    input:
        inp="results/mcpb/SIRT6_Zn_large_mk.inp"
    output:
        log="results/mcpb/SIRT6_Zn_large_mk.log",
        dat="results/mcpb/SIRT6_Zn_large_mk.dat"
    shell:
        "python scripts/run_gamess_fase2.py --inp {input.inp} --log {output.log} --dat {output.dat} --gamess-dir bin/gamess"

rule mcpb_seminario:
    input:
        opt_log="results/mcpb/SIRT6_Zn_small_opt.log",
        fc_log="results/mcpb/SIRT6_Zn_small_fc.log",
        mk_log="results/mcpb/SIRT6_Zn_large_mk.log",
        mk_dat="results/mcpb/SIRT6_Zn_large_mk.dat"
    output:
        frcmod="results/mcpb/SIRT6_Zn_mcpbpy.frcmod"
    shell:
        "cd results/mcpb && MCPB.py -i SIRT6_Zn.in -s 2"

rule mcpb_resp_charges:
    input:
        frcmod="results/mcpb/SIRT6_Zn_mcpbpy.frcmod"
    output:
        zn_mol2="results/mcpb/ZN1.mol2"
    shell:
        "cd results/mcpb && MCPB.py -i SIRT6_Zn.in -s 3"

rule mcpb_tleap_router:
    input:
        zn_mol2="results/mcpb/ZN1.mol2"
    output:
        tleap_in="results/mcpb/SIRT6_Zn_tleap.in"
    shell:
        "cd results/mcpb && MCPB.py -i SIRT6_Zn.in -s 4"

rule build_amber_topology:
    input:
        tleap_in="results/mcpb/SIRT6_Zn_tleap.in"
    output:
        prmtop="results/mcpb/SIRT6_Zn_solv.prmtop",
        inpcrd="results/mcpb/SIRT6_Zn_solv.inpcrd",
        dry_pdb="results/mcpb/SIRT6_Zn_dry.pdb"
    shell:
        "cd results/mcpb && tleap -s -f SIRT6_Zn_tleap.in > tleap_build.log"

rule apo_equilibration:
    input:
        prmtop="results/mcpb/SIRT6_Zn_solv.prmtop",
        inpcrd="results/mcpb/SIRT6_Zn_solv.inpcrd"
    output:
        equil_xml="results/mcpb/SIRT6_Zn_equil.xml",
        equil_dcd="results/mcpb/SIRT6_Zn_equil.dcd"
    shell:
        "python scripts/run_apo_equilibration.py"

rule apo_production:
    input:
        equil_xml="results/mcpb/SIRT6_Zn_equil.xml"
    output:
        prod_xml="results/mcpb/SIRT6_Zn_prod_final.xml",
        prod_dcd="results/mcpb/SIRT6_Zn_prod_50ns.dcd"
    shell:
        "python scripts/run_apo_production.py > results/mcpb/production_telemetry.log"

rule cluster_apo_trajectory:
    input:
        raw_pdb=f"data/raw/{PDB_ID}.pdb",
        dry_pdb="results/mcpb/SIRT6_Zn_dry.pdb",
        prmtop="results/mcpb/SIRT6_Zn_solv.prmtop",
        dcd="results/mcpb/SIRT6_Zn_prod_50ns.dcd"
    output:
        centroids=expand("results/clustering/centroid.c{i}.pdb", i=[0, 1, 2, 3, 4]),
        box="results/clustering/pocket_box.txt",
        summary="results/clustering/cluster_summary.dat",
        rmsd="results/clustering/rmsd_clean.dat"
    shell:
        "python scripts/run_clustering.py"

rule prep_docking:
    input:
        raw_pdb=f"data/raw/{PDB_ID}.pdb",
        dry_pdb="results/mcpb/SIRT6_Zn_dry.pdb",
        centroids=expand("results/clustering/centroid.c{i}.pdb", i=[0, 1, 2, 3])
    output:
        ref_sdf="results/docking/ligands/UBCS039_crystal_ref.sdf",
        ligs=expand("results/docking/ligands/{lig}.sdf", lig=["UBCS039", "MDL800", "NAM"]),
        recs=expand("results/docking/receptors/SIRT6_{rec}.pdb", rec=["dry", "c0", "c1", "c2", "c3"])
    shell:
        "python scripts/prep_docking.py"

rule ensemble_docking:
    input:
        ref_sdf="results/docking/ligands/UBCS039_crystal_ref.sdf",
        ligs=expand("results/docking/ligands/{lig}.sdf", lig=["UBCS039", "MDL800", "NAM"]),
        recs=expand("results/docking/receptors/SIRT6_{rec}.pdb", rec=["dry", "c0", "c1", "c2", "c3"]),
        summary="results/clustering/cluster_summary.dat"
    output:
        detailed="results/docking/ensemble_detailed.csv",
        ens_sum="results/docking/ensemble_summary.csv",
        poses=expand("results/docking/poses/{lig}_{rec}.sdf", lig=["UBCS039", "MDL800", "NAM"], rec=["dry", "c0", "c1", "c2", "c3"])
    shell:
        "python scripts/run_ensemble_docking.py"

rule prep_holo_mdl800:
    input:
        pose="results/docking/poses/MDL800_c1.sdf",
        c1="results/clustering/centroid.c1.pdb",
        mcpb_pdb="results/mcpb/SIRT6_Zn_mcpbpy.pdb",
        tleap_in="results/mcpb/SIRT6_Zn_tleap.in"
    output:
        mol2="results/holo_mdl800/MDL800.mol2",
        frcmod="results/holo_mdl800/MDL800.frcmod",
        prmtop="results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop",
        inpcrd="results/holo_mdl800/SIRT6_Holo_MDL800_solv.inpcrd",
        dry_prmtop="results/holo_mdl800/SIRT6_Holo_MDL800_dry.prmtop"
    shell:
        "python scripts/prep_holo_mdl800.py"

rule holo_equilibration:
    input:
        prmtop="results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop",
        inpcrd="results/holo_mdl800/SIRT6_Holo_MDL800_solv.inpcrd"
    output:
        equil_xml="results/holo_mdl800/SIRT6_Holo_equil.xml",
        equil_dcd="results/holo_mdl800/SIRT6_Holo_equil.dcd"
    shell:
        "python scripts/run_holo_equilibration.py"

rule holo_production:
    input:
        prmtop="results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop",
        equil_xml="results/holo_mdl800/SIRT6_Holo_equil.xml"
    output:
        prod_xml="results/holo_mdl800/SIRT6_Holo_prod_final.xml",
        prod_dcd="results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd"
    shell:
        "python scripts/run_holo_production.py >> results/holo_mdl800/holo_production_telemetry.log"

rule analyze_allosteric_rescue:
    input:
        apo_prmtop="results/mcpb/SIRT6_Zn_solv.prmtop",
        apo_dcd="results/mcpb/SIRT6_Zn_prod_50ns.dcd",
        holo_prmtop="results/holo_mdl800/SIRT6_Holo_MDL800_solv.prmtop",
        holo_dry_prmtop="results/holo_mdl800/SIRT6_Holo_MDL800_dry.prmtop",
        holo_dcd="results/holo_mdl800/SIRT6_Holo_prod_50ns.dcd",
        box="results/clustering/pocket_box.txt"
    output:
        summary="results/analysis/biophysical_summary.csv",
        ene="results/analysis/mdl800_per_residue_energy.csv",
        rep_pdb="results/analysis/SIRT6_Holo_MDL800_basin_rep.pdb",
        rmsf="results/analysis/rmsf_comparison.csv",
        fig="results/analysis/allosteric_rescue_dashboard.png"
    shell:
        "python scripts/analyze_allosteric_rescue.py"
