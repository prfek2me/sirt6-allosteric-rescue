import argparse
import os
import subprocess
import sys
import glob
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def isolate_chain_a(input_pdb, temp_pdb):
    valid_res = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
                 "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", 
                 "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ZN"}
    
    with open(input_pdb, 'r') as fin, open(temp_pdb, 'w') as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                alt_loc = line[16]
                resname = line[17:20].strip()
                chain = line[21]
                if chain != 'A' or (alt_loc not in [' ', 'A']): continue
                if alt_loc == 'A': line = line[:16] + ' ' + line[17:]
                if resname == "ZN": line = "HETATM" + line[6:17] + "ZN " + line[20:]
                elif resname not in valid_res: continue
                fout.write(line)
            elif line.startswith("TER") and 'chain' in locals() and chain == 'A':
                fout.write(line)

def protonate_structure(temp_in, temp_out):
    fixer = PDBFixer(filename=temp_in)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(temp_out, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

def enforce_cym_and_renumber(protonated_pdb, final_pdb, cym_resids):
    cym_count = set()
    current_seq_id = 0
    last_res_key = None
    with open(protonated_pdb, 'r') as fin, open(final_pdb, 'w') as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resid = line[22:26].strip()
                icode = line[26]
                atom_name = line[12:16].strip()
                res_key = resname + resid + icode
                if res_key != last_res_key:
                    current_seq_id += 1
                    last_res_key = res_key
                if resname == "CYS" and resid in cym_resids:
                    if atom_name == "HG": continue 
                    line = line[:17] + "CYM" + line[20:]
                    cym_count.add(resid)
                new_resid_str = f"{current_seq_id:>4}"
                line = line[:22] + new_resid_str + line[26:]
                fout.write(line)
            elif line.startswith("TER"):
                if len(line) >= 26:
                    new_resid_str = f"{current_seq_id:>4}"
                    line = line[:22] + new_resid_str + line[26:]
                fout.write(line)
            else:
                fout.write(line)
    print(f"[*] Cisteínas mutadas a CYM: {len(cym_count)}/4 ({', '.join(cym_count)})")

def get_zinc_atom_id(final_pdb):
    with open(final_pdb, 'r') as fin:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")) and "ZN" in line[17:20]:
                return line[6:11].strip()
    return None

def translate_g09_to_gms(prefix, inp_file, model_type):
    possible_files = glob.glob(f"{prefix}.*")
    com_file = next((f for f in possible_files if f.endswith(('.com', '.gau', '.mkf'))), None)
    if not com_file: return False

    with open(com_file, 'r') as f: lines = f.read().split('\n')
    coord_start = 0
    charge, mult = "0", "1"
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 2 and parts[0].lstrip('+-').isdigit() and parts[1].lstrip('+-').isdigit():
            charge, mult = parts
            coord_start = i + 1
            break
            
    atomic_nums = {'C':'6.0', 'H':'1.0', 'O':'8.0', 'N':'7.0', 'S':'16.0', 'ZN':'30.0'}
    coords = []
    for line in lines[coord_start:]:
        if not line.strip(): break
        parts = line.split()
        if len(parts) == 4:
            elem = parts[0][:2].upper() if parts[0][:2].upper() == 'ZN' else parts[0][0].upper()
            anum = atomic_nums.get(elem, '0.0')
            coords.append(f"{elem:<2} {anum:>4} {parts[1]:>12} {parts[2]:>12} {parts[3]:>12}")

    scftyp = "RHF" if mult == "1" else "UHF"
    base_header = f" $CONTRL SCFTYP={scftyp} MULT={mult} ICHARG={charge} MAXIT=200 COORD=UNIQUE "
    
    if model_type == "small":
        header = base_header + "RUNTYP=OPTIMIZE $END\n"
        data_name = "SIRT6_Zn_small_opt"
    else:
        header = base_header + "RUNTYP=ENERGY $END\n"
        header += " $ELPOT IEPOT=1 WHERE=PDC $END\n $PDC PTSEL=GEODESIC $END\n"
        data_name = "SIRT6_Zn_large_mk"
        
    header += " $DFT DFTTYP=B3LYP $END\n $SYSTEM MWORDS=250 $END\n $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n"
    header += f" $DATA\n{data_name}\nC1\n"
    
    with open(inp_file, 'w') as f:
        f.write(header)
        f.write('\n'.join(coords) + '\n')
        f.write(" $END\n")
    os.remove(com_file)
    print(f"[+] Traducción léxica completada: {os.path.basename(inp_file)}")
    return True

def run_mcpb_pipeline(final_prot_pdb, out_dir):
    abs_out = os.path.abspath(out_dir)
    os.makedirs(abs_out, exist_ok=True)
    subprocess.run(["cp", final_prot_pdb, abs_out])

    zn_atom_id = get_zinc_atom_id(final_prot_pdb)
    print(f"[*] Ion Zinc anclado al ATOM ID estricto: {zn_atom_id}")

    zn_pdb = os.path.join(abs_out, "ZN.pdb")
    with open(final_prot_pdb, 'r') as fin, open(zn_pdb, 'w') as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")) and "ZN" in line[17:20]:
                fout.write(line)
                break 

    subprocess.run(["metalpdb2mol2.py", "-i", "ZN.pdb", "-o", "ZN.mol2", "-c", "2"], cwd=abs_out, check=True, capture_output=True)

    mcpb_in = os.path.join(abs_out, "SIRT6_Zn.in")
    content = f"original_pdb {os.path.basename(final_prot_pdb)}\ngroup_name SIRT6_Zn\ncut_off 2.8\nion_ids {zn_atom_id}\nion_mol2files ZN.mol2\nsoftware_version g09\n"
    with open(mcpb_in, "w") as f: f.write(content)

    print("[*] Ejecutando MCPB.py -s 1 (Modo Spoofing)...")
    res = subprocess.run(["MCPB.py", "-i", "SIRT6_Zn.in", "-s", "1"], cwd=abs_out, capture_output=True, text=True)

    large_ok = translate_g09_to_gms(os.path.join(abs_out, "SIRT6_Zn_large_mk"), os.path.join(abs_out, "SIRT6_Zn_large_mk.inp"), "large")
    small_ok = translate_g09_to_gms(os.path.join(abs_out, "SIRT6_Zn_small_opt"), os.path.join(abs_out, "SIRT6_Zn_small_opt.inp"), "small")

    if not large_ok or not small_ok:
        print("\n[!] MCPB.py ABORTÓ SILENCIOSAMENTE. Radiografía:\n", res.stdout)
        sys.exit(1)

    with open(mcpb_in, "w") as f: f.write(content.replace("software_version g09", "software_version gms"))
    print(f"[✓] Clústeres nativos para GAMESS listos en {abs_out}/")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", required=True)
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument("--mcpb-dir", default="results/mcpb")
    args = parser.parse_args()

    os.makedirs(args.mcpb_dir, exist_ok=True)
    t_apo = os.path.join(args.mcpb_dir, "t_apo.pdb")
    t_prot = os.path.join(args.mcpb_dir, "t_prot.pdb")
    cym_resids = ["141", "144", "166", "177"]

    isolate_chain_a(args.input_pdb, t_apo)
    protonate_structure(t_apo, t_prot)
    enforce_cym_and_renumber(t_prot, args.out_pdb, cym_resids)
    run_mcpb_pipeline(args.out_pdb, args.mcpb_dir)

    os.remove(t_apo)
    os.remove(t_prot)
