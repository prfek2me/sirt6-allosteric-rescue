import argparse, os, subprocess, shutil, sys

def audit_and_relocate(inp_path, log_path, dat_path, gamess_dir):
    job_name = os.path.basename(inp_path).replace(".inp", "")
    expected_dat = os.path.join(os.path.abspath(gamess_dir), "userscr", f"{job_name}.dat")
    if os.path.exists(expected_dat): os.remove(expected_dat)

    print(f"[*] Detonando motor cuántico GAMESS: {job_name}...")
    cmd = f"csh {os.path.abspath(os.path.join(gamess_dir, 'rungms'))} {os.path.abspath(inp_path)} 00 1 1 0"
    with open(log_path, 'w') as f_log:
        subprocess.run(cmd, shell=True, stdout=f_log, stderr=subprocess.STDOUT, cwd=gamess_dir)

    with open(log_path, 'r') as f:
        if "TERMINATED NORMALLY" not in f.read():
            print(f"[!] COLAPSO CUÁNTICO en {log_path}")
            sys.exit(1)

    if not os.path.exists(expected_dat):
        print(f"[!] Fuga de matriz PUNCH en {expected_dat}")
        sys.exit(1)
    shutil.move(expected_dat, dat_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--dat", required=True)
    parser.add_argument("--gamess-dir", required=True)
    args = parser.parse_args()
    audit_and_relocate(args.inp, args.log, args.dat, args.gamess_dir)
