import sys

def extract(opt_log, fc_inp):
    with open(opt_log, 'r') as f: lines = f.readlines()
    start_idx = next(i for i in range(len(lines)-1, -1, -1) if "COORDINATES OF ALL ATOMS ARE (ANGS)" in lines[i])
    
    coords = []
    for line in lines[start_idx + 3:]:
        parts = line.split()
        if len(parts) == 5:
            try: coords.append(f"{parts[0]:<2} {parts[1]:>4} {parts[2]:>12} {parts[3]:>12} {parts[4]:>12}")
            except ValueError: break
        else: break
            
    header = " $CONTRL SCFTYP=RHF MULT=1 ICHARG=-2 MAXIT=200 COORD=UNIQUE RUNTYP=HESSIAN $END\n"
    header += " $DFT DFTTYP=B3LYP $END\n $SYSTEM MWORDS=250 $END\n $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n"
    header += " $DATA\nSIRT6_Zn_small_fc\nC1\n"
    
    with open(fc_inp, 'w') as f:
        f.write(header)
        f.write('\n'.join(coords) + '\n')
        f.write(" $END\n")

if __name__ == "__main__":
    extract(sys.argv[1], sys.argv[2])
