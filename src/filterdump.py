
#name = "visualize_M2000_N256_T0.2_r10_cang45.out"
#name = "visualizepull_M2000_N256_r10_cang60_t5000000_nve.out"
#name = "viscomp_M2000_N500.out"

name = "vis_sphere_stiff_M2000_N256_T0.1_r80.out"

filename = "../visfiles/" + name
epsilon = 0.0000000000001
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def RepresentsFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False


def get_von_mises(words, idx):
    if idx < 0: raise "Error"
    s_xx, s_yy, s_zz = float(words[idx]), float(words[idx+1]), float(words[idx+2])
    s_xy, s_xz, s_yz = float(words[idx+3]), float(words[idx+4]), float(words[idx+5])
    a = (s_xx-s_yy)**2 + (s_zz-s_xx)**2 + (s_yy-s_zz)**2
    b = 6*(s_xy**2 + s_xz**2 + s_yz**2)
    return ((a+b)/2)**(1/2)

def filter_for_ovito(filename, new_filename):
    newfile = open(new_filename, "w+")
    stress_idx = -1
    atom_vol = (4/3)*3.14*(1.12)**3
    with open(filename, 'r') as f:
        for line in f:
            words = line.split()
            if RepresentsInt(words[0]) and len(words) > 5: #check for small values
                for i in range(len(words)):
                    val = abs(float(words[i]))
                    if val < epsilon: words[i] = "0"
                von_mis_stress = str(get_von_mises(words, stress_idx)/atom_vol)
                words.append(von_mis_stress)
                newline = ' '.join(words) + "\n"
                newfile.write(newline)
            elif len(words) > 2 and words[1] == "ATOMS": 
                stress_idx = words.index("c_full_s[1]")
                words.append("von_mises")
                if words[5] != "x":
                    words[5] = "x"
                    words[6] = "y"
                    words[7] = "z"
                    words[8] = "fx"
                    words[9] = "fy"
                    words[10] = "fz"
                    newline = ' '.join(words) + "\n"
                    newfile.write(newline)
                else:
                    newfile.write(line)
            else:
                newfile.write(line)
    newfile.close()
            
def isolate_dist_changes(filename, new_filename):
    print("Isolating particles from " + filename + "\n to " + new_filename)
    newfile = open(new_filename, "w+")
    with open(filename, 'r') as f:
        new_lines = []
        atom_count = 0
        for line in f:
            if "TIMESTEP" in line: 
                #update file
                prev_line = ""
                for nline in new_lines:
                    if "NUMBER OF ATOMS" in prev_line:
                        nline = str(atom_count) + "\n"
                    newfile.write(nline)
                    prev_line = nline
                new_lines = []
                atom_count = 0
            words = line.split()
            if RepresentsInt(words[0]) and len(words) > 5:
                for i in range(len(words)-5, len(words)):
                    val = abs(float(words[i])) #check if number of bond changes isn't zero
                    if val != 0:
                        newline = ' '.join(words) + "\n"
                        new_lines.append(newline)
                        atom_count += 1
                        break                    
            elif len(words) > 2 and words[1] == "ATOMS":
                if words[5] != "x":
                    words[5] = "x"
                    words[6] = "y"
                    words[7] = "z"
                    words[8] = "fx"
                    words[9] = "fy"
                    words[10] = "fz"
                    newline = ' '.join(words) + "\n"
                    new_lines.append(newline)
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)
    newfile.close()


def main():
    M, N = 2000, 256
    T = 0.2
    r = 25
    cang = 45
    ovito = True
    if ovito:   
        name = "vis_sphere_stiff_M%d_N%d_T%g_r%d.out" %(M,N,T,r)
        filename = "../visfiles/" + name
        new_filename = "../visfiles/filt_" + name
        filter_for_ovito(filename, new_filename)
    else:
        Ts = [0.1, 0.2]
        drs = [0.5, 0.75, 1, 1.25, 1.5]
        for T in Ts:
            for dr in drs:
                name = "visualizechanges_M%d_N%d_T%g_r%d_cang%d_p%g.out" %(M,N,T,r,cang, dr)
                filename = "../visfiles/" + name
                new_filename = "../visfiles/filt_" + name
                isolate_dist_changes(filename, new_filename)
        


if __name__=="__main__":
    main()
