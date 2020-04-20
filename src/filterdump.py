
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


newfilename = "../visfiles/filt_" + name
newfile = open(newfilename, "w+")

with open(filename, 'r') as f:
    for line in f:
        words = line.split()
        if RepresentsInt(words[0]) and len(words) > 5:
            for i in range(len(words)):
                val = abs(float(words[i]))
                if val < epsilon: words[i] = "0"
            #fz, fy, fx = abs(float(words[-1])), abs(float(words[-2])), abs(float(words[-3]))
            #if fx < epsilon: words[-3] = "0"
            #if fy < epsilon: words[-2] = "0"
            #if fz < epsilon: words[-1] = "0"
            newline = ' '.join(words) + "\n"
            newfile.write(newline)
        elif len(words) > 2 and words[1] == "ATOMS":
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
            #if not ( RepresentsFloat(words[-1]) and RepresentsFloat(words[-2]) and RepresentsFloat(words[-3])):
            #    print(words[-3], words[-2], words[-3])
            
