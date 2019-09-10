
#name = "visualize_M1000_N256_r10_cang45_npt.out"
name = "viscomp_M1000_N256.out"

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
            fz, fy, fx = abs(float(words[-1])), abs(float(words[-2])), abs(float(words[-3]))
            if fx < epsilon: words[-3] = "0"
            if fy < epsilon: words[-2] = "0"
            if fz < epsilon: words[-1] = "0"
            newline = ' '.join(words) + "\n"
            newfile.write(newline)
        else:
            newfile.write(line)
            #if not ( RepresentsFloat(words[-1]) and RepresentsFloat(words[-2]) and RepresentsFloat(words[-3])):
            #    print(words[-3], words[-2], words[-3])
            
