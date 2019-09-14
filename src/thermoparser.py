import matplotlib.pyplot as plt
import math

filename = "../visfiles/compress_M1000_N256.txt"
vz = 0.0001
dt = 0.005
area = 5184
D = 10
Lz = 52.8 - 10
Lx, Ly = 72.1, 72.1


displ = []
fzs   = []
z0 = 24.9
with open(filename, 'r') as file:
    for line in file:
        line = file.readline()
        line = line.strip()
        print(line)
        arr = line.split()
        print(len(arr))
        if len(arr) != 15 or not arr[0].isdigit(): continue
        d =  float(arr[0]) * vz * dt
        Lz_cur = Lz - d
        zhi = float(arr[8])
        pzz = float(arr[11])
        if pzz > 3 or zhi > z0: continue
        lambdaz = zhi/z0
        strainz = abs(math.log(lambdaz))
        displ.append(strainz)
        fzs.append(pzz)
        '''lambdaz = Lz_cur / Lz
        lambdax = 1 / lambdaz**(1/2)
        lambday = lambdax
        Lx_cur = lambdax * Lx
        Ly_cur = lambday * Ly
        strainz = abs(math.log(lambdaz))
        displ.append(d)
        fz_top = float(arr[10])
        fz_bot = float(arr[13])
        fzs.append( fz_top/ (Lx_cur * Ly_cur))'''

plt.plot(displ, fzs)
plt.show()
