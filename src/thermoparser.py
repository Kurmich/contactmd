import matplotlib.pyplot as plt
import math

M, N = 2000, 500
filename = "../visfiles/deform_M%d_N%d.txt" %(M, N)

params = {}
params[(1000,256)] = [72.1, 72.1, 52.8]
params[(2000,256)] = [90.84, 90.84, 66.52]
params[(2000,500)] = [113.93, 113.93, 82]
vz = 0.0001
dt = 0.005
area = 5184
D = 10
Lz = params[(M,N)][2] - 10
Lx, Ly = params[(M,N)][0], params[(M,N)][1]


displ = []
fzs   = []
fys = []
fxs = []
z0 = 24.9
with open(filename, 'r') as file:
    for line in file:
        line = file.readline()
        line = line.strip()
        arr = line.split()
        
        if len(arr) != 15 or not arr[0].isdigit(): continue
        zhi = float(arr[8])
        z0 = max(z0, zhi)
        pzz = float(arr[11])
        if pzz > 3 or zhi > z0: continue
        lambdaz = zhi/z0
        strainz = abs(math.log(lambdaz))
        strain_cutoff = 0.06
        if strainz > 0.50 or strainz < strain_cutoff: continue
        print(line)
        print(len(arr))
        displ.append(strainz)
        fzs.append(pzz)
        '''
        print(line)
        print(len(arr))
        if len(arr) != 14 or not arr[0].isdigit(): continue
        d =  float(arr[0]) * vz * dt
        Lz_cur = Lz - d
        lambdaz = Lz_cur / Lz
        lambdax = 1 / lambdaz**(1/2)
        lambday = lambdax
        Lx_cur = lambdax * Lx
        Ly_cur = lambday * Ly
        strainz = abs(math.log(lambdaz))
        displ.append(lambdaz)
        fz_top = float(arr[10])
        fy = float(arr[9])
        fx = float(arr[8])
        fz_bot = float(arr[13])
        fzs.append( fz_top/ (Lx_cur * Ly_cur))
        fxs.append(fx / (Lz_cur * Ly_cur))
        fys.append(fy / (Lz_cur * Lx_cur))
        '''

plt.plot(displ, fzs)
#plt.plot(displ, fys)
#plt.plot(displ, fxs)
plt.show()
