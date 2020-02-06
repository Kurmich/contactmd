import matplotlib.pyplot as plt
import math

M, N = 2000, 256
T = 0.0001
R, cang = 10, 45
filename = "../visfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve_nzT_stats.txt" %(M, N, T, R, cang)
filename = "../visfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve.txt" %(M, N, T, R, cang)

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




def heat_stats(filename, scale):
    N = 100000
    max_time = 1000000000000
    dt = 0.01
    vz = 0.0001
    times, fzs, fys, pes  = [], [], [], []
    del_N = 20000
    with open(filename, 'r') as file:
        for line in file:
            words = line.strip().split()
            if len(words) == 0: continue
            if words[0] == "Step":
                N = len(words)
                keywords = words
                print(line)
            if N == len(words) and words[0].isdigit():
                time = int(words[0])
                if time > max_time: break
                times.append(time)
                fzs.append(float(words[-1]))
                fys.append(float(words[-2]))
                pes.append(float(words[-4]))
    print(N, len(times), times[-1], times[0])
    t0 = times[0]
    del_t = times[1] - times[0]
    heat = [ fzs[i]-fzs[i-1] for i in range(1, len(fzs))]
    dws  = [fzs[i] * vz *dt * del_t   for i in range(0, len(fzs))]
    cum_works = []
    cum_work = 0
    for dw in dws:
        cum_work += dw
        cum_works.append(cum_work)
    ds, qs = [], []
    us, ws = [], []
    print(len(cum_works))
    for i in range(2*del_N, len(cum_works), del_N):
        ds.append(times[i-del_N] * vz * dt)
        qs.append( -scale * ((pes[i] - pes[i-del_N]) -  (cum_works[i] - cum_works[i-del_N])) )
        us.append( scale * (pes[i] - pes[i-del_N]) )
        ws.append(scale*(cum_works[i] - cum_works[i-del_N]))
    newtimes = [t - t0 for t in times]
    plt.plot(ds, qs, label = "-Q", color = 'black')
    #plt.plot(ds, us, label = "$\Delta U$")
    #plt.plot(ds, ws, label = "W")
    plt.xlabel("d")
    plt.ylabel("$-Q, Frac$")
    #plt.ylabel("$Q, W, \Delta U$")
    plt.legend()
    plt.show() 


def motion_stats(filename):
    N = 100000
    max_time = 14000000000
    with open(filename, 'r') as file:
        for line in file:
            words = line.strip().split()
            if len(words) == 0: continue
            if words[0] == "Step":
                N = len(words)
                keywords = words
                times, x2s, x2s_ave, vacfs  = [], [], [], []
                print(line)
            if N == len(words) and words[0].isdigit():
                time = int(words[0])
                if time > max_time: break
                times.append(time)
                x2s_ave.append(float(words[-1]))
                x2s.append(float(words[-2]))
                vacfs.append(float(words[-3]))
    t0 = times[0]
    print(t0)
    newtimes = [t - t0 for t in times]
    tt = [math.log10(newtimes[i]) for i in range(1, len(newtimes))]
    x2s_log = [math.log10(x2s[i]) for i in range(1, len(x2s))]
    x2save_log = [math.log10(x2s_ave[i]) for i in range(1, len(x2s_ave))]
    #plt.plot(newtimes, x2s, label = "x^2")
    #plt.plot(newtimes, x2s_ave, label = "xave^2")
    plt.plot(tt, x2s_log, label = "$x^2$")
    plt.plot(tt, x2save_log, label = "$x_{ave}^2$")
    plt.plot(tt, vacfs[:-1], label = "$v_{acf}$")
    plt.xlabel("$log(t)$")
    plt.ylabel("$log(<x^2>)$")
    plt.legend()
    plt.show()
    
def main():
    M, N = 2000, 256
    T = 0.0001
    R, cang = 10, 45
    #filename = "../visfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve_nzT_stats.txt" %(M, N, T, R, cang)
    #motion_stats(filename)
    filename = "../visfiles//conetip_M%d_N%d_T%g_sphR%d_cang%d_nve.txt" %(M, N, T, R, cang)
    heat_stats(filename, 0.00000416)
    
    
if __name__ == "__main__":
    main()
'''
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
        ''''''
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
'''
plt.plot(displ, fzs)
#plt.plot(displ, fys)
#plt.plot(displ, fxs)
plt.show()

'''


