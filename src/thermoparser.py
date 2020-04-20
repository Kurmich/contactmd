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




def running_average(arr, Nevery, Nrepeat, period):
    avgs = []
    for i in range(period, len(arr)):
        idx = i
        count = 1
        cur_sum = 0
        while count <= Nrepeat:
            cur_sum += arr[idx]
            idx     -= Nevery
            count   += 1
        avgs.append(cur_sum / Nrepeat)
        
    for i in range(period, len(arr)):
        arr[i] = avgs[i-period]

def heat_stats(filename, scale, step):
    N = 100000
    max_time = 1000000000000
    dt = 0.01
    vz = 0.0001
    times, fzs, fys, pes  = [], [], [], []
    del_N = step * 20000
    Nevery = 200
    Nrepeat = 100
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
    Nevery = int(Nevery/del_t)
    running_average(fzs, Nevery, Nrepeat, del_N)
    running_average(pes, Nevery, Nrepeat, del_N)
    #running_average(arr, Nevery, Nrepeat, del_N)
    
    #heat = [ fzs[i]-fzs[i-1] for i in range(1, len(fzs))]
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
    for i in range(len(qs)):
        qs[i] += 0.005
    
    plt.plot(ds, qs, label = "-Q", color = 'black')
    #plt.plot(ds, us, label = "$\Delta U$")
    #plt.plot(ds, ws, label = "W")
    plt.xlabel("d")
    plt.ylabel("$-Q, Frac$")
    #plt.ylabel("$Q, W, \Delta U$")
    plt.legend()
    


def compression_stats(filename, T):
    N = 100000
    max_time = 14000000000
    izhi = 13
    ipzz = 16
    Lxs, Lys, Lzs = [], [], []
    pxx, pyy, pzz = [], [], []
    with open(filename, 'r') as file:
        for line in file:
            words = line.strip().split()
            if len(words) == 0: continue
            if words[0] == "Step":
                N = len(words)
                assert words[izhi] == "Zhi"
                assert words[ipzz] == "Pzz"
                #sprint(line)
            if N == len(words) and words[0].isdigit():
                time = int(words[0])
                if time > max_time: break
                vals = list(map(float, words))
                Lz = vals[izhi]   - vals[izhi-1]
                Ly = vals[izhi-2] - vals[izhi-3]
                Lx = vals[izhi-4] - vals[izhi-5]
                if vals[ipzz] < 0:
                    Lxs, Lys, Lzs = [], [], []
                    pxx, pyy, pzz = [], [], []
                if Lz < 15: continue
                Lzs.append(Lz)
                Lys.append(Ly)
                Lxs.append(Lx)
                pzz.append(vals[ipzz])
                pyy.append(vals[ipzz-1])
                pxx.append(vals[ipzz-2])
    t_strainz = [math.log(lz/Lzs[0]) for lz in Lzs]
    E = get_Youngs_modulus(t_strainz, pzz)
    print(E)
    plt.plot(t_strainz, pzz, label = "$T = %g$" %T)
    #plt.plot(Lzs, pyy, label = "$p_{yy}$")
    #plt.plot(Lzs, pxx, label = "$p_{xx}$")
    plt.xlabel("$\epsilon_t$", fontsize=16)
    plt.ylabel("$\sigma_t$", fontsize=16, rotation=0)
    plt.legend()


def get_Yield_stress(strains, stresses):
    return 0

def get_Youngs_modulus(strains, stresses):
    p0 = stresses[0]
    e0 = strains[0]
    idx = 1
    for eps in strains:
        if abs(eps) < 0.04:  #ADHOC
            idx += 1
        else:
            break
    p1 = stresses[idx]
    e1 = strains[idx]
    E = (p1-p0) / (e1-e0)
    return E

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


def plot_changes(ds, ccfrac, cefrac, contactd, delta_r):
    cfrac = [ccfrac[i] + cefrac[i] for i in range(len(ccfrac))]
    #fig = plt.figure()
    plt.title("Changes by >= %g vs displacement of the tip" %(delta_r))
    plt.plot(ds, ccfrac, label = "Compressed")
    plt.plot(ds, cefrac, label = "Extended")
    plt.plot(ds, cfrac,  label = "Changed" )
    plt.axvline(x=contactd, color = 'red', label = "Contact point: %g" %contactd)
    plt.xlabel("d")
    plt.ylabel("Fraction")
    plt.legend()
    #plt.show()

def plot_breaks(ds, bfrac, ffrac, contactd):
    plt.title("Broken and formed bond fractions  vs displacement of the tip")
    plt.plot(ds, bfrac, label = "Broke")
    plt.plot(ds, ffrac, label = "Formed")
    plt.xlabel("d")
    plt.ylabel("Fraction")
    plt.axvline(x=contactd, color = 'red', label = "Contact point: %g" %contactd)
    plt.legend()
    #plt.show()

def bond_change_stats(filename):
    comp_frac, ext_frac, bfrac, ffrac  = [], [], [], []
    comp_ext_frac, ext_comp_frac = [], []
    ds = []
    vz = 0.0001
    dt = 0.01
    delta_r = 0.03
    contactd = 7
    with open(filename, 'r') as file:
        for line in file:
            words = line.strip().split()
            if words[0] == "time":
                print(line)
                continue
            t, comp, ext, breaks, forms, ce, ec, pair_count = map(int, words)
            ds.append(vz * (t*dt))
            comp_frac.append(comp / pair_count)
            ext_frac.append(ext / pair_count)
            bfrac.append(breaks / pair_count)
            ffrac.append(forms / pair_count)
            comp_ext_frac.append(ce / pair_count)
            ext_comp_frac.append(ec / pair_count)
    
    plot_changes(ds, comp_frac, ext_frac, contactd, delta_r)
    
   
def main():
    M, N = 2000, 256
    #T = 0.1
    Ts = [0.0001, 0.1, 0.2]
    R, cang = 10, 45
    for T in Ts:
        filename = "../outputfiles/deform_M%d_N%d_T%g.txt" %(M, N, T)
        compression_stats(filename, T)
    plt.show()
    #filename = "../visfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve_nzT_stats.txt" %(M, N, T, R, cang)
    #motion_stats(filename)
    #filename = "../outputfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve_nzT.txt" %(M, N, T, R, cang)
    #heat_stats(filename, 0.8214 * 0.00000416 * 4.5, 1)
    #filename = "../outputfiles/stats_M2000_N256_T0.2_r10_cang45_p0.3.txt"
    #bond_change_stats(filename)
    #plt.show()
    
    
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


