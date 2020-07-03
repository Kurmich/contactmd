import matplotlib.pyplot as plt
import math
import copy
import numpy             as np

M, N = 2000, 256
T = 0.0001
R, cang = 10, 45
dt = 0.01
filename = "../visfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve_nzT_stats.txt" %(M, N, T, R, cang)
filename = "../visfiles/conetip_M%d_N%d_T%g_sphR%d_cang%d_nve.txt" %(M, N, T, R, cang)

params = {}
params[(1000,256)] = [72.1, 72.1, 52.8]
params[(2000,256)] = [90.84, 90.84, 66.52]
params[(2000,500)] = [113.93, 113.93, 82]
vz = 0.0001
dt = 0.01
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
    return ds, qs, us, ws
    for i in range(len(qs)):
        qs[i] += 0.0007
    
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
    t_strainz       = [math.log(lz/Lzs[0]) for lz in Lzs]
    pois_strains, pois_ratios = get_poissons_ratios(Lxs, Lys, Lzs)
    E = get_Youngs_modulus(t_strainz, pzz)
    plt.plot(t_strainz, pzz, label = "$T = %g, E = %3.1f$" %(T,-E))
    plt.plot( pois_strains, pois_ratios, label ="T = %g, $ \\nu $" %T)
    #plt.plot(Lzs, pyy, label = "$p_{yy}$")
    #plt.plot(Lzs, pxx, label = "$p_{xx}$")
    plt.xlabel("$\epsilon_t$", fontsize=16)
    plt.ylabel("$\sigma_t$", fontsize=16, rotation=0)
    plt.legend()


def get_Yield_stress(strains, stresses):
    return 0


def get_poissons_ratios(Lxs, Lys, Lzs):
    step            = 2
    t_strainz       = [math.log(Lzs[i]/Lzs[i-step]) for i in range(step, len(Lzs))]
    t_strainx       = [math.log(Lxs[i]/Lxs[i-step]) for i in range(step, len(Lxs))]
    poissons_ratios = [-t_strainx[i]/t_strainz[i]   for i in range(len(t_strainz))]
    strains         = [math.log(Lzs[i]/Lzs[0])      for i in range(step, len(Lzs))]
    #get average poissons rations for respective strains by averaging Nevery value
    Nevery = 100
    avg_strains, avg_poissons_ratios = [], []
    for i in range(Nevery, len(poissons_ratios), Nevery):
        total_pr = 0
        for j in range(i-Nevery, i):
            total_pr += poissons_ratios[i]
        avg_strains.append(strains[i-Nevery])
        avg_poissons_ratios.append(total_pr/Nevery)
    return avg_strains, avg_poissons_ratios

def get_Youngs_modulus(strains, stresses):
    idx        = 1
    max_strain = 0.03 #small enough
    for eps in strains:
        if abs(eps) < max_strain:
            idx += 1
        else:
            break
    #fit linear equation
    (E,b) = np.polyfit(strains[0:idx], stresses[0:idx], 1)
    return E

def autocorr_stats(filename):
    N = 100000
    max_time = 14000000000
    times, x2s, x2s_ave, vacfs  = [], [], [], []
    cur_times, cur_x2s, cur_x2s_ave, cur_vacfs  = [], [], [], []
    count = 0
    with open(filename, 'r') as file:
        for line in file:
            words = line.strip().split()
            if len(words) == 0: continue
            if words[0] == "Step":
                N = len(words)
                keywords = words
                if len(cur_times) > 1000:
                    t0 = cur_times[0]
                    cur_times = [t - t0 for t in cur_times]
                    count += 1
                    if len(times) == 0:
                        times   = copy.deepcopy(cur_times)
                        x2s     = copy.deepcopy(cur_x2s)
                        x2s_ave = copy.deepcopy(cur_x2s_ave)
                        vacfs   = copy.deepcopy(cur_vacfs)
                    else:
                        assert len(times) == len(cur_times), "Lengths of lists must match"
                        times     = [ times[i]   + cur_times[i]         for i in range(len(times))]
                        x2s       = [ x2s[i]     + cur_x2s[i]           for i in range(len(x2s))]
                        x2s_ave   = [ x2s_ave[i] + cur_x2s_ave[i]       for i in range(len(x2s_ave))]
                        vacfs     = [ vacfs[i]   + cur_vacfs[i]         for i in range(len(vacfs))]
                cur_times, cur_x2s, cur_x2s_ave, cur_vacfs  = [], [], [], []
                print(line)
            if N == len(words) and words[0].isdigit():
                time = int(words[0])
                if time > max_time: break
                cur_times.append(time)
                cur_x2s_ave.append(float(words[-1]))
                cur_x2s.append(float(words[-2]))
                cur_vacfs.append(float(words[-3]))
    t0 = times[0]
    print(t0)
    print(count)
    print(len(times), len(x2s))
    times     = [ times[i]   /count       for i in range(len(times))]
    x2s       = [ x2s[i]     /count       for i in range(len(x2s))]
    x2s_ave   = [ x2s_ave[i] /count       for i in range(len(x2s_ave))]
    vacfs     = [ vacfs[i]   /count       for i in range(len(vacfs))]
    tt = [math.log10(times[i]) for i in range(1, len(times))]
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
    contactd = 2.8 + 2.5
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
    return ds, comp_frac, ext_frac
    plot_changes(ds, comp_frac, ext_frac, contactd, delta_r)
    

def impose_heat_bonds(heatfile, ljbondfile):
    ds1, qs, us, ws = heat_stats(heatfile, 1, 1)
    ds, comp_frac, ext_frac = bond_change_stats(ljbondfile)
    for i in range(len(ds1)):
        if ds1[i] == ds[0]:
            qs, us, ws = qs[i:], us[i:], ws[i:]
            break
    N = min(len(qs), len(comp_frac))
    cfrac = [comp_frac[i] + ext_frac[i] for i in range(N)]
    coeffs = np.polyfit(cfrac, qs[:N], 1)
    scale, shift = coeffs[0], coeffs[1]
    print(scale, shift)
    cfrac = [scale * cf + shift for cf in cfrac]
    d0 = ds[0]
    ds = [ds[i]-d0 for i in range(len(ds))]
    plt.plot(ds[:N], qs[:N], label = '-Q', color = 'black')
    plt.plot(ds[:N], cfrac,  label = 'Scaled changes', color = 'green')
    plt.xlabel("d")
    plt.ylabel("Fraction, -Q")
    plt.legend()
    plt.show()



def get_scale_shift(l1, l2):
    '''Compute scaling and shifting needed for two vectors to have 'best' fit.
    Ref: https://stackoverflow.com/questions/13563453/finding-the-best-scale-shift-between-two-vectors'''
    N      = min(len(l1), len(l2))
    l1     = np.array(l1[:N])
    l2     = np.array(l2[:N])
    l1     = np.abs(l1) / np.sum(l1)
    l2     = np.abs(l2) / np.sum(l2)
    m1, v1 = np.mean(l1), np.var(l1)
    m2, v2 = np.mean(l2), np.var(l2)
    scale  = (v2/v1)**(1/2)
    shift  = m2 - scale * m1
    return scale, shift

def main():
    M, N = 2000, 256
    T = 0.1
    Ts = [0.1]
    drs = [0.3, 0.6, 0.9, 1.2]
    R, cang = 25, 45
    for T in Ts:
        filename = "../outputfiles/deform_M%d_N%d_T%g.txt" %(M, N, T)
        compression_stats(filename, T)
    plt.show()
    '''filename = "../outputfiles/autocorr_stiff_M%d_N%d_T%g.txt" %(M, N, T)
    autocorr_stats(filename)'''
    '''heat_filename = "../outputfiles/spheretip_stiff_M%d_N%d_T%g_sphR%d_nve.txt" %(M, N, T, R)
    bond_filename = "../outputfiles/stats_stiff_M%d_N%d_T%g_r%d_p%g.out" %(M, N, T, R, 0.3)
    impose_heat_bonds(heat_filename, bond_filename)'''
    '''for T in Ts:
        for dr in drs:
            bond_filename = "../outputfiles/stats_stiff_M%d_N%d_T%g_r%d_p%g.out" %(M, N, T, R, dr)
            ds, comp_frac, ext_frac = bond_change_stats(bond_filename)
            c_frac = [comp_frac[i] + ext_frac[i] for i in range(len(comp_frac))]
            ds = [d - ds[0] for d in ds]
            plt.plot(ds, c_frac, label = "$T = %g, \Delta R = %g$" %(T, dr))
    plt.xlabel("d")
    plt.ylabel("Fraction")
    plt.legend()
    plt.show()'''
    #heat_stats(heat_filename, 0.8214 * 0.0000045 * 2, 1)
    #bond_change_stats(bond_filename)
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


