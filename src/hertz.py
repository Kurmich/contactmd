#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 23:16:30 2020

@author: kurmich
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from dumpparser import *
import argparse
from analyzer import FileNames, ConesimSettings
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
epsilon = 0.000001
forces = [500.0]
poisson = 0.5
G_shear_mod = 16.0
E     = 27 # youngs modulus
E_star = E/(1 - poisson**2) #4 * G_shear_mod
R = 25
sigma = 1.0 #2**(1/6)
d = 2.0**(1.0/6.0) * sigma
atom_N = [1, 2, 3, 4]






def binary_search_hi(atom_forces, max_r):
    """Returns the index of atom located at radius closest to max_r (i.e. atom.radius <= max_r)"""
    lo = 0
    hi = len(atom_forces)
    while hi - lo > 1:
        mid = (hi + lo)//2
        if atom_forces[mid].radius <= max_r:
            lo = mid
        else:
            hi = mid
    return lo


def binary_search_lo(atom_forces, min_r):
    """Returns the index of atom located at radius closest to max_r (i.e. atom.radius <= max_r)"""
    lo = -1
    hi = len(atom_forces)-1
    while hi - lo > 1:
        mid = (hi + lo)//2
        if atom_forces[mid].radius >= min_r:
            hi = mid
        else:
            lo = mid
    return hi


def get_total_fxfyfz(atom_forces, ilo = 0, ihi = -1):
    if ihi == -1: ihi = len(atom_forces)
    fx, fy, fz = 0, 0, 0
    for i in range(ilo, ihi):
        af  = atom_forces[i] 
        fx += af.fx
        fy += af.fy
        fz += af.fz
    return fx, fy, fz
'''
def get_avg_pressure(atom_forces, r):
    count = len(atom_forces)
    normal_pressures = np.zeros(count)
    i = 0
    for atom_f in atom_forces:
        normal_pressures[i] = atom_f.fz/(math.pi*((d/2)**2))
        i += 1
    if len(normal_pressures) == 0:
        print("NO ATOMS TO COMPUTE AVG PRESSURE")
        return 0
    total_pressure = np.sum(normal_pressures)
    print("Radius r: %f total pressure: %f number of atoms: %d" %(r, total_pressure, count))
    return abs(total_pressure/count) #IS ABS VALUE ALWAYS VALID?
'''

def get_hertz_pressures(atom_forces):
    _,_,fz_tot = get_total_fxfyfz(atom_forces)
    r_max = ( 3 * fz_tot * R / (4*E_star) )**(1/3)   #contact radius
    p_max = 3 * fz_tot / (2 * math.pi * r_max**2 ) #maximum pressure
    
    bins = np.arange(0, r_max, r_max/100)
    bins = np.append(bins, r_max)
    
    p_max_red = p_max / E_star
    pres_red  = [p_max_red * math.sqrt(1 - (r/r_max)**2) for r in bins]
    return bins, pres_red


def get_hertz_penetration(atom_forces):
    _,_,fz_tot = get_total_fxfyfz(atom_forces)
    r_max = ( 3 * fz_tot * R / (4*E_star) )**(1/3)   #contact radius
    return r_max**2 / R

def get_avg_pressures(atom_forces):
    """Assumption: atom_forces are sorted by radial distance from z-axis"""
    
    max_r = atom_forces[-1].radius
    r_limit = max_r #* math.cos(math.pi/4)
    bins = np.arange(sigma/2, r_limit + sigma/4, sigma )    
    #calculate average pressures
    avg_pres = []
    for bin_r in bins:
        rlo     = bin_r - sigma/2
        rhi     = bin_r + sigma/2
        ilo     = binary_search_lo(atom_forces, rlo)
        ihi     = binary_search_hi(atom_forces, rhi)
        _,_, fz = get_total_fxfyfz(atom_forces, ilo, ihi)
        area    = 2 * math.pi * bin_r * sigma
        pzz     = fz/area 
        #append average 
        avg_pres.append(pzz)
        #print(bin_r, low_idx, up_idx)
    reduced_avg_pres = [p/E_star for p in avg_pres]
    
    return bins, reduced_avg_pres
   
    


def plot_avg_pressure(atom_forces):
    bins, reduced_avg_pres = get_avg_pressures(atom_forces)
    hertz_bins, hertz_pres = get_hertz_pressures(atom_forces)
    _,_,fz_tot = get_total_fxfyfz(atom_forces)
    xs, ys = [], []
    for af in atom_forces:
        xs.append(af.x)
        ys.append(af.y)
    #plt.scatter(xs, ys)
    plt.show()
    plt.plot(hertz_bins, hertz_pres, 'r')
    plt.plot(bins, reduced_avg_pres, 'bo')
    plt.suptitle('Normal pressure distribution in contact zone.\n Normal load: %d' %(fz_tot), fontsize = 16) #adhoc here on choice of forces
    plt.ylabel(r'$p/E^{\ast}$', fontsize = 16)
    plt.xlabel(r'$r/\sigma$', fontsize = 16)
    plt.show()

def select_from_frame(frame, type):
    '''Get given type of molecules from set of all molecules'''
    idx = frame[:, 1] == type
    return frame[idx, :]


def print_total_load(frame, type):
    frame = select_from_frame(frame, type)
    print(np.sum(frame[:, 7]))
    
def plot_nforce_vs_cont_area():
    max_radii = []
    t1 = np.arange(0.0, 0.06, 0.01)
    for force in forces:
        filename = 'visualize_%d.out' %force
        res = get_interactions(filename)
        max_r = 0
        for type, atom_forces in res.items():
            max_r = max(max_r, atom_forces[-1].radius)
        max_radii.append(max_r)
    new_forces = [(3*force/(4*E_star*R**2))**(1/3) for force in forces]
    new_radii = [r/R for r in max_radii]
    print(t1)
    plt.plot(t1, t1, 'r')
    plt.plot(new_forces, new_radii, 'bo')
    plt.suptitle('Contact radius a vs normal load N', fontsize = 20)
    plt.ylabel(r'$a/R$', fontsize = 16)
    plt.xlabel(r'$(3N/4E^{\ast}R^2)^{1/3}$', fontsize = 16)
    plt.show()
    
fig, ax = plt.subplots()
line1, = ax.plot([], [], 'bo')
ax.set(xlim=(0, 50), ylim=(0, 0.6))
ax.set_ylabel(r'$p/E^{\ast}$', fontsize = 16)
ax.set_xlabel(r'$r/\sigma$', fontsize = 16)
depth_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
line2, = ax.plot([], [], 'r', label = 'Hertz')
def init():
    line1.set_data([],[])
    line2.set_data([],[])
    depth_text.set_text('')
    return line1, line2, depth_text

def animate(i):
    bins, reduced_avg_pres = get_avg_pressures(anim_atom_forces[i])
    hertz_bins, hertz_pres = get_hertz_pressures(anim_atom_forces[i])
    line1.set_data(bins, reduced_avg_pres)
    line2.set_data(hertz_bins, hertz_pres)
    d_h = get_hertz_penetration(anim_atom_forces[i])
    d = (times[i] - times[0]) * css.vz * css.dt
    depth_text.set_text(r'$d$ = %.1f, $d_h$ = %.1f' % (d, d_h))
    #ax.set_title("Penetration d = %g" %d)
    return line1, line2, depth_text

anim_atom_forces = []  
times = None  
def animate_hertz_contact():
    global anim_atom_forces, times
    t_init, t_final = 5, 75
    tip_type = 2
    glass = 1
    types = [tip_type]
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)
    all_res, bounds, times = get_interactions(filenames.vis, t_init, t_final, types, interacting = True, sortkey = radius_sort)
    t0 = times[0]
    print(t0)
    for j in range(len(all_res)):
        res = all_res[j]
        #(times[j] - t0)
        atom_forces = res[tip_type]
        anim_atom_forces.append(atom_forces)
    print(len(anim_atom_forces))
    anim = FuncAnimation(fig, animate, init_func=init, interval=150, frames=len(anim_atom_forces)-1, repeat_delay=3000, blit=True)
    plt.legend()
    plt.draw()
    plt.show()
    anim.save('hertz_pressure_M%d_N%d_T%g_r%d.mp4' %(css.M, css.N, css.T, css.r), dpi=200, writer=writer)
    
filenames, css = None, None
radius_sort = idsort = lambda af: af.radius
def main():
    parser = argparse.ArgumentParser(description = "Contact analysis")
    parser.add_argument('--M',    type=int,   default = 2000,   help = '# of chains in a melt')
    parser.add_argument('--N',    type=int,   default = 256,    help = '# of monomers per chain')
    parser.add_argument('--T',    type=float, default = 0.0001, help = 'Temperature of the system')
    parser.add_argument('--poissons_r',  type=float,   default = 0.5,   help = 'Poissons ratio of the substrate')
    parser.add_argument('--E',    type=float,   default = 256,    help = 'Youngs modulus')
    parser.add_argument('--E_star',    type=float, default = 0.0001, help = 'Contact modulus')
    parser.add_argument('--vz',   type=float, default = 0.0001, help = 'Indenting tip velocity')
    parser.add_argument('--dt',   type=float, default = 0.01,   help = 'Integration time step of the simluation')
    parser.add_argument('--r',    type=float, default = 80,     help = 'Reduced radius of curvature')
    parser.add_argument('--cang', type=float, default = 45,     help = 'Angle the cone surface makes with the horizontal plane')
    parser.add_argument('--stiff',     action = 'store_true', default = False, help = 'True if polymer stiff (i.e. there is an angle style defined)')
    parser.add_argument('--conetip',   action = 'store_true', default = False, help = 'True if tip is of spherical shape')
    args = parser.parse_args()
    Temp = 0.0001
    args.T = Temp
    is_stiff = True
    global css, filenames
    args.r = 25
    global R
    R = args.r
    print(args.M, args.N, args.T, args.r, args.cang, args.stiff, args.conetip)
    css = ConesimSettings(args.M, args.N, args.T, args.r, args.cang, args.vz, args.dt)
    #css.set_analysisvals(1, 50, 1)
    args.stiff = True
    filenames = FileNames(args.M, args.N, args.T, args.r, args.cang, args.stiff, args.conetip)
    dt = css.dt
    #cang = 45
    tip_type = 2
    glass = 1
    types = [tip_type]
    atype = glass
    #rc = 1.5
    vz = css.vz
    d0 = 0 #2.2 
    t_init, t_final = 46, 47
    #t_step = css.t_step
    filename = filenames.vis
    animate_hertz_contact()
    return
    all_res, bounds, times = get_interactions(filename, t_init, t_final, types, interacting = True, sortkey = radius_sort)
    atom_forces = all_res[0][tip_type]
    plot_avg_pressure(atom_forces)
    
    


if __name__=="__main__":
    main()