#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:48:22 2020

@author: kurmich
"""


import math
import numpy as np
import matplotlib.pyplot as plt
from dumpparser import *
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from scipy.interpolate import griddata
from boxcell import *

class BinStatistics:
    def __init__(self):
        self.all_bins = []
        self.all_drs  = []
        self.all_dfs  = []
        self.all_dens = []
        self.all_ts   = []
    def test(self):
        return len(self.all_bins) == len(self.all_ts)

def plot_layer_density(atoms):
    atoms = filter_by_height(atoms, -20, 20)
    N = len(atoms)
    z_vals = np.zeros(N)
    for i in range(N):
        z_vals[i] = atoms[i].z
    hist, bin_edges = np.histogram(z_vals, bins = 120)
    bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])
    plt.close()
    plt.plot(bincenters, hist, '-', marker = 'o', fillstyle = 'none', markerfacecolor = 'r')
    plt.xlabel(r'$z(\sigma)$')
    plt.title("Number of atoms vs z.")
    
    #ax.hist(z_vals, bins='auto')
    plt.legend()
    plt.show()
    
    
def filter_by_height(atoms, zlo, zhi):
    filtered_atoms = []
    for atom in atoms:
        if atom.z > zlo and atom.z < zhi:
            filtered_atoms.append(atom)
    return filtered_atoms



def save_fluctuations(binned_stats, filename):
    N = len(binned_stats.all_bins)
    with open(filename, 'w') as f:
        for i in range(N):
            f.write("TIME: %d\n" %binned_stats.all_ts[i])
            f.write("z dr df density\n")
            z_bins, drs      = binned_stats.all_bins[i], binned_stats.all_drs[i] 
            dfs, n_densities = binned_stats.all_dfs[i],    binned_stats.all_dens[i]
            for j in range(len(z_bins)):
                f.write("%g %g %g %g\n" %(z_bins[j], drs[j], dfs[j], n_densities[j]))



def read_fluctuations(filename):
    z_bins, drs      = [], []
    dfs, n_densities = [], []
    binned_stats = BinStatistics()
    with open(filename, 'r') as f:
        for line in f:
            if "TIME:" in line:
                if len(z_bins) != 0:
                    binned_stats.all_bins.append(z_bins)
                    binned_stats.all_drs.append(drs)
                    binned_stats.all_dfs.append(dfs)
                    binned_stats.all_dens.append(n_densities)
                words = line.split(' ')
                binned_stats.all_ts.append(float(words[-1]))
                z_bins, drs      = [], []
                dfs, n_densities = [], []
            elif "z" not in line:
                words = line.split(' ')
                z_bins.append(float(words[0]))
                drs.append(float(words[1]))
                dfs.append(float(words[2]))
                n_densities.append(float(words[3]))
                
    if len(z_bins) != 0:
        binned_stats.all_bins.append(z_bins)
        binned_stats.all_drs.append(drs)
        binned_stats.all_dfs.append(dfs)
        binned_stats.all_dens.append(n_densities)
    assert binned_stats.test()
    return binned_stats

def visualize_fluctuations(css, filenames, del_z):
    glass = 1
    types = [glass]
    atype = glass
    #rc = 1.5
    vz = css.vz
    t_init, t_final = css.t_init, css.t_final
    t_step = css.t_step
    step = 1 #step for calculations
    filename             = filenames.vis
    flucatuations_filename  = "fluctuations_M%d_N%d_T%g_r%d_cang%d.out" %(css.M, css.N, css.T, css.r, css.cang)
    #remove files if they already exist
    #remove_file(vischanges_filename)

    binned_stats = BinStatistics()
    for t_start in range(t_init, t_final, t_step):
        t_end = t_start + t_step + step - 1
        all_res, all_bounds, times         = get_interactions(filename, t_start, t_end, types, interacting = False)
        set_force_position_changes(all_res, binned_stats, atype, all_bounds, del_z, step)
        binned_stats.all_ts.extend(times)
        print("t start: %d\n" %t_start, flush = True)
    save_fluctuations(binned_stats, flucatuations_filename)
    

memory_atom_forces_p  = None
def set_force_position_changes(all_res, binned_stats, atype, all_bounds, del_z, step):
    '''ASSUMPTION: that all atom_forces are sorted by their IDs'''
    drs_batch       = []
    dfs_batch       = []
    z_bins_batch    = []
    n_dens_batch    = []
    global memory_atom_forces_p
    N = len(all_res)
    assert N >= step
    for i in range(step, N):
        atom_forces                 = all_res[i][atype]
        atom_forces_p               = all_res[i-step][atype]
        z_idx_tuples, bin_indices   = get_indices_and_pos(atom_forces_p, del_z)
        z_bins, drs, dfs            = get_pos_force_fluctuations(atom_forces, atom_forces_p, z_idx_tuples, bin_indices)
        total                       = bin_indices[-1] - bin_indices[0] + 1
        n_density                   = [(bin_indices[i] - bin_indices[i-1] + 1)/total for i in range(1, len(bin_indices))]
        z_bins_batch.append(z_bins)
        drs_batch.append(drs)
        dfs_batch.append(dfs)
        n_dens_batch.append(n_density)
        if i == step and memory_atom_forces_p is not None:
            print("Assigning from memory")
            atom_forces_p = memory_atom_forces_p
        elif i == N-1:
            print("Assigning memory")
            memory_atom_forces_p = atom_forces
    binned_stats.all_bins.extend(z_bins_batch)
    binned_stats.all_drs.extend(drs_batch)
    binned_stats.all_dfs.extend(dfs_batch)
    binned_stats.all_dens.extend(n_dens_batch)


def get_indices_and_pos(atom_forces, del_z):
    '''Returns: 1)  z-index tuples sorted w.r.t. z positions. 
                    Here index is index in atom_forces array where the atom at position z is located.
                2)  Indices for bins separated by del_z.
    '''
    N            = len(atom_forces)
    z_idx_tuples = sorted([(atom_forces[k].z, k) for k in range(N)], key = lambda tup: tup[0]) #array of tuples (z position, index in the array)
    bin_indices  = [0]
    zlo ,_       = z_idx_tuples[0]
    zhi          = zlo + del_z
    for i in range(N):
        z,_ = z_idx_tuples[i]
        if z < zhi: continue
        bin_indices.append(i)
        zlo = zhi
        zhi = zlo + del_z
    if bin_indices[-1] != N-1: bin_indices.append(N-1)
    return z_idx_tuples, bin_indices

def get_pos_force_fluctuations(atom_forces, atom_forces_p, z_idx_tuples, bin_indices):        
    N = len(bin_indices)
    drs, dfs, z_bins = np.zeros(N-1), np.zeros(N-1), np.zeros(N-1)
    for i in range(1, N):
        ilo, ihi = bin_indices[i-1], bin_indices[i]
        avg_dr = 0
        avg_df = 0
        for j in range(ilo, ihi+1):
            _, idx =  z_idx_tuples[j]
            af     = atom_forces[idx]
            af_p   = atom_forces_p[idx]
            assert af.id == af_p.id 
            avg_dr += ( (af.x - af_p.x)**2 + (af.y - af_p.y)**2 + (af.z - af_p.z)**2 )**(1/2)
            avg_df += ( (af.fx - af_p.fx)**2 + (af.fy - af_p.fy)**2 + (af.fz - af_p.fz)**2 )**(1/2)
        avg_dr /= (ihi-ilo+1)
        avg_df /= (ihi-ilo+1)
        #print(z_idx_tuples[ilo][0], avg_dr, avg_df)
        drs[i-1], dfs[i-1], z_bins[i-1] = avg_dr, avg_df, z_idx_tuples[ilo][0]
    return z_bins, drs, dfs




def get_avg_points_and_vals(atom_forces, bounds, rc):
    cells = construct_cell_list_2D(atom_forces, bounds, rc)
    xs, ys, values = [], [], []
    for cell in cells:
        x_sum, y_sum, val_sum = 0, 0, 0
        N = len(cell.elements)
        if N == 0: continue
        for idx in cell.elements:
            x_sum   += atom_forces[idx].x
            y_sum   += atom_forces[idx].y
            val_sum += atom_forces[idx].fz
        xs.append(x_sum/N)
        ys.append(y_sum/N)
        values.append(val_sum)
    return np.array( [ xs, ys ] ).T, values
            
def autocorrfz(atom_forces, bounds, zlo, zhi):
    plt.close()  
    print(len(atom_forces))
    atom_forces         = filter_by_height(atom_forces, zlo, zhi)
    print(len(atom_forces))
    values              = [ af.fz for af in atom_forces ]
    fz_av               = np.mean(values)
    fz_sum              = np.sum(values)
    print("avg sum:" , fz_av, fz_sum, np.max(values), np.min(values))
    
    points              = np.array( [ [af.x for af in atom_forces], [af.y for af in atom_forces] ] ).T
    print( points.shape )
    points, values      = get_avg_points_and_vals(atom_forces, bounds, 1.5)
    fz_av               = np.mean(values)
    fz_sum              = np.sum(values)
    print("avg sum:" ,fz_av, fz_sum)
    nx                  = 2 * int(bounds.xhi - bounds.xlo)
    ny                  = 2 * int(bounds.yhi - bounds.ylo)
    xs                  = np.linspace(bounds.xlo - 0.5, bounds.xhi + 0.5, nx)
    ys                  = np.linspace(bounds.ylo - 0.5, bounds.yhi + 0.5, ny)
    grid_x, grid_y      = np.meshgrid(xs, ys)
    grid_fz             = griddata(points, values, (grid_x, grid_y), method='linear', fill_value = fz_av)
    #print(np.isnan(grid_fz).sum(), np.isnan(values).sum())
    force_fft           = np.fft.fft2(grid_fz)
    pow_spec            = force_fft * np.conj(force_fft) / len(force_fft)
    force_autocorr      = np.fft.ifft2(pow_spec)
    #print(force_fft)
   # print(pow_spec)
    #print(force_autocorr)
    #plt.imshow(grid_fz.T, extent=(bounds.xlo - 0.5, bounds.xhi + 0.5, bounds.ylo - 0.5, bounds.yhi + 0.5))
    #plt.imshow(np.abs(np.fft.ifftshift(pow_spec_mult)).T, extent=(bounds.xlo - 0.5, bounds.xhi + 0.5, bounds.ylo - 0.5, bounds.yhi + 0.5))
    plt.imshow(np.abs(np.fft.ifftshift(force_autocorr)).T, extent=(bounds.xlo - 0.5, bounds.xhi + 0.5, bounds.ylo - 0.5, bounds.yhi + 0.5))
    #plt.plot(points[:,0], points[:,1], 'k.', ms=1)
    plt.colorbar()
    plt.show()

        

fig, ax = plt.subplots()
line1, = ax.plot([], [], 'bo', label = "Average Displacements")
ax.set(xlim=(-40, 40), ylim=(0, 0.5))
#ax.set_ylabel(r'$stats$', fontsize = 16)
ax.set_xlabel(r'$z(\sigma)$', fontsize = 16)
depth_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
line2, = ax.plot([], [], 'ro', label = 'Number densities')
def init():
    line1.set_data([],[])
    line2.set_data([],[])
    depth_text.set_text('')
    return line1, line2, depth_text

def animate(i):
    z_bins, drs, n_densities = binstats.all_bins[i], binstats.all_dfs[i], binstats.all_dens[i]
    time = binstats.all_ts[i]
    line1.set_data(z_bins, drs)
    line2.set_data(z_bins, n_densities)
    depth_text.set_text(r'$t$ = %d' % (time))
    #ax.set_title("Penetration d = %g" %d)
    return line1, line2, depth_text

anim_atom_forces = []  
times = None  
def animate_fluctuations():
    global anim_atom_forces, times
    t_init, t_final = 5, 75
    tip_type = 2
    glass = 1
    types = [tip_type]
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)
    anim = FuncAnimation(fig, animate, init_func=init, interval=200, frames=len(binstats.all_dens)-1, repeat_delay=3000, blit=True)
    plt.legend()
    plt.draw()
    plt.show()
    #anim.save('hertz_pressure_M%d_N%d_T%g_r%d.mp4' %(css.M, css.N, css.T, css.r), dpi=200, writer=writer)
    

#binstats = read_fluctuations("fluctuations_M2000_N256_T0.0001_r0_cang0.out")
#animate_fluctuations()