#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 17:19:51 2020

@author: kurmich
"""

import secrets
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import math
import numpy as np
import settings

dt = settings.dt

def interacting_particles(atom_forces):
    interacting = []
    epsilon = settings.epsilon
    for af in atom_forces:
        if abs(af.fx) > epsilon or abs(af.fy) > epsilon or abs(af.fz) > epsilon:
            interacting.append(af)
    N = len(interacting)
    points = np.zeros([N, 2])
    for i in range(N):
        points[i, 0], points[i, 1] = interacting[i].x, interacting[i].y
    return interacting, points


def shoelace_area(xs, ys):
    '''Let 'vertices' be an array of N pairs (x,y), indexed from 0
        Let 'area' = 0.0
        for i = 0 to N-1, do
          Let j = (i+1) mod N
          Let area = area + vertices[i].x * vertices[j].y
          Let area = area - vertices[i].y * vertices[j].x
        end for'''
    N = len(xs)
    area = 0
    for i in range(N):
        j = (i+1) % N
        area += (xs[i] * ys[j] - ys[i] * xs[j])
    return abs(area)/2

def orientation(atom_f1, atom_f2, atom_f3):
    #http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    val = (atom_f2.y - atom_f1.y) * (atom_f3.x - atom_f2.x) - (atom_f3.y - atom_f2.y) * (atom_f2.x - atom_f1.x)
    if val == 0: #collinear
        return 0
    elif val > 0: #clockwise
        return 1
    else:
        return 2 # counterlockwise
def further(atom_f1, atom_f2, atom_f3):
    d1 = (atom_f2.y - atom_f1.y) * (atom_f2.y - atom_f1.y) + (atom_f2.x - atom_f1.x) * (atom_f2.x - atom_f1.x)
    d2 = (atom_f3.y - atom_f1.y) * (atom_f3.y - atom_f1.y) + (atom_f3.x - atom_f1.x) * (atom_f3.x - atom_f1.x)
    if d2 > d1:
        return 1 #return 1 if atom 2 is further than atom 1
    return -1

def jarvis(atom_forces, visualize = False):
    leftmost = atom_forces[0]
    for af in atom_forces:
        if af.x < leftmost.x:
            leftmost = af
    hull = [leftmost]
    max_r = 0
    cur, cand = leftmost, None
    xs, ys = [], []
    remaining = [af for af in atom_forces if af not in hull]
    while True:
        cand = secrets.choice(atom_forces)
        for af in atom_forces:
            if af != cand:
                if orientation(cur, cand, af) == 2:
                    #print("here")
                    cand = af
                elif orientation(cur, cand, af) == 0 and further(cur, cand, af) == 1:
                    cand = af
        if cand == leftmost: break
        cur = cand
        max_r = max(max_r, cur.radius)
        xs.append(cur.x)
        ys.append(cur.y)
        hull.append(cand)

    area = shoelace_area(xs, ys)

    if visualize:
        xs, ys = [], []
        for af in hull:
            xs.append(af.x)
            ys.append(af.y)
        print("vis area %g" %(shoelace_area(xs, ys)))
        xs.append(hull[0].x)
        ys.append(hull[0].y)
        plt.plot(xs, ys, color = 'r')
        xs, ys = [], []
        for af in atom_forces:
            xs.append(af.x)
            ys.append(af.y)
        plt.scatter(xs, ys)
        plt.show()
    return hull, max_r, area
#def get_edge_avg(atom_forces):

def get_total_forces(atomic_forces):
    fx, fy, fz = 0, 0, 0
    for af in atomic_forces:
        fx += af.fx
        fy += af.fy
        fz += af.fz
    print(fx, fy, fz)
    return fx, fy, fz


def get_contact_depth(interacting_afs):
    if len(interacting_afs) == 0: return 0
    zlo, zhi = 100000000000, -100000000000
    for af in interacting_afs:
        zlo = min(zlo, af.z)
        zhi = max(zhi, af.z)
    return zhi - zlo


def get_indentation_params(all_res, times, v, type):
    dt         = settings.dt
    areas      = []
    ts         = []
    stresses_z = []
    fzs        = []
    num_intrs  = []
    t0         = times[0]
    times = [ t - t0 for t in times]
    assert len(all_res) == len(times),"%d %d" %(len(all_res), len(times))
    first_contact = True
    for i in range(len(all_res)):
        t = times[i]
        print("New timestep: %d" %t)
        res = all_res[i]
        interacting_af, points = interacting_particles(res[type])
        count = len(interacting_af)
        if count < 3:
            fx, fy, fz = get_total_forces(interacting_af)
            area = count * math.pi * (1.12/2)**2
            hc = get_contact_depth(interacting_af)
            if first_contact:
                print("First contact time: ", t)
                first_contact = False
        else:
            if first_contact:
                print("First contact time: ", t)
                first_contact = False
            fx, fy, fz = get_total_forces(interacting_af)
            hull = ConvexHull(points)
            #plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
            #plt.show()
            print("scipy area: %g" %hull.area)
            area = shoelace_area(points[hull.vertices,0], points[hull.vertices,1])
            hc = get_contact_depth(interacting_af)
            #hull, max_r, area = jarvis(interacting_af, visualize=True)
            #area = math.pi * max_r**2
        d = v * dt * t
        if area != 0: stress_zz = fz / area
        else:         continue  
        
        areas.append(area)
        fzs.append(fz)
        num_intrs.append(count)
        ts.append(t)
        stresses_z.append(stress_zz)
        print("Displacement d: %g Contact Depth: %g Num of particles: %d Total Fz: %g Area: %g StressZ: %g " %(d, hc, count, fz, area, stress_zz))
    return ts, num_intrs, areas, fzs, stresses_z

def plot_stresszz_d(all_res, times, v, type):
    dt         = settings.dt
    ts, num_intrs, areas, fzs, stresses_z = get_indentation_params(all_res, times, v, type)
    ds = [v * t * dt for t in ts]
    lj_times = [ t*dt for t in times ]
    fig, ax = plt.subplots(2, 2)
    ax[0][0].plot(ds, num_intrs, 'g')
    ax[0][0].set_ylabel("Atoms in contact", rotation = 90)
    ax[0][1].plot(ds, areas, 'r')
    ax[0][1].set_ylabel("$A(a^2)$", rotation = 90)
    ax[1][0].plot(ds, fzs, 'k')
    ax[1][0].set_ylabel("$F_z(u_0/a)$", rotation = 90)
    ax[1][1].plot(ds, stresses_z, 'b')
    ax[1][1].set_ylabel("$\sigma_{zz}(u_0/a^3)$", rotation = 90)
    for i in range(2):
        for j in range(2):
            if True:
                 ax[i][j].set_xlabel(r'$t(t_{LJ})$')
                 #ax[i][j].set_xscale('log')
            else:
                ax[i][j].set_xlabel("$d/a$")
           
    fig.suptitle("Indentation stats")
    plt.show()
    fig2, ax2 = plt.subplots(1, 1)
    tshift = 10**7 + 2* 10**5
    idx = 0
    for i in range(len(ts)):
        if ts[i] - tshift >= 0:
            idx = i
            break
    ds = ds[i:]
    ds = [d-ds[0] for d in ds]
    ax2.plot(ds, num_intrs[i:])
    ax2.set_xscale('log')
    plt.show()
    
    