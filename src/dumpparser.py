#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 21:42:57 2020

@author: kurmich
"""

from boxcell import *
from atom import AtomicForces
epsilon = 0.000001

idsort = lambda af: af.id

def get_interactions(filename, t_start, t_end, types, interacting = False, sortkey = idsort):
    assert t_start <= t_end
    r_time = False
    r_atom_count = False
    r_boundary = False
    r_atoms = False
    skip = False
    dim = 0
    t = 0
    max_r = 0
    pbcX, pbcY, pbcZ = False, False, False #boundary periodicity flags
    all_res = []
    all_bounds = []
    times = []
    cur_b = []
    res = dict()
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                times.append(time)
                #print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                res = dict()
                #print("# of atoms: %d" %count)
                r_atom_count = False
            elif r_boundary:
                if dim > 0:
                    words = line.strip().split()
                    lo, hi = map(float, words)
                    if   dim == 3: bounds.set_x_bounds(lo, hi)
                    elif dim == 2: bounds.set_y_bounds(lo, hi)
                    elif dim == 1: bounds.set_z_bounds(lo, hi)
                dim -= 1
                if dim == 0:
                    all_bounds.append(bounds)
                    r_boundary = False
            elif r_atoms:
                if 'ITEM: TIMESTEP' in line:
                    r_atoms = False
                    r_time = True
                else:

                    if skip:
                        continue
                    line = line.strip()
                    '''Return a list of the words in the string, using sep as the delimiter string. 
                    If maxsplit is given, at most maxsplit splits are done (thus, the list will have at most maxsplit+1 elements). 
                    If maxsplit is not specified or -1, then there is no limit on the number of splits (all possible splits are made).'''
                    words = line.split(' ', 9)
                    if len(words) > 9:
                        a_id, mol, atype, x, y, z, fx, fy, fz, arributes  = words
                    else:
                        a_id, mol, atype, x, y, z, fx, fy, fz  = words
                    a_id, mol, atype, x, y, z, fx, fy, fz = int(a_id), int(mol), int(atype), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    if atype not in types: continue
                    if interacting and abs(fx) < epsilon and abs(fy) < epsilon and abs(fz) < epsilon: continue
                    #choose interacting atoms
                    radius = (x**2 + y**2)**(1/2)
                    if atype not in res: res[atype] = []
                    res[atype].append(AtomicForces(a_id, mol, atype, x, y, z, fx, fy, fz))
                    max_r = max(max_r, radius)
                    
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if len(res) != 0:
                    for type, atom_forces in res.items():
                        l = sorted(atom_forces, key = sortkey)
                        res[type] = l
                    all_res.append(res)
                elif len(times) > 0:
                    times.pop()
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                #print("Next is timestep")
            elif 'ITEM: NUMBER OF ATOMS' in line:
                r_atom_count = True
                #print("Next is number of atoms")
            elif 'ITEM: BOX BOUNDS' in line:
                words = line.strip().split()
                pbcX = True if words[-3] == "pp" else False
                pbcY = True if words[-2] == "pp" else False
                pbcZ = True if words[-1] == "pp" else False
                r_boundary = True
                bounds = Boundary(line)
                dim = 3
                #print("Next 3 lines are bondaries")
            elif 'ITEM: ATOMS' in line:
                r_atoms = True
                atr_line = line[line.index('S')+1:]
                atr_words = atr_line.split(' ')
                #print("Atom coordinate lines are coming")
    return all_res, all_bounds, times


def get_frames(filename, t_start, t_end, types = None):
    assert t_start <= t_end
    r_time = False
    r_atom_count = False
    r_boundary = False
    r_atoms = False
    skip = False
    d = 0
    t = 0
    max_r = 0
    cur_frame = None
    real_count = 0
    frames = []
    #type_count = {}
    #for type in types:
     #   type_count[type] = 0
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                if not skip:
                    cur_frame = np.zeros([count, 8])
                print("# of atoms: %d" %count)
                r_atom_count = False
            elif r_boundary:
                d -= 1
                r_boundary = False if d == 0 else True
            elif r_atoms:
                if 'ITEM: TIMESTEP' in line:
                    r_atoms = False
                else:
                    #print("reading atoms")
                   # print(len(line.split(' ')))
                    if skip:
                        continue
                        
                    id, mol, type, x, y, z, fx, fy, fz, _  = line.split(' ')
                    id, mol, type, x, y, z, fx, fy, fz = int(id), int(mol), int(type), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    #if types is not None and type in types and time == t_start:
                    #    typecount[type] += 1
                    #if type not in types: continue
                    cur_frame[id-1, 0], cur_frame[id-1, 1] = mol, type
                    cur_frame[id-1, 2], cur_frame[id-1, 3], cur_frame[id-1, 4] = x, y, z
                    cur_frame[id-1, 5], cur_frame[id-1, 6], cur_frame[id-1, 7] = fx, fy, fz
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if cur_frame is not None and not skip:
                    frames.append(cur_frame)
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                print("Next is timestep")
            elif 'ITEM: NUMBER OF ATOMS' in line:
                r_atom_count = True
                print("Next is number of atoms")
            elif 'ITEM: BOX BOUNDS pp pp mm' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ATOMS' in line:
                r_atoms = True
                print("Atom coordinate lines are coming")
    return  frames


def get_pair_interactions(filename, t_start, t_end):
    assert t_start <= t_end
    r_time = False
    r_entry_count = False
    r_boundary = False
    r_entries = False
    skip = False
    d = 0
    t = 0
    pair_counts = []
    all_pair_ids = []
    res = []
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_entry_count:
                count = int(line)
                pair_counts.append(count)
                res = []
                print("# of pair ids: %d" %count)
                r_entry_count = False
            elif r_boundary:
                d -= 1
                r_boundary = False if d == 0 else True
            elif r_entries:
                if 'ITEM: TIMESTEP' in line:
                    r_entries = False
                else:
                    #print("reading atoms")
                   # print(len(line.split(' ')))
                    if skip:
                        continue
                    line = line.strip()
                    id1, id2  = line.split(' ')
                    id1, id2 = int(id1), int(id2)
                    res.append((id1, id2))
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if len(res) != 0:
                    all_pair_ids.append(res)
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                print("Next is timestep")
            elif 'ITEM: NUMBER OF ENTRIES' in line:
                r_entry_count = True
                print("Next is number of entries")
            elif 'ITEM: BOX BOUNDS' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ENTRIES' in line:
                r_entries = True
                print("Pair id lines are coming")
    return all_pair_ids, pair_counts


