# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 18:38:22 2020

@author: kkurman1
"""
import numpy as np

class Boundary: 
    def __init__(self, textline):
        words = textline.strip().split()
        self.textline = textline
        self.pbcX = True if words[-3] == "pp" else False
        self.pbcY = True if words[-2] == "pp" else False
        self.pbcZ = True if words[-1] == "pp" else False
    def set_x_bounds(self, xlo, xhi):
        self.xlo, self.xhi = xlo, xhi
        self.Lx = xhi - xlo
    
    def set_y_bounds(self, ylo, yhi):
        self.ylo, self.yhi = ylo, yhi
        self.Ly = yhi - ylo
    
    def set_z_bounds(self, zlo, zhi):
        self.zlo, self.zhi = zlo, zhi
        self.Lz = zhi - zlo
    
    def get_bounds_text(self):
        text = self.textline
        text = text +  "%e %e\n%e %e\n%e %e\n" %(self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi)
        return text
        
class BoxCell:
    def __init__(self, cell_id):
        self.cell_id = cell_id
        self.elements = []
    def add_element(self, elem):
        self.elements.append(elem)


def wrap_coordinates(atom_forces, bounds):
    for af in atom_forces:
        if bounds.pbcX:
            if af.x > bounds.xhi: af.x -= bounds.Lx
            if af.x < bounds.xlo: af.x += bounds.Lx
        if bounds.pbcY:
            if af.y > bounds.yhi: af.y -= bounds.Ly
            if af.y < bounds.ylo: af.y += bounds.Ly
        if bounds.pbcZ:
            if af.z > bounds.zhi: af.z -= bounds.Lz
            if af.z < bounds.zlo: af.z += bounds.Lz

def construct_cell_list(atom_forces, bounds, rc):
    #numnber of cells in x, y and z dimensions
    Nx, Ny, Nz = int(bounds.Lx//rc), int(bounds.Ly//rc), int(bounds.Lz//rc)
    print(bounds.get_bounds_text())
    print(bounds.Lx, bounds.Ly, bounds.Lz)
    N = Nx * Ny * Nz #total number of cells
    rcx, rcy, rcz = bounds.Lx/Nx, bounds.Ly/Ny, bounds.Lz/Nz
    #create cells
    cells = []
    print(Nx, Ny, Nz, N, rcx, rcy, rcz)
    for i in range(N): cells.append(BoxCell(i))
    wrap_coordinates(atom_forces, bounds)
    for i in range(len(atom_forces)):
        af = atom_forces[i]
        #make all coordinates positive
        x = (af.x - bounds.xlo)
        y = (af.y - bounds.ylo)
        z = (af.z - bounds.zlo)
        #get cell indices in x, y, z directions
        cx, cy, cz = x//rcx, y//rcy, z//rcz
        if cx < 0 or cy < 0 or cz < 0: print("Error, cx: %d cy: %d cz: %d" %(cx, cy, cz))
        #linearize cell index to get index in an array
        c_idx = int( cz + Nz*(cy + Ny*cx) )
        cells[c_idx].elements.append(i)
    #cells stores indices of atom_forces in an array
    return cells

def get_displ_pbr(x_next, x_prev, L):
        #get displacement vector pointing form x_prev to x_next for periodic boundary condition
        #periodicity is L
        sign = np.sign(x_next - x_prev)
        dx = abs(x_next - x_prev)
        if dx > L/2:
            dx = -(L - dx)
        return sign * dx

def get_dxdydz(a_next, a_cur, bounds):
    dx = a_next.x - a_cur.x
    dy = a_next.y - a_cur.y
    dz = a_next.z - a_cur.z
    if bounds.pbcX:  dx = get_displ_pbr(a_next.x, a_cur.x, bounds.Lx)
    if bounds.pbcY:  dy = get_displ_pbr(a_next.y, a_cur.y, bounds.Ly)
    if bounds.pbcZ:  dz = get_displ_pbr(a_next.z, a_cur.z, bounds.Lz)
    return dx, dy, dz


def create_neighbor_lists(atom_forces, bounds, rc):
    cells = construct_cell_list(atom_forces, bounds, rc)
    Nx, Ny, Nz = int(bounds.Lx//rc), int(bounds.Ly//rc), int(bounds.Lz//rc)
    #N = Nx * Ny * Nz #total number of cells
    #scan cells
    rc_sq = rc**2
    pair_count = 0
    for cx in range(Nx):
        for cy in range(Ny):
            for cz in range(Nz):
                c_idx = cz + Nz * (cy + Ny * cx)
                mid_cell = cells[c_idx]
                #scan neighbor cells including itself
                for ix in range(cx-1, cx+2):
                    for iy in range(cy-1, cy+2):
                        for iz in range(cz-1, cz+2):
                            #shift cell index to range [0, N-1]
                            if not bounds.pbcZ and (iz < 0 or iz > Nz): continue
                            nbr_c_idx = (iz + Nz)%Nz + ((iy + Ny)%Ny) * Nz + ((ix + Nx)%Nx) * Nz * Ny
                            #print("mid: %d nbr: %d" %(c_idx, nbr_c_idx))
                            nbr_cell = cells[nbr_c_idx]
                            for af_idx in mid_cell.elements:
                                for nbr_af_idx in nbr_cell.elements:
                                    af     = atom_forces[af_idx]
                                    nbr_af = atom_forces[nbr_af_idx]
                                    if af.id < nbr_af.id:
                                        dx, dy, dz = get_dxdydz(af, nbr_af, bounds)
                                        #dx = get_displ_pbr(af.x, nbr_af.x, bounds.Lx)
                                        #dy = get_displ_pbr(af.y, nbr_af.y, bounds.Ly)
                                        #dz = get_displ_pbr(af.z, nbr_af.z, bounds.Lz)
                                        rsq = dx**2 + dy**2 + dz**2
                                        if rsq < rc_sq:
                                            af.neighbors.append(nbr_af)
                                            nbr_af.neighbors.append(af)
                                            pair_count += 1
                                        if rsq > rc_sq * 15:
                                            print("mid: %d nbr: %d" %(c_idx, nbr_c_idx))
                                            print(af.x, af.y, af.z, nbr_af.x, nbr_af.y, nbr_af.z, rsq**(1/2))
    print("Pair count: %d" %pair_count)
    del cells
    return pair_count
                                            
def test_neighbor_lists(atom_forces, pair_ids, atype, bounds, rc):
    pair_count = create_neighbor_lists(atom_forces, bounds, rc)
    M = len(atom_forces)
    mismatch_count = 0
    for (id1, id2) in pair_ids:
        if id1 > M or id2 > M: continue
        af1 = atom_forces[id1-1]
        af2 = atom_forces[id2-1]
        pair_count -= 1
        miss = True
        for nbr_af in af1.neighbors:
            if nbr_af.id == af2.id:
                miss = False
                break
        if miss:
            dx, dy, dz = get_dxdydz(af1, af2, bounds)
            #dx = get_displ_pbr(af1.x, af2.x, bounds.Lx)
            #dy = get_displ_pbr(af1.y, af2.y, bounds.Ly)
            #dz = get_displ_pbr(af1.z, af2.z, bounds.Lz)
            print((dx**2 + dy**2 +dz**2)**(1/2))
            print(af1.x, af1.y, af1.z, af2.x, af2.y, af2.z)
            mismatch_count += 1
            #print("Mismatch %d %d" %(id1, id2))
        miss = True
        for nbr_af in af2.neighbors:
            if nbr_af.id == af1.id:
                miss = False
                break
        #if miss: 
        #    print("Mismatch")
    print("# of mismatches: %d total: %d # pairs lammps-cell list (diff): %d" %(mismatch_count, len(pair_ids), pair_count))
'''
b = Boundary("ITEM: BOX BOUNDS pp pp ff")
b.set_x_bounds(1,2)
b.set_y_bounds(4,5)
b.set_z_bounds(6,7)
print(b.get_bounds_text(),b.pbcX,b.pbcY, b.pbcZ)
'''       