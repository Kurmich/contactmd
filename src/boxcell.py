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
        text = text +  "\n%e %e\n%e %e\n%e %e\n" %(self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi)
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
    Nx, Ny, Nz = bounds.Lx//rc, bounds.Ly//rc, bounds.Lz//rc
    N = Nx * Ny * Nz #total number of cells
    rcx, rcy, rcz = bounds.Lx/Nx, bounds.Ly/Ny, bounds.Lz/Nz
    #create cells
    cells = []
    for i in range(N): cells.append(BoxCell(i))
    wrap_coordinates(atom_forces, bounds)
    for af in atom_forces:
        cx, cy, cz = af.x//rcx, af.y//rcy, af.z//rcz
        c_idx = cz + Nz * cy + Nz * Ny * cx
        cells[c_idx].append(af)
    return cells

def get_displ_pbr(x_next, x_prev, L):
        #get displacement vector pointing form x_prev to x_next for periodic boundary condition
        #periodicity is L
        sign = np.sign(x_next - x_prev)
        dx = abs(x_next - x_prev)
        if dx > L/2:
            dx = -(L - dx)
        return sign * dx

def create_neighbor_lists(atom_forces, bounds, rc):
    cells = construct_cell_list(atom_forces, bounds, rc)
    Nx, Ny, Nz = bounds.Lx//rc, bounds.Ly//rc, bounds.Lz//rc
    #N = Nx * Ny * Nz #total number of cells
    #scan cells
    rc_sq = rc**2
    for cx in range(Nx):
        for cy in range(Ny):
            for cz in range(Nz):
                c_idx = cz + Nz * cy + Nz * Ny * cx
                mid_cell = cells[c_idx]
                #scan neighbor cells including itself
                for ix in range(cx-1, cx+2):
                    for iy in range(cy-1, cy+2):
                        for iz in range(cz-1, cz+2):
                            #shift cell index to range [0, N-1]
                            nbr_c_idx = (iz + Nz)%Nz + ((iy + Ny)%Ny) * Nz + ((ix + Nx)%Nx) * Nz * Ny
                            nbr_cell = cells[nbr_c_idx]
                            for af in mid_cell.elements:
                                for nbr_af in nbr_cell.elements:
                                    if af.id < nbr_af.id:
                                        dx = get_displ_pbr(af.x, nbr_af.x, bounds.Lx)
                                        dy = get_displ_pbr(af.y, nbr_af.y, bounds.Ly)
                                        dz = get_displ_pbr(af.z, nbr_af.z, bounds.Lz)
                                        rsq = dx**2 + dy**2 + dz**2
                                        if rsq < rc_sq:
                                            af.neighbors.append(nbr_af)
                                            nbr_af.append(af)
'''
b = Boundary("ITEM: BOX BOUNDS pp pp ff")
b.set_x_bounds(1,2)
b.set_y_bounds(4,5)
b.set_z_bounds(6,7)
print(b.get_bounds_text(),b.pbcX,b.pbcY, b.pbcZ)
'''       