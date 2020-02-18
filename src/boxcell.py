# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 18:38:22 2020

@author: kkurman1
"""
class Boundary: 
        
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
        text = "ITEM: BOX BOUNDS pp pp ff\n%e %e\n%e %e\n%e %e\n" %(self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi)
        return text
        
class BoxCell:
    def __init__(self, cell_id):
        self.cell_id = cell_id
        self.elements = []
    def add_element(self, elem):
        self.elements.append(elem)

       