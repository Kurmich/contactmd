# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 18:38:22 2020

@author: kkurman1
"""

class BoxCell:
    def __init__(self, cell_id):
        self.cell_id = cell_id
        self.elements = []
    def add_element(self, elem):
        self.elements.append(elem)
        
