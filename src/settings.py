#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 17:36:31 2020

@author: kurmic"""
a = 5
epsilon  = 0.000001
dt       = 0.01    #integration time step used in simulations
rc       = 1.89 + epsilon #cutoff distance for lennard Jones interactions
vis_data_path = '../visfiles/'
out_data_path = '../outputfiles/'
#idx_id, idx_mol, idx_type, idx_x, idx_y, idx_z, idx_fx, idx_fy, idx_fz, _ = line.split(' ')

#vis_data_path = '/home-3/kkurman1@jhu.edu/GlassExperiments/stiff/indentation/visfiles/'
#out_data_path =  '/home-3/kkurman1@jhu.edu/GlassExperiments/stiff/indentation/outputfiles/'

#vis_data_path = '/home/kkurman1/GlassExperiments/stiff/indentation/visfiles/'
#out_data_path =  '/home/kkurman1/GlassExperiments/stiff/indentation//outputfiles/'