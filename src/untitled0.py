# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 11:55:46 2019

@author: kurmanbek
"""

import numpy as np
import matplotlib.pyplot as plt
import math


N = 100
R = 3

beta_min = -math.pi / 4
beta_max =  math.pi / 4
betas = np.linspace(beta_min, beta_max, 500, endpoint = True)

theta_min = 0
theta_max = math.pi / 2
thetas = np.linspace(theta_min, theta_max, 1000)

rs = np.linspace(0, 0.99*R, 1000)


mu_n = 1
I = 1
dI = I / N
coef = mu_n * I / (2 * math.pi) # this actually doesn't matter since we are taking the ratio in the end

alpha_l  = lambda r, theta, beta: math.atan2( (R * math.sin(beta) - r * math.sin(theta)), (R * math.cos(beta) + r * math.cos(theta)) )

alpha_r  = lambda r, theta, beta: math.atan2( (R * math.sin(beta) - r * math.sin(theta)), (R * math.cos(beta) - r * math.cos(theta)) )

a      = lambda r, theta, beta: math.sqrt( (R * math.sin(beta) - r * math.sin(theta))**2 + (R * math.cos(beta) - r * math.cos(theta))**2 )

b      = lambda r, theta, beta: math.sqrt( (R * math.sin(beta) - r * math.sin(theta))**2 + (R * math.cos(beta) + r * math.cos(theta))**2 )

def dB_l(r, theta, beta):
    alpha = alpha_l(r, theta, beta)
    dB_lx =  -  coef *  math.sin(alpha) / b(r, theta, beta)  
    dB_ly = -  coef * math.cos(alpha)/ b(r, theta, beta) 
    return dB_lx, dB_ly

def dB_r(r, theta, beta):
    alpha = alpha_r(r, theta, beta)
    dB_rx =   coef *  math.sin(alpha)/ a(r, theta, beta) 
    dB_ry = - coef * math.cos(alpha) / a(r, theta, beta) 
    return dB_rx, dB_ry


Bx2, By2 = 0, 0
count = 0
for r in rs:
    for theta in thetas:
        Bx, By = 0, 0
        for beta in betas:
            dB_lx, dB_ly = dB_l(r, theta, beta)
            dB_rx, dB_ry = dB_r(r, theta, beta)
            Bx += (dB_rx + dB_lx)
            By += (dB_ry + dB_ly)
        #print(Bx, By)
        Bx2 += (Bx**2) 
        By2 += (By**2)

#since we are taking ratio all cofficients do not matter

print(Bx2/By2) # = 0.01389203151903227
