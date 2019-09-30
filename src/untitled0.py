# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 11:55:46 2019

@author: kurmanbek
"""

import numpy as np
import matplotlib.pyplot as plt
import math

B = 1
q = 1
Q = 1
m = 1
c = 1
k = 1 #8.9875517873681764 * 10**9
R0 = 3

C1 =  k * q * Q / m
C2 =   q * B / (m*c)
ang = 0
ang_vel = 0
ang_acc = 0
acc = -C1 / (R0 * R0)
pos = R0
vel = 0
dt = 0.01
times = np.arange(0, 147, dt)
rs = []
thetas = []
a = (2 * q * Q/ (m * R0)) *(2 * m *c / (q*B))**2
tmp = 3 * (3 * math.sqrt(3 * (4*a**3 - 13*a**2 + 32*a)) + 9 * a + 16)**(1/3)
rmin = R0 * ( tmp / (3 * 2**(1/3)) - 1/3 - 2**(1/3) * (3 * a - 4) / tmp) /3


for t in times:
    new_pos = pos + vel * dt + acc*dt*dt * 0.5
    new_ang = ang + ang_vel * dt + ang_acc*ang_acc*dt*dt * 0.5
    
    
    new_ang_vel = 0.5 * C2 * (1 - (R0/new_pos)**2)
    
    new_acc = - C1 / (new_pos * new_pos) -  C2 * new_ang_vel * pos
    
    new_vel = vel + (acc + new_acc) * dt * 0.5
    new_ang_acc =  C2 * R0*R0 * new_vel / (new_pos**3)
    
    ang = new_ang
    ang_vel = new_ang_vel
    ang_acc = new_ang_acc
    
    pos = new_pos
    acc = new_acc
    vel = new_vel
    
    rs.append(pos)
    thetas.append(ang)

xs = [rs[i] * math.cos(thetas[i]) / R0 for i in range(len(times)) ]
ys = [rs[i] * math.sin(thetas[i]) / R0 for i in range(len(times)) ]
plt.plot( xs, ys, label = 'Trajectory')

xs = [rmin * math.cos(thetas[i]) / R0 for i in range(len(times)) ]
ys = [rmin * math.sin(thetas[i]) / R0 for i in range(len(times)) ]
plt.plot( xs, ys, label = r'$R_{min}$')
plt.xlabel(r'$x/R_0$')
plt.ylabel(r'$y/R_0$')
plt.legend()
plt.show()