# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:36:11 2020

@author: kkurman1
"""

import matplotlib.pyplot as plt
import numpy as np
T = 0.0001
file2 = "../rate_T%g_deld0.4.txt" %T
file1 = "../rate_T%g_deld0.2.txt" %T 
file3 = "../rate_T%g_deld0.6.txt" %T   


displacements = [0.2, 0.4, 0.6]
'''
N = len(displacements)
for i in range(1):
    d1 = displacements[i]
    file1 = "../rate_T%g_deld%g.txt" %(T, d1)
    r1 = []
    with open(file1, "r") as f:
            for line in f:
                r1.append(float(line.strip()))
    for j in range(i+1, N):
        d2 = displacements[j]
        
        
        ratio = int(d2/d1)
        r1_new = []
        for k in range(len(r1) - (ratio-1)):
            tmp = 0
            for kk in range(ratio):
                tmp += r1[kk]
            r1_new.append(tmp/ratio)
        
        file2 = "../rate_T%g_deld%g.txt" %(T, d2)
        r2 = []
        with open(file2, "r") as f:
            for line in f:
                r2.append(float(line.strip()))
        coeff = np.polyfit(r1_new, r2[:-(ratio-1)], 1)
        poly1d_fn = np.poly1d(coeff)
        plt.scatter(r1_new, r2[:-(ratio-1)])
        print(len(r2), len(r1_new))
        plt.plot( r1_new, poly1d_fn(r1_new), 'k--', label = "m = %5.4f, b = %5.4f" %(coeff[0], coeff[1]))
        plt.legend()
        plt.show()
    






r3 = []
with open(file3, "r") as f:
  for line in f:
    r3.append(float(line.strip()))


#plt.plot(r1, label = "d = 0.2")
#plt.plot(r2, label = "d = 0.4")

r1_new = []
for i in range(len(r1) - 1):
    r1_new.append((r1[i] + r1[i+1])/2)
coeff = np.polyfit(r1_new, r2[:-1], 1)
poly1d_fn = np.poly1d(coeff)
plt.scatter(r1_new, r2[:-1])
print(len(r2), len(r1_new))
plt.plot( r1_new, poly1d_fn(r1_new), 'k--', label = "m = %5.4f, b = %5.4f" %(coeff[0], coeff[1]))
plt.legend()
plt.show()
'''

r1 = []
with open(file1, "r") as f:
  for line in f:
    r1.append(float(line.strip()))



r2 = []
with open(file2, "r") as f:
  for line in f:
    r2.append(float(line.strip()))
 
r3 = []
with open(file3, "r") as f:
  for line in f:
    r3.append(float(line.strip()))


#plt.plot(r1, label = "d = 0.2")
#plt.plot(r2, label = "d = 0.4")

r1_new = []
x = 0
for i in range(len(r1) - 1):
    r1_new.append((r1[i] + r1[i+1])/2)
coeff = np.polyfit(r1_new, r2[:-1], 1)
poly1d_fn = np.poly1d(coeff)
plt.scatter(r1_new[x:], r2[x:-1])
print(len(r2), len(r1_new))
plt.plot( r1_new[x:], poly1d_fn(r1_new[x:]), 'k--', label = "m = %3.2f, b = %3.2f, displ = 0.4" %(coeff[0], coeff[1]), color = 'green')
plt.legend()
plt.show()


r1_new2 = []
for i in range(len(r1) - 2):
    r1_new2.append((r1[i] + r1[i+1] + r1[i+2])/3)
coeff = np.polyfit(r1_new2, r3[:-2], 1)
poly1d_fn = np.poly1d(coeff)
plt.scatter(r1_new2, r3[:-2])
print(len(r3), len(r1_new2))
plt.plot( r1_new2, poly1d_fn(r1_new2), 'k--', label = "m = %3.2f, b = %3.2f, displ = 0.6" %(coeff[0], coeff[1]))
plt.legend()
plt.show()