#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 13:37:05 2020

@author: kurmich
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

class Spectrum:
    def __init__(self, qL, qr, qs, H, C):
        self.qL     = qL #shortest wave vector (Longest wavelength)
        self.qr     = qr #roll-off wave vector (Intermediate wavelength)
        self.qs     = qs #longest wave vector  (shortest wavelength)
        self.H      = H  #hurst exponent
        self.C      = C
        
epsilon = 10**(-10)       
class RoughSurface:
    def __init__(self, Lx, Ly, spectrum):
        self.Lx, self.Ly = Lx, Ly
        self.spectrum    = spectrum
        self.__make_psd()
        
    def __make_psd(self):
        D   = self.spectrum.qs + 1
        self.psd = np.zeros((D,D))
        qL_sq = self.spectrum.qL**2
        qs_sq = self.spectrum.qs**2
        qr_sq = self.spectrum.qr**2
        for r in range(D):
            for c in range(D):
                q_sq = r**2 + c**2
                if q_sq < qL_sq  or  q_sq > qs_sq:
                    continue
                if q_sq >= qL_sq and q_sq < qr_sq:
                    self.psd[r, c] = self.spectrum.C
                else:
                    self.psd[r, c] = self.spectrum.C / (q_sq/qr_sq)**(self.spectrum.H + 1)
    def random_phase_surface(self):
        phases = self.get_rand_uniform_phases()
        H = np.exp(phases*1j) * np.sqrt(self.psd)
        self.heights = (self.Lx**2 + self.Ly**2)**(1/2) * np.fft.ifft2(H)
        print(self.heights)
        #plt.imshow(np.abs(self.heights))
        #plt.show()
    def get_rand_uniform_phases(self):
        M, N = self.psd.shape
        M_half = int(M/2) + 1
        N_half = int(N/2) + 1
        phases = np.random.uniform(0, 2*np.pi, self.psd.shape )
        #restrictions to make surface heights real
        for r in range(M_half):
            for c in range(N_half):
                phases[r, c] = -phases[M-r-1, N-c-1]
        #for r in range(M_half):
        #    phases[r, 0] = -phases[M-r-1, 0]
        #for c in range(N_half):
        #    phases[0, c] = -phases[0, N-c-1]
        phases[0, 0] = phases[0, N_half] = phases[M_half, 0] = phases[M_half, N_half] = 0
        print(phases)
        return phases
    def vis_surface(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        Z = self.heights.imag
        Nx, Ny = Z.shape
        x = np.arange(0, Nx, 1)
        y = np.arange(0, Ny, 1)
        X, Y = np.meshgrid(x, y)
        
        print(Z.shape, X.shape, Y.shape)
        ax.plot_surface(X, Y, Z,cmap=cm.jet)
        

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        
        plt.show()
        
                    
    def vis_psd(self):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        cs, qs = [], []
        Dr, Dc = self.psd.shape
        for r in range(Dr):
            for c in range(Dc):
                if self.psd[r, c] == 0: continue
                qs.append((r**2 + c**2)**(1/2))
                cs.append(self.psd[r, c])
        ax1.imshow(self.psd)
        #ax1.colorbar()
        ax2.scatter(qs, cs)
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        plt.show()

spectrum   = Spectrum(2, 6, 10, 0.7, 1)      
r = RoughSurface(10, 10, spectrum)
r.random_phase_surface()
r.vis_surface()
#r.vis_psd()