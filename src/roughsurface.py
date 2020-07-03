#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 13:37:05 2020

@author: kurmich
"""
import numpy as np
import math
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
    def __init__(self, L_vec, N_vec, spectrum):
        self.Lx, self.Ly = L_vec
        self.Nx, self.Ny = N_vec
        self.spectrum    = spectrum
        self.__make_psd()
        
    def __make_psd(self):
        qL_sq  = self.spectrum.qL**2
        qs_sq = self.spectrum.qs**2
        qr_sq = self.spectrum.qr**2
        freq_x   = np.fft.fftfreq(self.Nx, 1/self.Nx) * 2 * np.pi / self.Lx
        freq_y   = np.fft.fftfreq(self.Ny, 1/self.Ny) * 2 * np.pi / self.Ly
        xx, yy = np.meshgrid(freq_x, freq_y)
        q_sq_mat   = xx**2 + yy**2
        self.psd = np.zeros((self.Nx, self.Ny))
        for i in range(self.Nx):
            for j in range(self.Ny):
                q_sq = q_sq_mat[i, j]
                print(i, j , q_sq)
                if q_sq < qL_sq  or  q_sq > qs_sq:
                    continue
                if q_sq >= qL_sq and q_sq < qr_sq:
                    self.psd[i, j] = self.spectrum.C
                else:
                    self.psd[i, j] = self.spectrum.C / (q_sq/qr_sq)**(self.spectrum.H + 1)
        print(self.psd)
        self.q_sq_mat = q_sq_mat
        print(xx.shape, yy.shape, freq_x.size, freq_y.size)
        
    def random_phase_surface(self):
        phases = self.get_rand_uniform_phases()
        H = np.exp(phases*1j) * np.sqrt(self.psd)
        self.heights = (self.Lx**2 + self.Ly**2)**(1/2) * np.fft.ifft2(H)
        print(self.heights)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.imshow(phases)
        ax2.imshow(np.fft.fftshift(phases))
        plt.show()
    def get_rand_uniform_phases(self):
        M, N = self.psd.shape
        phases_centered = np.zeros((M,N))
        M_half = int(M/2)
        N_half = int(N/2)
        print(M_half, N_half)
        phases_centered[0:(M_half+1), :] =  np.random.uniform(0, 2*np.pi, (M_half + 1, N) )
        X = M - M_half - 1
        phases_centered[M_half+1:M, :] = -np.flip(np.flip(phases_centered[0:X, :], axis = 0), axis=1)
        lo, hi = 0, N-1
        while hi >= lo:
            phases_centered[M_half, lo] = - phases_centered[M_half, hi]
            lo += 1
            hi -= 1
            if lo == hi: phases_centered[M_half, lo] = 0
        return np.fft.ifftshift(phases_centered)
    def vis_surface(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        Z = self.heights.real
        Nx, Ny = Z.shape
        x = np.arange(0, Nx, 1)
        y = np.arange(0, Ny, 1)
        X, Y = np.meshgrid(x, y)
        
        print(Z.shape, X.shape, Y.shape)
        ax.plot_surface(X, Y, Z,cmap=cm.jet)
        

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        plt.show()
        
                    
    def vis_psd(self):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        cs, qs = [], []
        Dr, Dc = self.psd.shape
        for r in range(Dr):
            for c in range(Dc):
                if self.psd[r, c] == 0: continue
                qs.append( self.q_sq_mat[r,c] )
                cs.append(self.psd[r, c])
        ax1.imshow(np.fft.fftshift(self.psd))
        #ax1.colorbar()
        ax2.scatter(qs, cs)
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        plt.show()

a = 2**(1/6)/2
L = 11
qs = 2 * np.pi / a
qL = 2 * np.pi / L
print(qs, qL)
spectrum   = Spectrum(qL, qL*2, qs, 0.7, 1)   
L_vec = (11,11) 
N_vec = (11,11)  
r = RoughSurface(L_vec, N_vec, spectrum)
r.random_phase_surface()
r.vis_surface()
#r.vis_psd()