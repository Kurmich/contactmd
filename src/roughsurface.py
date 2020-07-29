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
import argparse
from data import data

class Spectrum:
    def __init__(self, qL, qr, qs, H, C):
        assert qL <= qr and qr < qs, "qs > qr > qL is required. Currently qs: %g, qr: %g, qL: %g" %(qs, qr, qL)
        self.qL          = qL #shortest wave vector (Longest wavelength)
        self.qr          = qr #roll-off wave vector (Intermediate wavelength)
        self.qs          = qs #longest wave vector  (shortest wavelength)
        self.H           = H  #hurst exponent
        self.C           = C #/ qL**(2*(H+1))
        self.description = "qL: %g qr: %g qs: %g H: %g C: %g" %(qL, qr, qs, H, self.C)
    def get_info(self):
        return self.description
        
epsilon = 10**(-10)       
class RoughSurface:
    def __init__(self, L_vec, N_vec, spectrum):
        self.Lx, self.Ly = L_vec
        self.Nx, self.Ny = N_vec
        self.spectrum    = spectrum
        self.__make_psd()
        
    def __make_psd(self):
        self.rms_slope = 0
        qL_sq  = self.spectrum.qL**2
        qs_sq = self.spectrum.qs**2
        qr_sq = self.spectrum.qr**2
        freq_x   = np.fft.fftfreq(self.Nx, 1/self.Nx) #* 2 * np.pi / self.Lx
        freq_y   = np.fft.fftfreq(self.Ny, 1/self.Ny) #* 2 * np.pi / self.Ly
        xx, yy = np.meshgrid(freq_x, freq_y, indexing='ij')
        
        q_sq_mat   = xx**2 + yy**2
        self.psd = np.zeros((self.Nx, self.Ny))
        for i in range(self.Nx):
            for j in range(self.Ny):
                q_sq = q_sq_mat[i, j]
                if q_sq < qL_sq  or  q_sq > qs_sq:
                    continue
                if q_sq >= qL_sq and q_sq < qr_sq:
                    print("here")
                    self.psd[i, j] = self.spectrum.C
                else:
                    self.psd[i, j] = self.spectrum.C / (q_sq/qr_sq)**(self.spectrum.H + 1)
                self.rms_slope += q_sq * self.psd[i, j] * (4 * np.pi**2)
        #self.rms_slope /= ( self.Nx * self.Ny )
        self.roughness = np.sum(self.psd)#/( self.Nx * self.Ny )
        print("rms slope: %g rms roughness %g" %(np.sqrt(self.rms_slope), np.sqrt(self.roughness)) )
        print("psd max", np.max(self.psd))
        self.q_sq_mat = q_sq_mat
    
    def filtering_algorithm_surface(self):
        gauss_mat = self.get_rand_matrix(distribution = 'gaussian')
        H =  np.fft.fft2(gauss_mat) * np.sqrt(self.psd)
        self.heights =  np.fft.ifft2(H) *  (self.Nx * self.Ny)**(1/2) # (self.Nx**2 + self.Ny**2)**(1/4) 
        print("std: ", np.std(self.heights) )
    
    def random_phase_surface(self):
        phases = self.get_rand_matrix()
        H = np.exp(phases*1j) * np.sqrt(self.psd)
        self.heights = (self.Nx**2 + self.Ny**2)**(1/2) * np.fft.ifft2(H)
        
    def get_rand_matrix(self, distribution = 'uniform', visual = False):
        M, N = self.psd.shape
        phases_centered = np.zeros((M,N))
        M_half = int(M/2)
        #N_half = int(N/2)
        #print(M_half, N_half)
        if distribution == 'uniform':
            phases_centered[0:(M_half+1), :] =  np.random.uniform(0, 2*np.pi, (M_half + 1, N) )
        elif distribution == 'gaussian':
            mu, sigma = 0, 1
            phases_centered[0:(M_half+1), :] =  np.random.normal(mu, sigma, (M_half + 1, N) )
        #phases_centered[:, N_half] = -10
        X = M - M_half - 1
        lo, hi = 0, N-1
        if M%2 == 1:
            lo = 0
            if N%2 == 0:
                phases_centered[M_half+1:M, 1:] = -np.flip(np.flip(phases_centered[0:X, 1:], axis = 0), axis=1)
                phases_centered[:, 0] = 0
                lo = 1
            else:
                phases_centered[M_half+1:M, :]  = -np.flip(np.flip(phases_centered[0:X, :], axis = 0),  axis=1)
        else:
            lo = 1
            if N%2 == 0:
                phases_centered[M_half+1:M, 1:] = -np.flip(np.flip(phases_centered[1:X+1, 1:], axis = 0), axis=1)
                phases_centered[:, 0] = 0               
            else:
                phases_centered[M_half+1:M, :]  = -np.flip(np.flip(phases_centered[1:X+1, :], axis = 0),  axis=1)
                lo = 0
            phases_centered[0, :] = 0
        while hi >= lo:
                phases_centered[M_half, lo] = - phases_centered[M_half, hi]
                lo += 1
                hi -= 1
                if lo == hi: phases_centered[M_half, lo] = 0
            
        if visual:
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(phases_centered)
            ax2.imshow(np.fft.ifftshift(phases_centered))
            plt.show()
        return np.fft.ifftshift(phases_centered)
    def vis_surface(self):
        assert np.sum(self.heights.imag) < 0.01, "Generated surface heights have substantial imaginary part"
        
        Z = self.heights.real
        
        Nx, Ny = Z.shape
        x = np.arange(0, Nx, 1)
        y = np.arange(0, Ny, 1)
        X, Y = np.meshgrid(x, y, indexing='ij')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        print(Z.shape, X.shape, Y.shape)
        ax.plot_surface(X, Y, Z,cmap=cm.jet)
        

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        plt.show()
    def test_rms_slope(self):
        rms_slope_sq = 0
        count = 0
        Z = self.heights.real
        Nx, Ny = Z.shape
        x = np.arange(0, Nx, 1)
        y = np.arange(0, Ny, 1)
        X, Y = np.meshgrid(x, y, indexing='ij')
        for i in range(Nx):
            for j in range(Ny):
                if i+1 < Nx:
                    dx = X[i+1, j] - X[i, j]
                    dz = Z[i+1, j] - Z[i, j]
                    rms_slope_sq += (dz/dx)**2
                    count += 1
                if j+1 < Ny:
                    dy = Y[i, j+1] - Y[i, j]
                    dz = Z[i, j+1] - Z[i, j]
                    rms_slope_sq += (dz/dy)**2
                    count += 1
                if i+1 < Nx and j+1 < Ny:
                    d = (Y[i+1, j+1] - Y[i, j])**2 + (X[i+1, j+1] - X[i, j])**2
                    dz = Z[i+1, j+1] - Z[i, j]
                    rms_slope_sq += (dz)**2/d
                    count += 1
        rms_slope_sq /= count
        rms_roughness = np.mean(np.square(Z))
        print("rms slope squared: %g rms roughness: %g" %(rms_slope_sq,rms_roughness))
        return True
    def vis_psd(self):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        cs, qs = [], []
        Dr, Dc = self.psd.shape
        for r in range(Dr):
            for c in range(Dc):
                if self.psd[r, c] == 0: continue
                qs.append( (self.q_sq_mat[r,c])**(1/2) )
                cs.append(self.psd[r, c])
        ax1.imshow(np.fft.fftshift(self.psd), extent=(-np.pi, np.pi, -np.pi, np.pi))
        #ax1.colorbar()
        ax2.scatter(qs, cs)
        ax2.set_xlabel("log q")
        ax2.set_ylabel("log C")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        plt.show()
    def write_lammps_file(self, filename):
        d = data()
        d.title = "Lammps data file; Rough surface stats - " + self.spectrum.get_info()
        d.headers = {}
        d.sections = {}
        dd = 2**(1/6)
        dx, dy = dd/2, dd/2
        atom_type       = 1
        M, N      = self.heights.shape
        Zs       = self.heights.real
        xlo, ylo = - (M*dx)/2, -(N*dy)/2
        d.headers["zlo zhi"] = (np.min(Zs), np.max(Zs))
        d.headers["xlo xhi"] = (xlo, xlo + M*dx)
        d.headers["ylo yhi"] = (ylo, ylo + N*dy)
        d.headers["atoms"]   = M*N
        d.headers["atom types"] = atom_type
        atom_lines      = []
        velocity_lines  = []
        atom_id         = 1
        for i in range(M):
            for j in range(N):
                atom_line = "%d %d %d %g %g %g %d %d %d\n" % (atom_id, 0, atom_type, xlo+i*dx, ylo+j*dy, Zs[i, j], 0, 0, 0)
                atom_lines.append(atom_line)
                velocity_lines.append("%d 0 0 0\n" %atom_id)
                atom_id   += 1        
        d.sections["Atoms"]       = atom_lines
        d.sections["Velocities"]  = velocity_lines
        d.sections["Masses"]      = ["%d %g\n" %(atom_type, 1)]
        d.write(filename)


def main():
    parser = argparse.ArgumentParser(description = "Rough surface generation")
    parser.add_argument('--qs', type=float, default = 11.2,    help = 'Longest wave vector  (shortest wavelength).')
    parser.add_argument('--qr', type=float, default = 4,    help = 'Roll-off wave vector (Intermediate wavelength).')
    parser.add_argument('--qL', type=float, default = 3,    help = 'Shortest wave vector (Longest wavelength).')
    parser.add_argument('--H',  type=float, default = 0.5,    help = 'Hurst exponent')
    parser.add_argument('--C',  type=float, default = 1,      help = 'Power spectrum constant')
    parser.add_argument('--Lx', type=float, default = 96.7,     help = 'Length in x axis.')
    parser.add_argument('--Ly', type=float, default = 96.7,     help = 'Length in y axis.')
    parser.add_argument('--Nx', type=int  , default = 100,     help = 'Length in x axis.')
    parser.add_argument('--Ny', type=int  , default = 100,     help = 'Length in y axis.')
    args = parser.parse_args()
    args.C = 1
    args.qL = 1#2 * np.pi / args.Lx
    args.qr = 2#2 * args.qL
    args.qs = 16#2 * np.pi
    args.H = 0.8
    spectrum   = Spectrum(args.qL, args.qr, args.qs, args.H, args.C)
    d = 2**(1/6)
    dx, dy = d/2, d/2
    args.Lx, args.Ly = 150, 150
    
    args.Nx, args.Ny = int(args.Lx/dx), int(args.Ly/dy)
    N_vec = (args.Nx, args.Ny)
    L_vec = (args.Lx, args.Ly)
    r = RoughSurface(L_vec, N_vec, spectrum)
    r.filtering_algorithm_surface()
    #r.random_phase_surface()
    r.vis_surface()
    r.vis_psd()
    #r.test_rms_slope()
    r.write_lammps_file("../lammpsinput/rough_Lx%d_Ly%d.data" %(args.Lx, args.Ly))

if __name__=="__main__":
    main()
'''
a = 2**(1/6)/2
L = 100
qs = 2 * np.pi / a
qL = 2 * np.pi / L
H = 0.5
print(qs, qL, 2*qL)
spectrum   = Spectrum(qL, qL*2, qs, H, 1)   
L_vec = (L,L) 
N_vec = (2*L,2*L)  
r = RoughSurface(L_vec, N_vec, spectrum)
r.filtering_algorithm_surface()
#r.random_phase_surface()
r.vis_surface()
#r.vis_psd()'''