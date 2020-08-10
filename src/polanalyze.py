import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from data import data
from atom import Monomer
from boxcell import *
from topology import *

forces = [500.0]
poisson = 0.5
G_shear_mod = 16.0
E_star = 4 * G_shear_mod
R = 1000.0
sigma = 1.0 #2**(1/6)
d = 2.0**(1.0/6.0) * sigma
atom_N = [1, 2, 3, 4]

hkeywords = ["atoms", "ellipsoids", "lines", "triangles", "bodies",
             "bonds", "angles", "dihedrals", "impropers",
             "atom types", "bond types", "angle types", "dihedral types",
             "improper types", "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"]
skeywords = [["Masses", "atom types"],
             ["Atoms", "atoms"], ["Ellipsoids", "ellipsoids"],
             ["Lines", "lines"], ["Triangles", "triangles"], ["Bodies", "bodies"],
             ["Bonds", "bonds"],
             ["Angles", "angles"], ["Dihedrals", "dihedrals"],
             ["Impropers", "impropers"], ["Velocities", "atoms"],
             ["Pair Coeffs", "atom types"],
             ["Bond Coeffs", "bond types"], ["Angle Coeffs", "angle types"],
             ["Dihedral Coeffs", "dihedral types"],
             ["Improper Coeffs", "improper types"],
             ["BondBond Coeffs", "angle types"],
             ["BondAngle Coeffs", "angle types"],
             ["MiddleBondTorsion Coeffs", "dihedral types"],
             ["EndBondTorsion Coeffs", "dihedral types"],
             ["AngleTorsion Coeffs", "dihedral types"],
             ["AngleAngleTorsion Coeffs", "dihedral types"],
             ["BondBond13 Coeffs", "dihedral types"],
             ["AngleAngle Coeffs", "improper types"],
             ["Molecules", "atoms"]]


class PolymerMelt:
    def __init__(self, polymers, headers, sections):
        self.polymers  = polymers
        self.headers   = headers
        self.sections  = sections
        self.Lx = headers["xlo xhi"][1] - headers["xlo xhi"][0]
        self.Ly = headers["ylo yhi"][1] - headers["ylo yhi"][0]
        self.Lz = headers["zlo zhi"][1] - headers["zlo zhi"][0]
        self.chain_count   = len(polymers)
        self.monomer_count = len(polymers[0].monomers)
        self.sorted_ids    = False

    def plot_polymer(self, pol_id):
        polymer = self.polymers[pol_id]
        xs, ys, zs = [], [], []
        fig = plt.figure()
        ax = Axes3D(fig)
        for monomer in polymer.monomers:
            xs.append(monomer.x)
            ys.append(monomer.y)
            zs.append(monomer.z)
        for i in range(len(xs)-1):
            x_next = xs[i] + self.get_displ_pbr(xs[i+1], xs[i], self.Lx)
            y_next = ys[i] + self.get_displ_pbr(ys[i+1], ys[i], self.Ly)
            z_next = zs[i] + self.get_displ_pbr(zs[i+1], zs[i], self.Lz)
            ax.plot([xs[i], x_next], [ys[i], y_next], [zs[i], z_next])
        ax.scatter(xs, ys, zs)
        plt.show()
    def get_polymer(self, pol_id):
        return self.polymers[pol_id]

    def sort_ids(self):
        if self.sorted_ids:
            return
        mon_start_id = 0 
        for pi in range(len(self.polymers)):
            polymer = self.polymers[pi]
            polymer.pol_id = pi + 1
            for mi in range(len(polymer.monomers)):
                monomer = polymer.monomers[mi]
                monomer.id = mon_start_id + mi + 1
            mon_start_id += len(polymer.monomers)
        self.sorted_ids = True

    def write_sections(self):
        if not self.sorted_ids:
            self.sort_ids()
        atom_lines      = []
        vel_lines       = []
        bond_lines      = []
        angle_lines     = []
        bond_id         = 1
        angle_id        = 1
        bond_type       = 1 #THIS IS ADHOC ONLY FOR SIGNLE BOND TYPE
        angle_type      = 1 #THIS IS ADHOC ONLY FOR SIGNLE ANGLE TYPE
        for polymer in self.polymers:
            for mi in range(len(polymer.monomers)):
                mon       = polymer.monomers[mi]
                pol_id    = polymer.pol_id if mon.mol_id == -1 else mon.mol_id
                atom_line = "%d %d %d %g %g %g %d %d %d\n" %(mon.id, pol_id, mon.type, mon.x, mon.y, mon.z, 0, 0, 0)
                vel_line  = "%d %g %g %g\n" %(mon.id, mon.vx, mon.vy, mon.vz)
                atom_lines.append(atom_line)
                vel_lines.append(vel_line)
                if mi < len(polymer.monomers)-1:
                    next_mon      = polymer.monomers[mi+1]
                    bond_line     = "%d %d %d %d\n" %(bond_id, bond_type, mon.id, next_mon.id)
                    bond_id       += 1
                    bond_lines.append(bond_line)
                if mi < len(polymer.monomers)-2:
                    next_mon      = polymer.monomers[mi+1]
                    next_next_mon = polymer.monomers[mi+2]
                    angle_line    = "%d %d %d %d %d\n" %(angle_id, angle_type, mon.id, next_mon.id, next_next_mon.id)
                    angle_id      += 1
                    angle_lines.append(angle_line)
                    
        #add angle headers
        self.headers["angles"]       = angle_id - 1 #number of angles
        self.headers["angle types"]  = angle_type
        self.headers["bonds"]        = bond_id - 1  #number of bonds
        
        
        self.sections["Atoms"]       = atom_lines
        self.sections["Velocities"]  = vel_lines
        self.sections["Bonds"]       = bond_lines
        self.sections["Angles"]      = angle_lines
    def write_lammps_file(self, filename):
        d = data()
        d.title = "Lammps data file; Equilibrated Polymer Melt"
        d.headers = self.headers
        d.sections = self.sections
        d.write(filename)
    
    def mountain_mol_ids(self):
        '''This is required for Monte-Carlo swaps during equilibration'''
        for polymer in self.polymers:
            N = len(polymer.monomers)
            lo, hi = 0, N-1
            i = 1
            while lo <= hi:
                polymer.monomers[lo].mol_id = i
                polymer.monomers[hi].mol_id = i
                lo += 1
                hi -= 1
                i  += 1

    def plot_mean_square(self, c, lbl):
        b = 0.97        
        r2 = [0 for i in range(self.monomer_count-1)]
        xvals = [i for i in range(1, self.monomer_count)]
        for pol_id in range(self.chain_count):
            yvals = self.get_mean_square2(pol_id)
            for i in range(len(yvals)):
                r2[i] += yvals[i]
        for i in range(len(r2)):
            r2[i] /= (self.chain_count)
        for y in r2:
            print(y)
        plt.semilogx(xvals, r2, c, label = lbl)
    def get_mean_square2(self, pol_id):
        pmr = self.polymers[pol_id] #select the polymer
        r2s = [] #array of mean square internal distances for polymer pol_id
        #initialize pointwise displacements
        dxs = [0] 
        dys = [0] 
        dzs = [0] 
        for i in range(self.monomer_count-1):
            mon_next = pmr.monomers[i+1]
            mon_cur  = pmr.monomers[i]
            #add pointswise distpacements between monomer i and monomer i+1
            dxs.append( self.get_displ_pbr(mon_next.x, mon_cur.x, self.Lx) )
            dys.append( self.get_displ_pbr(mon_next.y, mon_cur.y, self.Ly) )
            dzs.append( self.get_displ_pbr(mon_next.z, mon_cur.z, self.Lz) )
        #update arrays of pointwise displacements to cumulative displacements (in place update)
        self.cumsum(dxs)
        self.cumsum(dys)
        self.cumsum(dzs)
        for n in range(1, self.monomer_count):
            cur_r2 = 0
            r2 = 0
            for i in range(self.monomer_count-n):
                j = i + n 
                cur_r2 = (dxs[j] - dxs[i])**2 + (dys[j] - dys[i])**2 + (dzs[j]-dzs[i])**2
                r2 += cur_r2

            r2 /= ((self.monomer_count - n) * n)
            r2s.append(r2)

        return r2s
    def cumsum(self, arr):
        #update array to cumulative sum
        cur = arr[0]
        for i in range(1, len(arr)):
            cur += arr[i]
            arr[i] = cur
        
    def get_displ_pbr(self, x_next, x_prev, L):
        #get displacement vector pointing form x_prev to x_next for periodic boundary condition
        #periodicity is L
        sign = np.sign(x_next - x_prev)
        dx = abs(x_next - x_prev)
        if dx > L/2:
            dx = -(L - dx)
        return sign * dx
    
    def get_Florys_ratio(self):
        count      = 0
        sum_cosine = 0
        for polymer in self.polymers:
            for mi in range(len(polymer.monomers)-2):
                mon_cur       = polymer.monomers[mi]
                mon_next      = polymer.monomers[mi+1]
                mon_next_next = polymer.monomers[mi+2]
                dx1 = self.get_displ_pbr(mon_next.x, mon_cur.x, self.Lx)
                dy1 = self.get_displ_pbr(mon_next.y, mon_cur.y, self.Ly)
                dz1 = self.get_displ_pbr(mon_next.z, mon_cur.z, self.Lz)
                dx2 = self.get_displ_pbr(mon_next_next.x, mon_next.x, self.Lx)
                dy2 = self.get_displ_pbr(mon_next_next.y, mon_next.y, self.Ly)
                dz2 = self.get_displ_pbr(mon_next_next.z, mon_next.z, self.Lz)
                num = dx1 * dx2 + dy1 * dy2 + dz1 * dz2
                den = ((dx1**2 + dy1**2 + dz1**2) * (dx2**2 + dy2**2 + dz2**2))**(1/2)
                sum_cosine += (num/den)
                count   += 1
        exp_cosine = sum_cosine/count #expectation value of bond angle
        c_inf = (1 + exp_cosine) / (1 - exp_cosine)
        print("Florys characteristic ratio: %g" %c_inf)
        return c_inf
    def bond_isotropy(self):
        count  = 0
        cos_sq_x  = 0
        cos_sq_y  = 0
        cos_sq_z  = 0
        for polymer in self.polymers:
            for mi in range(len(polymer.monomers)-1):
                mon_cur       = polymer.monomers[mi]
                mon_next      = polymer.monomers[mi+1]
                dx            = self.get_displ_pbr(mon_next.x, mon_cur.x, self.Lx)
                dy            = self.get_displ_pbr(mon_next.y, mon_cur.y, self.Ly)
                dz            = self.get_displ_pbr(mon_next.z, mon_cur.z, self.Lz)
                bond_len_sq   = (dx**2 + dy**2 + dz**2)
                cos_sq_x      += (dx**2) / bond_len_sq
                cos_sq_y      += (dy**2) / bond_len_sq
                cos_sq_z      += (dz**2) / bond_len_sq
            count   += (len(polymer.monomers)-1)
        cos_sq_x  /= count
        cos_sq_y  /= count
        cos_sq_z  /= count
        fx   = (3*cos_sq_x - 1)/2
        fy   = (3*cos_sq_y - 1)/2
        fz   = (3*cos_sq_z - 1)/2
        print("f_x: %g f_y: %g f_z: %g" %(fx, fy, fz))
        
    def partition_chains(self, chain_len):
        assert chain_len < self.monomer_count and self.monomer_count % chain_len == 0
        branch_count = self.monomer_count / chain_len
        new_polymers = []
        new_pol_id = 1
        for polymer in self.polymers:
            for mi in range(0, len(polymer.monomers), chain_len):
                oligomer             = Polymer(new_pol_id)
                oligomer.monomers    = polymer.monomers[mi:mi+chain_len]
                new_pol_id           += 1
                new_polymers.append(oligomer)
        new_monomer_count = len(new_polymers[0].monomers)
        assert len(new_polymers) == self.chain_count * branch_count
        assert new_monomer_count == chain_len
        
        return PolymerMelt(new_polymers, self.headers, self.sections)
                

class Polymer:
    def __init__(self, pol_id):
        self.pol_id = pol_id
        self.monomers = []
    def add_monomer(self, monomer):
        self.monomers.append(monomer)
        
class Node:
    def __init__(self):
        self.pol_id = -1
        self.monomer = Monomer(0,0,0,0,0,0,0,0)
        self.neighbors = []
    def __lt__(self, other):
        return self.monomer.id < other.monomer.id
    def __eq__(self, other):
        return self.monomer.id == other.monomer.id


class Graph:
    def __init__(self, node_count):
        self.nodes = []
        for i in range(node_count):
            self.nodes.append( Node() )
    
    def update_mon_pos(self, mon_id, type, x, y, z):
        monomer = self.nodes[mon_id].monomer
        monomer.id = mon_id
        monomer.type = type
        monomer.x, monomer.y, monomer.z = x, y, z

    def update_mon_vel(self, mon_id, vx, vy, vz):
        monomer = self.nodes[mon_id].monomer
        monomer.set_velocity(vx, vy, vz)

    def add_bond(self, mon_id1, mon_id2):
        self.nodes[mon_id1].neighbors.append(self.nodes[mon_id2])
        self.nodes[mon_id2].neighbors.append(self.nodes[mon_id1])
        #print("adding bonds %d %d" %(mol_id1, mol_id2))

    def cluster_nodes(self):
        cur_id = 0
        self.clusters = []
        for node in self.nodes:
            if node.pol_id < 0:
                group = []
                self.dfs(node, cur_id, group)
                cur_id += 1
                self.clusters.append(group)
        '''for p in self.polymers:
            print(len(p.monomers))'''
        print("Number of clusters: %d" %cur_id)

    def dfs(self, root, id, group):
        if root.pol_id != -1:
            return
        root.pol_id = id
        group.append(root)
        for nbr in root.neighbors:
            self.dfs(nbr, id, group)

    def group_polymers(self):
        polymers = []
        if self.clusters is None:
            self.cluster_nodes()
        pol_id = 0
        for group in self.clusters:
            self.arrange_ids(group)
            p = Polymer(pol_id)
            #print("next")
            for node in sorted(group):
                #print(len(node.neighbors))
                p.add_monomer(node.monomer)
            polymers.append(p)
            pol_id += 1
        print(pol_id)
        return polymers
    def arrange_ids(self, group):
        #arrange polymer id's
        start_node = None
        for node in group:
            if len(node.neighbors) == 1:
                start_node = node
                break
        self.create_ids(start_node, start_node, 0)
    def create_ids(self, parent, node, id):
        node.monomer.id = id
        for nbr in node.neighbors:
            if nbr != parent:
                self.create_ids(node, nbr, id+1)


def get_graph(filename, chain_count, monomer_count):
    total_monomers = monomer_count * chain_count
    graph = Graph(total_monomers)
    with open(filename) as file:
        title = file.readline()
        print("Title of the data file: %s" %title)
        headers = {}
        while 1:
            line = file.readline()
            if line.isspace():
                continue
            found = False
            for keyword in hkeywords:
                if line.find(keyword) >= 0:
                    found = True
                    words = line.split()
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or keyword == "zlo zhi":
                        headers[keyword] = (float(words[0]), float(words[1]))
                    elif keyword == "xy xz yz":
                        headers[keyword] = (float(words[0]), float(words[1]), float(words[2]))
                    else:
                        headers[keyword] = int(words[0])
            if not found:
                break

        sections = {}
        while 1:
            found = False
            for pair in skeywords:
                keyword, length = pair[0], pair[1]
                if keyword in line:
                    found = True
                    if length not in headers.keys():
                        raise Exception("data section %s has no matching header value" % line)
                    file.readline()
                    list = []
                    print(keyword)
                    for i in range(headers[length]):
                        list.append(file.readline())
                    sections[keyword] = list
            if not found:
                raise Exception("invalid section %s in data file" % line)
            file.readline()
            line = file.readline()
            if not line:
                break
            line = line.strip()

        for keyword in sections.keys():
            for line in sections[keyword]:
                if keyword == "Atoms":
                    mon_id, _, type, x, y, z, _, _, _ = line.split(' ')
                    mon_id, type, x, y, z = int(mon_id), int(type), float(x), float(y), float(z)
                    mon_id -= 1
                    graph.update_mon_pos(mon_id, type, x, y, z)
                elif keyword == "Bonds":
                    bond_id, bond_type, mon_id1, mon_id2 = line.split(' ')
                    bond_id, bond_type, mon_id1, mon_id2 = int(bond_id), int(bond_type), int(mon_id1)-1, int(mon_id2)-1
                    graph.add_bond(mon_id1,  mon_id2)
                elif keyword == "Velocities":
                    mon_id, vx, vy, vz = line.split(' ')
                    mon_id, vx, vy, vz = int(mon_id),  float(vx), float(vy), float(vz)
                    mon_id -= 1
                    graph.update_mon_vel(mon_id,  vx, vy, vz)

    graph.cluster_nodes()
    return graph, headers, sections

def read_goal(filename):
    xs, ys = [], []
    with open(filename, 'r') as f:
        for line in f:
#            print(line)
            x, y = line.split(' ')
            x, y = float(x), float(y)
            xs.append(x)
            ys.append(y)
    return xs, ys

#def get_edge_avg(atom_forces):


def vis_layers(M, N, T):
    filename = "../lammpsinput/data_quenched_stiff_M%d_N%d_T%g_nve_smooth" %(M,N,T)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    all_monomers = []
    for polymer in pol_melt.polymers:
        all_monomers.extend(polymer.monomers)
    plot_layer_density(all_monomers)
    

def add_angles(M, N):
    #'''  
    filename = "../lammpsinput/melt_stiff_wallz_M%d_N%d.data" %(M,N)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    pol_melt.mountain_mol_ids()
    pol_melt.get_Florys_ratio()
    pol_melt.write_sections()
    pol_melt.write_lammps_file("../lammpsinput/melt_wallz_stiff_M%d_N%d.data" %(M,N))
    pol_melt.plot_mean_square('b', 'Initial')
    
def clean_quenched_file(M, N, T):
    xs, ys = read_goal("goal.txt")
    filename = "../lammpsinput/data_quenched_stiff_M%d_N%d_T%g_nve_smooth" %(M,N,T)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    pol_melt.write_sections()
    pol_melt.get_Florys_ratio()
    pol_melt.write_lammps_file("../lammpsinput/clean_quenched_stiff_M%d_N%d_T%g_smooth.data" %(M,N,T))
#    d = data("eq_M%d_N%d.data" %(M, N))

    pol_melt.plot_mean_square('r', 'Equilibrated')
    plt.suptitle('Mean Square Internal Distances M: %d N: %d T: %g' %(M, N, T), fontsize = 20)
    plt.plot(xs,ys, 'g', label='Target function')
    plt.ylabel(r'$<R^2(n)>/n$', fontsize = 16)
    plt.xlabel(r'$n$', fontsize = 16)
    plt.legend()
    plt.show()
    pol_melt.plot_polymer(3)


def check_equilibration(M, N):
    xs, ys = read_goal("goal.txt")
    filename = "../lammpsinput/data_eq_stiff_M%d_N%d" %(M,N)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    pol_melt.mountain_mol_ids()
    pol_melt.get_Florys_ratio()
    pol_melt.bond_isotropy()
    pol_melt.plot_mean_square('r', 'Equilibrated')
    plt.suptitle('Mean Square Internal Distances M: %d N: %d T: %g' %(M, N, 1), fontsize = 20)
    plt.plot(xs,ys, 'g', label='Target function')
    plt.ylabel(r'$<R^2(n)>/n$', fontsize = 16)
    plt.xlabel(r'$n$', fontsize = 16)
    plt.legend()
    plt.show()
    
def check_isotropy(M, N, T):
    xs, ys = read_goal("goal.txt")
    filename = "../lammpsinput/data_quenched_stiff_M%d_N%d_T%g_nve_smooth" %(M,N,T)
    #filename = "../lammpsinput/data_eq_stiff_M%d_N%d" %(M,N)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    #pol_melt.mountain_mol_ids()
    pol_melt.bond_isotropy()
    
def divide_chains(M, N, T, chain_len):
    filename = "../lammpsinput/data_quenched_stiff_M%d_N%d_T%g_nve_smooth" %(M,N,T)
    #filename = "../lammpsinput/data_eq_stiff_M%d_N%d" %(M,N)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    new_pol_melt = pol_melt.partition_chains(chain_len)
    chain_count = int(M*(N/chain_len))
    new_pol_melt.write_sections()
    new_pol_melt.write_lammps_file("../lammpsinput/clean_quenched_stiff_M%d_N%d_T%g_smooth.data" %(chain_count,chain_len,T))
    
    

def main():
    xs, ys = read_goal("goal.txt")
    M = 2000
    N = 256
    T = 0.0001
    #check_isotropy(M, N, T)
    divide_chains(M, N, T, 8)
    #check_equilibration(M, N)
    #vis_layers(M, N, T)
    #add_angles(M, N)
    #check_equilibration(M, N)
    #clean_quenched_file(M, N, T)
#    return
    '''
   
  
    filename = "../lammpsinput/melt_stiff_wallz_M%d_N%d.data" %(M,N)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    pol_melt.write_sections()
    pol_melt.write_lammps_file("../lammpsinput/melt_wallz_stiff_M%d_N%d.data" %(M,N))
    pol_melt.plot_mean_square('b', 'Initial')'''
    
    '''
    filename = "../lammpsinput/data_quenched_M%d_N%d_T%g" %(M,N,T)
    graph, headers, sections = get_graph(filename, M, N)
    polymers = graph.group_polymers()
    pol_melt = PolymerMelt(polymers, headers, sections)
    pol_melt.write_sections()
    pol_melt.write_lammps_file("../lammpsinput/clean_quenched_M%d_N%d_T%g_smooth.data" %(M,N,T))
#    d = data("eq_M%d_N%d.data" %(M, N))

    pol_melt.plot_mean_square('r', 'Equilibrated')
    plt.suptitle('Mean Square Internal Distances M: %d N: %d T: %g' %(M, N, T), fontsize = 20)
    plt.plot(xs,ys, 'g', label='Target function')
    plt.ylabel(r'$<R^2(n)>/n$', fontsize = 16)
    plt.xlabel(r'$n$', fontsize = 16)
    plt.legend()
    plt.show()
    pol_melt.plot_polymer(3)
    '''
if __name__=="__main__":
    main()

#plt.plot(timestep, forcesx)
#plt.suptitle('Friction force vs timestep', fontsize = 20)
#plt.ylabel('Fx', fontsize = 20)
#plt.xlabel('t', fontsize = 20)
#plt.axhline(0, color = 'black')
#plt.show()
