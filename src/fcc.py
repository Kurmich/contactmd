import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from atom import Atom, GFMD, AtomicForces
from data import data
import argparse

#class Substrate:


class Crystal:

    def __init__(self, d, exposed_surface, type, mass, atom_style):
        self.d = d #equilibrium distance between two nearest neighbors
        self.type = type
        self.mass = mass
        self.atom_style = atom_style
        self.exposed_surface = exposed_surface
        self.atom_planes = []
        self.init_borders()
        self.num_atoms = 0
        self.name = "plane"

    def print_lattice_params(self):
        print("Nearest neighbor distance: %f, atom type: %d atom mass: %d exposed surface: %s" %(self.d, self.type, self.mass, self.exposed_surface))
        print("Generator vector dx: %f dy: %f dz: %f" %(self.dx, self.dy, self.dz))

    def print_lattice_borders(self):
        print("Borders min x: %f max x: %f" %(self.x_min, self.x_max))
        print("Borders min y: %f max y: %f" %(self.y_min, self.y_max))
        print("Borders min z: %f max z: %f" %(self.z_min, self.z_max))

    def init_borders(self):
        self.x_max, self.y_max, self.z_max = float('-inf'), float('-inf'), float('-inf')
        self.x_min, self.y_min, self.z_min = float('inf'), float('inf'), float('inf')

    def update_borders(self, x, y, z):
        #update max borders                                                                                                                                                                              
        self.x_max = max(self.x_max, x)
        self.y_max = max(self.y_max, y)
        self.z_max = max(self.z_max, z)
        #update min borders                                                                                                                                                                              
        self.x_min = min(self.x_min, x)
        self.y_min = min(self.y_min, y)
        self.z_min = min(self.z_min, z)

    def shift_borders(self, dx, dy, dz):
        #shift max values                                                                                                                                                                                
        self.x_max -= dx
        self.y_max -= dy
        self.z_max -= dz
        #shift min values                                                                                                                                                                                
        self.x_min -= dx
        self.y_min -= dy
        self.z_min -= dz

    def shift_cm_to_origin(self):
        """shift coordintes of atoms to have center of mass of at the origin"""
        count = 0
        xcm, ycm, zcm = 0, 0, 0
        for atom_plane in self.atom_planes:
            for atom in atom_plane:
                count += 1
                xcm += atom.x
                ycm += atom.y
                zcm += atom.z
        xcm /= count
        ycm /= count
        zcm /= count
        self.shift_all_atom_coords(xcm, ycm, zcm)

    def top_plane_to_z0(self):
        diff = self.z_max
        self.shift_all_atom_coords(0, 0, diff)

    def shift_all_atom_coords(self, dx, dy, dz):
        """Transform atom coordiates from (x,y,z) -> (x-dx, y-dy, z-dz)"""
        self.shift_borders(dx, dy, dz)
        for atom_plane in self.atom_planes:
            for atom in atom_plane:
                atom.shift_coordinates(dx, dy, dz)

    def set_num_of_planes(self, num_planes):
        """set number of planes of the crystal to a specific number"""
        while len(self.atom_planes) > num_planes:
            self.atom_planes.pop()
        self.num_atoms = 0
        self.init_borders()
        for atom_plane in self.atom_planes:
            for atom in atom_plane:
                self.update_borders(atom.x, atom.y, atom.z)
                self.num_atoms += 1

    def rot_z(self, alpha):
        """Roatate about z axis by angle alpha (must be supplied in degrees)"""
        alpha = math.radians(alpha)
        for atom_plane in self.atom_planes:
            for atom in atom_plane:
                new_x = math.cos(alpha) * atom.x - math.sin(alpha) * atom.y
                new_y = math.sin(alpha) * atom.x + math.cos(alpha) * atom.y
                new_z = atom.z
                atom.update_coordinates(new_x, new_y, new_z)

    def bend_to_sphere(self, radius, concave_up = True):
        #bend the atoms around center
        self.name = "sphere"
        n = 0
        z_cut = 30 * self.d
        zcm = self.z_max + radius
        self.init_borders()
        self.num_atoms = 0
        #self.z_max, self.z_min = float('-inf'), float('inf')
        new_atom_planes = []
        for atom_plane in self.atom_planes:
            radius -= n*self.dz
            new_atom_plane = []
            for atom in atom_plane:
                h = math.sqrt(atom.x**2 + atom.y**2)
                if h > radius:
                    continue
                theta = math.asin(h/radius)
                z_new = zcm - radius * math.cos(theta)
                if z_new > z_cut:
                    continue
                atom.z = z_new
                new_atom_plane.append(atom)
                self.update_borders(atom.x, atom.y, atom.z)
                self.num_atoms += 1
                """
                z_sqrd = (radius-n*self.dz)**2 - atom.x**2 - atom.y**2
                z_new = math.sqrt(z_sqrd)
                if concave_up == True:
                    z_new = -z_new
                atom.z = z_new
                self.z_min = min(self.z_min, z_new)
                self.z_max = max(self.z_max, z_new)"""
            new_atom_planes.append(new_atom_plane)
            n += 1
        self.atom_planes = new_atom_planes
        return self

    def bend_to_conesphere(self, radius, angle):
        self.name = "cone"
        n = 0
        h = 25
        self.num_atoms = 0
        old_z_max = self.z_max
        zcm = self.z_max + radius
        angle_rad = math.radians(angle)
        new_atom_planes = []
        for atom_plane in self.atom_planes:
            radius -= n*self.dz
            r = radius
            dz = r - r * math.cos(angle_rad)
            rlo = r * math.sin(angle_rad)
            rhi = rlo + h /  math.tan( angle_rad )
            maxh = dz + h
            new_atom_plane = []
            for atom in atom_plane:
                d = math.sqrt(atom.x**2 + atom.y**2)
                if d < rlo:
                    theta = math.asin(d/radius)
                    z_new = zcm - radius * math.cos(theta)
                    atom.z = z_new
                    new_atom_plane.append(atom)
                else:
                    if d > self.x_max:
                        continue
                    z_new = old_z_max + dz + (d - rlo) * math.tan(angle_rad)
                    atom.z = z_new
                    new_atom_plane.append(atom)  
                self.update_borders(atom.x, atom.y, atom.z)
                self.num_atoms += 1
            n += 1
            new_atom_planes.append(new_atom_plane)
        self.atom_planes = new_atom_planes
        return self

    def visualize(self):
        """visualize lattice in a 3d plot"""
        fig = plt.figure()
        ax = Axes3D(fig)
        x, y, z = [], [], []
        for atom_plane in self.atom_planes:
            for atom in atom_plane:
                x.append(atom.x)
                y.append(atom.y)
                z.append(atom.z)
        ax.scatter(x, y, z)
        ax.set_xlabel(r'x(a)')
        ax.set_ylabel(r'y(a)')
        ax.set_zlabel(r'z(a)')
        plt.show()
    def write_lammps_file(self):
        d          = data()
        d.title    = "Crsytal with exposed surface: %s and interactomic distance: %g" %(self.exposed_surface, self.d)
        d.sections = {}
        d.headers  = {}
        d.headers["xlo xhi"] = (self.x_min , self.x_max)
        d.headers["ylo yhi"] = (self.y_min , self.y_max)
        d.headers["zlo zhi"] = (self.z_min , self.z_max)
        d.headers["atoms"]   = self.num_atoms
        unique_types = set()
        atom_lines      = []
        velocity_lines  = []
        for atom_plane in self.atom_planes:
            for atom in atom_plane:
                atom_lines.append( atom.get_contents() + "\n" )
                velocity_lines.append(atom.get_velocities() + "\n")
                unique_types.add(atom.type)
        d.headers["atom types"]   = len(unique_types)
        d.sections["Atoms"]       = atom_lines
        d.sections["Velocities"]  = velocity_lines
        Lx = int(self.x_max - self.x_min)
        filename = "%s_Lx%d.data" %(self.name, Lx)
        d.write(filename)
        





class FCC(Crystal):
    def __init__(self, d, exposed_surface, type, mass, atom_style):
        Crystal.__init__(self, d, exposed_surface, type, mass, atom_style)
        if exposed_surface == '001':
            self.set_expd_surf_to_001()
        
    def set_expd_surf_to_001(self):
        self.dx = self.d 
        self.dy = self.d 
        self.dz = self.d * math.cos( math.radians( 45 ) )
    def create_atoms(self, Nx, Ny, Nz):
        if self.exposed_surface == '001':
            self.create_atoms_001(Nx, Ny, Nz)
    def create_atoms_001(self, Nx, Ny, Nz):
        #create atoms
        cur_id = 0
        for nz in range(Nz):
            atom_plane = []
            for ny in range(Ny):
                for nx in range(Nx):
                    self.num_atoms += 1
                    x, y, z = self.get_coords_001(nx, ny, nz)
                    self.update_borders(x, y, z)
                    if self.atom_style == 'gfmd':
                        new_atom = GFMD(cur_id, self.type, x, y, z)
                        new_atom.set_grid_indices(nx, ny, nz)
                    elif self.atom_style == 'bond' or self.atom_style == 'atomic':
                        new_atom = Atom(cur_id, self.type, x, y, z)
                    atom_plane.append(new_atom)
            self.atom_planes.append(atom_plane)
        print("For Nx: %d Ny: %d Nz: %d created %d atoms in total" %(Nx, Ny, Nz, self.num_atoms))
    def get_coords_001(self, nx, ny, nz):
        """Get coordintes (x,y,z) for grid indices (nx, ny, nz)"""
        z = -nz * self.dz #crystal grows "downwards on z axis"
        shift_x, shift_y = 0, 0
        if nz%2 == 1: #check
            shift_x = self.dz * math.cos( math.radians( 45 ) )
            shift_y = shift_x
        y = ny * self.dy + shift_y
        x = nx * self.dx + shift_x
        return x, y, z



def main():
    parser = argparse.ArgumentParser(description = "Rough surface generation")
    parser.add_argument('--d',    type=float,    default = 2**(1/6),    help = 'Interatomic distance.')
    parser.add_argument('--surface', type=str,   default = '001',    help = 'Exposed surface of a crystal.')
    parser.add_argument('--atom_style', type=str, default = 'atomic',    help = 'Style: gfmd, bond, atomic')
    parser.add_argument('--atype',  type=int, default = 2**(1/6),    help = 'Atom type')
    parser.add_argument('--mass',  type=float, default = 1.0,    help = 'Atom mass')
    parser.add_argument('--r',  type=float, default = 25,      help = 'Radius of a sphere.')
    parser.add_argument('--Lx', type=float, default = 100,     help = 'Length in x axis.')
    parser.add_argument('--Ly', type=float, default = 100,     help = 'Length in y axis.')
    parser.add_argument('--Lz', type=float, default = 1,       help = 'Length in z axis.')
    parser.add_argument('--cang',    type=float, default = 45,     help = 'Angle the cone surface makes with the horizontal plane')
    parser.add_argument('--sphere', action = 'store_true', default = False, help = 'True if analysis of hardness is LJ bond lengths is needed')
    parser.add_argument('--cone', action = 'store_true', default = False, help = 'True if analysis of end-end polymer lengths is LJ bond lengths is needed')
    args = parser.parse_args()
    
    #set parameters for an fcc crystal
    args.d /= 2
    fcc     = FCC(args.d, args.surface, args.atype, args.mass, args.atom_style)
    fcc.print_lattice_params()
    
    #set dimentions
    Nx, Ny, Nz = int(math.floor(args.Lx/args.d)), int(math.floor(args.Ly/args.d)), int(math.ceil(args.Lz/args.d))
    fcc.create_atoms(Nx, Ny, Nz)
    fcc.shift_cm_to_origin()
    
    if args.sphere:
        fcc.bend_to_sphere(args.r)
    elif args.cone:
        fcc.bend_to_conesphere(args.r, args.cang)
        pass
        
    fcc.print_lattice_borders()
    fcc.visualize()
    fcc.write_lammps_file()

if __name__ == "__main__":
    main()
