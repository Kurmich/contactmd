
from fcc import *
import random
from scipy.optimize import minimize

def create_substrate(crystal_type, exposed_surface, nearest_nbr_d,  atom_type, atom_mass, atom_style, Dx, Dy, Dz):
    d = nearest_nbr_d
    if crystal_type == 'fcc':
       # a = (2**(1/2)) * d #cube side length
        fcc = FCC(d, exposed_surface , atom_type, atom_mass, atom_style)
        Nx, Ny, Nz = int(math.ceil(Dx/d)), int(math.ceil(Dy/d)), int(math.ceil(Dz/d)) #number of atoms to create in each dimension
        if Nx%2 == 0:
            Nx += 1
            Ny += 1
        fcc.create_atoms(Nx, Ny, Nz)
        fcc.print_lattice_params()
        fcc.shift_cm_to_origin()
    substrate = fcc
    return substrate

def create_spherical_tip(substrate, radius, layers = 1):
    substrate.set_num_of_planes(layers)
    tip = substrate.bend_to_sphere(radius)
    return tip

def bring_to_contact(cr_below, cr_above, d):
    z_diff = cr_above.z_min - (cr_below.z_max + d)
    cr_above.shift_all_atom_coords(0, 0, z_diff)


def create_oligomers( substrate, ratio, n_layers, sep_z, type, n_bonds ):
    top_atom_plane = substrate.atom_planes[0]
    chained_atoms = []
    #create layers for atoms connected with spring
    for i in range(n_layers):
        chained_atoms.append([])
    prev_y = top_atom_plane[0].y #here we assume that for loop generates atoms with different x coordinates first
    idx_x = 1
    idx_y = 1
    id = 0
    substr_bond = Bond(1, [1,2])
    chain_bond = Bond(2, [3,4])
    mol_id = 2
    for atom in top_atom_plane:
        #toggle lines in a plane
        if prev_y != atom.y:
            prev_y = atom.y
            idx_x  = 1
            idx_y += 1
        #for each ratio atom create a chained layer above it
        if (idx_x % ratio == 0 and idx_y % ratio == 0) or (idx_y % ratio == 1 and (idx_x+1) % ratio == 0):
            atom.set_molecule_id(mol_id)
            for i in range(n_layers):
                id += 1
                new_z = atom.z + (i+1) * sep_z
                new_atom = Atom(id, type, atom.x, atom.y, new_z)
                new_atom.set_molecule_id(mol_id)
                x_idx, y_idx, z_idx = int(idx_x/ratio - 1), int(idx_y/ratio - 1), i
                new_atom.set_grid_indices(x_idx, y_idx, z_idx)
                #add a bond between the substrate atom and a first atom in the oligomer chain 
                if i == 0:
                    substr_bond.add_bonded_atoms(new_atom, atom)
                chained_atoms[i].append(new_atom)
                substrate.update_borders(new_atom.x, new_atom.y, new_atom.z)
            mol_id += 1
        idx_x += 1
    for i in range(len(chained_atoms[0])):
        for j in range(1, n_layers):
            chain_bond.add_bonded_atoms(chained_atoms[j][i], chained_atoms[j-1][i])
    return chained_atoms, [substr_bond, chain_bond]


def fene_pot(r, K, R_0, E, sigma):
    return -0.5 * K * (R_0**2) * math.log(1 - (r/R_0)**2) + 4 * E * ( (sigma/r)**12 - (sigma/r)**6 ) + E

class Bond:
    def __init__(self, type, coeffs):
        self.type = type
        self.coeffs = coeffs
        self.bonded_atoms = set()

    def add_bonded_atoms(self, atom1, atom2):
        if atom1.id < atom2.id:
            self.bonded_atoms.add((atom1, atom2))
        else:
            self.bonded_atoms.add((atom2, atom1))


class StructureUnitParams:
    def __init__(self, id, structure_type, crystal_type, exposed_surface, nearest_nbr_d, atom_type, atom_mass, atom_style, Dx, Dy, Dz):
        self.id = id
        self.structure_type = structure_type
        self.crystal_type = crystal_type
        self.exposed_surface = exposed_surface
        self.nearest_nbr_d = nearest_nbr_d
        self.atom_type = atom_type
        self.atom_mass = atom_mass
        self.atom_style = atom_style
        self.Dx, self.Dy, self.Dz = Dx, Dy, Dz
        self.layers = -1
    def set_num_of_layers(self, layers):
        self.layers = layers
    def set_radius(self, radius):
        self.radius = radius
    def set_cone_angle(self, ang_deg):
        self.cone_angle = ang_deg
    def __eq__(self, other):
        return self.atom_type == other.atom_type
    def __lt__(self, other):
        return self.atom_type < other.atom_type

class SimulationStructure:
    def __init__(self, units, params):
        self.units = units
        self.params = sorted(params)
        self.sim_structures = dict()
        self.bonds = []
        self.oligomers = []
        self.intrs = SimulationStructure.Interactions()
        self.x_max, self.y_max, self.z_max = float('-inf'), float('-inf'), float('-inf')
        self.x_min, self.y_min, self.z_min = float('inf'), float('inf'), float('inf')
        self.create_all_structures()
    class Interactions:
        def __init__(self):
            self.atom_count = 0
            self.atom_types = 0
            self.bond_count = 0
            self.bond_types = 0
            self.angle_count = 0
            self.angle_types = 0
            self.dihedral_count = 0
            self.dihedral_types = 0
            self.improper_count = 0
            self.improper_types = 0
    class OligomerTypeMass:
        def __init__(self, atom_type, atom_mass):
            self.atom_type = atom_type
            self.atom_mass = atom_mass

    def update_max_borders(self, x, y, z):
        #update max borders
        self.x_max = max(self.x_max, x)
        self.y_max = max(self.y_max, y)
        self.z_max = max(self.z_max, z)
    def update_min_borders(self, x, y, z):
        #update min borders
        self.x_min = min(self.x_min, x)
        self.y_min = min(self.y_min, y)
        self.z_min = min(self.z_min, z)

    def create_all_structures(self):
        unique_atoms = set()
        for p in self.params:
            if p.structure_type == 'substrate' or p.structure_type == 'lubricant':
                print("Adding substrate.")
                unit = create_substrate(p.crystal_type, p.exposed_surface, p.nearest_nbr_d, p.atom_type, p.atom_mass, p.atom_style,  p.Dx, p.Dy, p.Dz)
                if p.layers != -1:
                    unit.set_num_of_planes(p.layers)
                unit.top_plane_to_z0()
            elif p.structure_type == 'spherical tip':
                print("Adding spherical tip.")
                unit = create_substrate(p.crystal_type, p.exposed_surface, p.nearest_nbr_d, p.atom_type, p.atom_mass, p.atom_style,  p.Dx, p.Dy, p.Dz)
                unit.set_num_of_planes(p.layers)
                unit = unit.bend_to_sphere(p.radius)
            elif p.structure_type == 'cone tip':
                print("Adding spherical tip.")
                unit = create_substrate(p.crystal_type, p.exposed_surface, p.nearest_nbr_d, p.atom_type, p.atom_mass, p.atom_style,  p.Dx, p.Dy, p.Dz)
                unit.set_num_of_planes(p.layers)
                unit = unit.bend_to_conesphere(p.radius, p.cone_angle)
    
            self.intrs.atom_count += unit.num_atoms
            self.sim_structures[p.id] = unit
            unique_atoms.add(p.atom_type)
        self.intrs.atom_types = len(unique_atoms)
        #visualize_system(list(self.sim_structures.values()))
    def bring_to_contact(self, id_1, id_2):
        cr_above = self.sim_structures[id_1]
        cr_below = self.sim_structures[id_2]
        d = max(cr_above.d, cr_below.d)
        z_diff = cr_above.z_min - (cr_below.z_max + d)
        cr_above.shift_all_atom_coords(0, 0, z_diff)
    def make_substrate_glassy(self):
        substrate = self.sim_structures['substrate']
        ratio = 2
        n_layers = 4
        sep_z = (2**(1/6))
        type = len(self.sim_structures) + 1
        n_bonds = 2
        mass = 1.0
        K, R_0, E, sigma = 7.5, 1.5, 0.25, 1.0
        res = minimize(lambda x: fene_pot(x, K, R_0, E, sigma), x0=sep_z)
        # Check if the optimization was successful
        print("Minimization result")
        print(res.success)
        print(res.x[0])
        if res.success and res.x[0] > 0:
            sep_z = res.x[0]
        #substrate.print_lattice_borders()
        oligomers, bonds = create_oligomers( substrate, ratio, n_layers, sep_z, type, n_bonds )
        self.intrs.atom_types  += 1
        self.params.append(SimulationStructure.OligomerTypeMass(type, mass)) #add
        self.bonds += bonds #merge
        self.intrs.bond_types = len(self.bonds)
        self.intrs.bound_count = 0
        for bond in self.bonds:
            self.intrs.bond_count += len(bond.bonded_atoms)
            print(self.intrs.bond_count)
        self.oligomers = oligomers
        for layer in oligomers:
            self.intrs.atom_count += len(layer) 
        #substrate.print_lattice_borders()
        #self.sim_structures['substrate'].print_lattice_borders()
        return oligomers

    def perturb_oligomer_pos(self, a, b, axes):
        """Perturb oligomer atom coordinates specified by axes"""
        assert a < b, "First range variable ( %f ) should be less than the second variable ( %f )" %(a, b)
        if not isinstance(axes, str):
            raise TypeError("Axes should be given as a string e.g. 'xyz', 'xy'")
        axes = axes.lower()
        for layer in self.oligomers:
            for atom in layer:
                if 'x'  in axes:
                    atom.x = random.uniform(a, b) + atom.x
                if 'y' in axes:
                    atom.y = random.uniform(a, b) + atom.y
                if 'z' in axes:
                    atom.z = random.uniform(a, b) + atom.z
    def perturb_oligomer_angle(self, theta_max, phi_max):
        "Randomly rotate oligomers"
        #phi_max = phi_max % 180
        #theta_max = theta_max % 90
        for i in range(len(self.oligomers[0])):
            #get random angles for rotation
            phi = math.radians( random.uniform(0, phi_max) )
            theta = math.radians( random.uniform(0, theta_max))
            #get reference points x, y, z about which to perform rotations
            bottom_atom = self.oligomers[0][i]
            sep_d = self.oligomers[1][i].z - bottom_atom.z
            assert sep_d > 0, "Invalid oligomers"
            x, y, z = bottom_atom.x, bottom_atom.y, bottom_atom.z - sep_d
            #rotate atoms in a single oligomer chain
            for j in range(0, len(self.oligomers)):
                atom = self.oligomers[j][i]
                r = math.sqrt( (atom.x - x)**2 + (atom.y - y)**2 + (atom.z - z)**2 )
                atom.z = z + r * math.cos( theta )
                atom.x = x + r * math.sin( theta ) * math.cos( phi )
                atom.y = y + r * math.sin( theta ) * math.cos( phi )
               
    def rotate_unit(self, id, alpha, beta = 0, phi = 0):
        """Rotate a unit by alpha (degrees) about z axis, beta (deg) about x axis, phi (deg) about y axis, all counterclockwise"""
        unit = self.sim_structures[id]
        unit.rot_z(alpha)


    def set_borders(self):
        for key, unit in self.sim_structures.items():
            self.update_max_borders(unit.x_max, unit.y_max, unit.z_max)
            self.update_min_borders(unit.x_min, unit.y_min, unit.z_min)

    def visualize(self):
        """Visualize the simulation structure via 3d plot"""
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        i = 0
        fig = plt.figure()
        ax = Axes3D(fig)

        print("Plotting crystal structures")
        for key, unit in self.sim_structures.items():
            x, y, z = [], [], []
            for atom_plane in unit.atom_planes:
                for atom in atom_plane:
                    x.append(atom.x)
                    y.append(atom.y)
                    z.append(atom.z)
            ax.scatter(x, y, z, colors[i])
        i += 1
    
        if len(self.oligomers) > 0:
            print("Plotting oligomers")
            x, y, z = [], [], []
            for layer in self.oligomers:
                for atom in layer:
                    x.append(atom.x)
                    y.append(atom.y)
                    z.append(atom.z)
            ax.scatter(x, y, z, colors[i])
        i += 1



        if len(self.bonds) > 0:
            print("Plotting bonds")
            for bond in self.bonds:
                for pair in bond.bonded_atoms:
                    atom1, atom2 = pair[0], pair[1]
                    ax.plot([atom1.x, atom2.x], [atom1.y, atom2.y], [atom1.z, atom2.z], colors[i])
                i += 1

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    def make_ids_unique(self):
        print("Creating unique ids for atoms")
        id = 1
        for key, unit in self.sim_structures.items():
            for atom_plane in unit.atom_planes:
                    for atom in atom_plane:
                        atom.set_id(id)
                        id += 1
        for layer in self.oligomers:
            for atom in layer:
                atom.set_id(id)
                id += 1
    def create_lammps_input(self, filename):
        self.make_ids_unique()
        with open(filename, 'w') as f:
            f.write("LAMMPS data file; \n \n")
            self.print_interaction_counts(f)
            self.print_boundaries(f)
            self.print_masses(f)
            self.print_atoms(f)
            self.print_velocities(f)
            if self.intrs.bond_count > 0:
                self.print_bonds(f)
            
    def print_interaction_counts(self, file):
        file.write("%d atoms\n%d bonds\n%d angles\n%d dihedrals\n%d impropers\n\n" %(self.intrs.atom_count, self.intrs.bond_count, self.intrs.angle_count, self.intrs.dihedral_count, self.intrs.improper_count))

        if self.intrs.atom_types > 0:
            file.write("%d atom types \n" %self.intrs.atom_types)
        if self.intrs.bond_types > 0:
            file.write("%d bond types \n" %self.intrs.bond_types)
        if self.intrs.angle_types > 0:
            file.write("%d angle types \n" %self.intrs.angle_types)
        if self.intrs.dihedral_types > 0:
            file.write("%d  dihedral types \n" %self.intrs.dihedral_types)
        if self.intrs.improper_types > 0:
            file.write("%d improper_types \n" %self.intrs.improper_types)
        file.write("\n")
    def print_boundaries(self, file):
        """
        -0.5 0.5 xlo xhi       (for periodic systems this is box size,
        -0.5 0.5 ylo yhi        for non-periodic it is min/max extent of atoms)
        -0.5 0.5 zlo zhi       (do not include this line for 2-d simulations)"""
        shift = list(self.sim_structures.values())[0].d/2
        file.write("%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\n" %(self.x_min-shift, self.x_max+shift, self.y_min-shift, self.y_max+shift, self.z_min-shift, self.z_max+shift))
    def print_masses(self, file):
        """Masses

        1 mass
        ...
        N mass         (N = # of atom types)"""
        unique = set()
        file.write("Masses\n\n")
        for param in self.params:
            if param.atom_type not in unique:
                unique.add(param.atom_type)
                file.write("%d %f \n" %(param.atom_type, param.atom_mass))
        file.write("\n")
    def print_atoms(self, file):
        """Atoms

        1 molecule-tag atom-type q x y z nx ny nz  (nx,ny,nz are optional -
        ...                                    see "true flag" input command)
        ...                
        N molecule-tag atom-type q x y z nx ny nz  (N = # of atoms)"""
        file.write("Atoms\n\n")
        atom_count = 0
        for key, unit in self.sim_structures.items():
            for atom_plane in unit.atom_planes:
                for atom in atom_plane:
                    atom_count += 1
                    file.write(atom.get_contents())
                    if len(self.oligomers) > 0:
                        if atom.mol_id == -1:
                            file.write(" %d\n" %(-atom.mol_id))
                        else:
                            file.write(" %d\n" %(atom.mol_id))
                    else:
                        file.write("\n")

        for layer in self.oligomers:
            for atom in layer:
                atom_count += 1
                file.write(atom.get_contents())
                if atom.mol_id == -1:
                    file.write(" %d\n" %(-atom.mol_id))
                else:
                    file.write(" %d\n" %(atom.mol_id))
        file.write("\n")
        assert atom_count == self.intrs.atom_count, "Total number of atoms and number of printed atoms do not match"

    def print_velocities(self, file):
        file.write("Velocities\n\n")
        atom_count = 0
        vx, vy, vz = 0, 0, 0
        for key, unit in self.sim_structures.items():
            for atom_plane in unit.atom_planes:
                for atom in atom_plane:
                    atom_count += 1
                    file.write("%d %g %g %g" %(atom.id, vx, vy, vz))
                    if len(self.oligomers) > 0:
                        if atom.mol_id == -1:
                            file.write(" %d\n" %(-atom.mol_id))
                        else:
                            file.write(" %d\n" %(atom.mol_id))
                    else:
                        file.write("\n")

        for layer in self.oligomers:
            for atom in layer:
                atom_count += 1
                file.write(atom.get_contents())
                if atom.mol_id == -1:
                    file.write(" %d\n" %(-atom.mol_id))
                else:
                    file.write(" %d\n" %(atom.mol_id))
        file.write("\n")
        assert atom_count == self.intrs.atom_count, "Total number of atoms and number of printed atoms do not match"


    def print_bonds(self, file):
        """Bonds

        1 bond-type atom-1 atom-2
        ...
        N bond-type atom-1 atom-2         (N = # of bonds)"""
        file.write("Bonds\n\n")
        bond_count = 0
        for bond in self.bonds:
            type = bond.type
            for pair in bond.bonded_atoms:
                bond_count += 1
                atom1, atom2 = pair[0], pair[1]
                file.write("%d %d %d %d\n" %(bond_count, type, atom1.id, atom2.id))
        #file.write("\n")
        assert bond_count == self.intrs.bond_count, "Total number of bonds and number of printed bonds do not match"


def main():
    d = 2**(1/6)
    Dx, Dy, Dz = 50, 50, 3  #195, 195, 3
    Dz_tip = (2**(1/2)) * d
    radius = 10
    units = 'lj'
    atom_style = 'bond'
    cone_angle = 75
#    substrate = create_substrate('fcc', '001', d,  Dx, Dy, Dz)
    substrate_unit = StructureUnitParams('substrate', 'substrate', 'fcc', '001', d, 1, 1.0, atom_style,  Dx, Dy, Dz)
    d /= 2
    spherical_tip_unit = StructureUnitParams('spherical tip', 'spherical tip', 'fcc', '001', d, 1, 1.0, atom_style,  Dx, Dy, Dz_tip)
    cone_tip_unit = StructureUnitParams('cone tip', 'cone tip', 'fcc', '001', d, 1, 1.0, atom_style,  Dx, Dy, Dz_tip)
    cone_tip_unit.set_cone_angle(cone_angle)
#    lubricant = StructureUnitParams('lubricant', 'lubricant', 'fcc', '001', d, 3, 1.0,  Dx/2, Dy/2, Dz)
#    lubricant.set_num_of_layers(1)
    substrate_unit.set_num_of_layers(1)
    spherical_tip_unit.set_num_of_layers(1)
    spherical_tip_unit.set_radius(radius)
    cone_tip_unit.set_num_of_layers(1)
    cone_tip_unit.set_radius(radius)
    sim_str = SimulationStructure(units, [cone_tip_unit]) # SimulationStructure(units, [substrate_unit, spherical_tip_unit, lubricant])
#    sim_str = SimulationStructure(units, [substrate_unit, spherical_tip_unit]) # SimulationStructure(units, [substrate_unit, spherical_tip_unit, lubricant])
#    sim_str.bring_to_contact('lubricant','substrate')
 #   oligomer = sim_str.make_substrate_glassy()
#    sim_str.bring_to_contact('spherical tip','substrate')
#    sim_str.perturb_oligomer_angle(5, 180)
#    sim_str.perturb_oligomer_pos(-0.20, 0.20, 'xy')
    #sim_str.make_substrate_glassy()
#    sim_str.rotate_unit('spherical tip', 30)
    #sim_str.bring_to_contact('lubricant','substrate')
#    sim_str.make_ids_unique()
    #visualize_system(list(sim_str.sim_structures.values()), oligomer, sim_str.bonds)
    sim_str.set_borders()
    sim_str.visualize()
#    sim_str.create_lammps_input("../lammpsinput/flattip_Dx%d.dat" %(Dx))
    sim_str.create_lammps_input("../lammpsinput/tip_r%d_Dx%d_cang%d.dat" %(radius, Dx, cone_angle))
#    tip_as_substrate = create_substrate('fcc', '001', d,  Dx, Dy, Dz_tip)
#    tip = create_spherical_tip(tip_as_substrate, radius)
#    bring_to_contact(substrate, tip, d)
#    visualize_system([substrate, tip])

if __name__ == "__main__":
    main()
