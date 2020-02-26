class Atom:
    def __init__(self, id, type, x, y, z):
        """x, y, z - instantaneous positions of atoms"""
        self.id = id
        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.mol_id = 1000000
    def set_id(self, id):
        self.id = id
    def set_molecule_id(self, mol_id):
        self.mol_id = mol_id
    def shift_coordinates(self, dx, dy, dz ):
        """shift coordinates by dx, dy, dz"""
        self.x -= dx
        self.y -= dy
        self.z -= dz
    def update_coordinates(self, new_x, new_y, new_z):
        """Update coordinates of the atom from (x, y, z) -> (new_x, new_y, new_z)"""
        self.x = new_x
        self.y = new_y
        self.z = new_z
    def get_contents(self):
        return "%d %d %d %f %f %f %d %d %d" %(self.id, self.mol_id, self.type, self.x, self.y, self.z, 0, 0, 0)



class GFMD(Atom):
    def __init__(self, id, type, x, y, z):
        Atom.__init__(self, id, type, x, y, z)
        self.set_equilibrium(x, y, z)

    def set_equilibrium(self, x_eq, y_eq, z_eq):
        """set the equilibrium position of atom"""
        self.x_eq = x_eq
        self.y_eq = y_eq
        self.z_eq = z_eq

    def set_grid_indices(self, x_idx, y_idx, z_idx):
        """set grid index of the atom"""
        #check types and ranges
        assert type(x_idx) is int, "x_idx is not an integer: %r" % x_idx
        assert type(y_idx) is int, "y_idx is not an integer: %r" % y_idx
        assert type(z_idx) is int, "z_idx is not an integer: %r" % z_idx
        assert x_idx >= 0 and y_idx >= 0 and z_idx >= 0, "got, negative grid indexall indices should be positive: x_idx: %d, y_idx: %d, z_idx: %d" % (x_idx, y_idx, z_idx)
        #assing grid indices
        self.x_idx = x_idx
        self.y_idx = y_idx
        self.z_idx = z_idx
    def get_contents(self):
        return "%d %d %f %f %f %d %d %d %f %f %f" %(self.id, self.type, self.x, self.y, self.z, self.x_idx, self.y_idx, self.z_idx, self.x_eq, self.y_eq, self.z_eq)
    def shift_coordinates(self, dx, dy, dz ):
        """shift coordinates by dx, dy, dz"""
        self.x -= dx
        self.y -= dy
        self.z -= dz
        self.x_eq -= dx
        self.y_eq -= dy
        self.z_eq -= dz
    def update_coordinates(self, new_x, new_y, new_z):
        """Update coordinates of the atom from (x, y, z) -> (new_x, new_y, new_z)"""
        self.x = new_x
        self.y = new_y
        self.z = new_z
        self.x_eq = new_x
        self.y_eq = new_y
        self.z_eq = new_z


class AtomicForces(Atom):
    def __init__(self, a_id, mol, type, x, y, z, fx, fy, fz):
        Atom.__init__(self, a_id, type, x, y, z)
        self.mol = mol
        self.fx, self.fy, self.fz = fx, fy, fz
        self.radius               = self.get_radius(x, y, 0)
        self.atr                  = {}
        self.neighbors            = []
        self.ext_ids              = set()
        self.comp_ids             = set()

    def get_radius(self, x, y, z):
        return (x**2 + y**2 + z**2)**(1/2)
    def __lt__(self, other):
        return self.radius < other.radius
    def __eq__(self, other):
        return self.radius == other.radius   

def main():
    a = Atom(1, 1, 1, 2, 3)
    a.set_equilibrium(1, 2, 3)
    a.set_grid_indices(0, 0, 1)
    print(a.get_contents())

if __name__ == "__main__":
    main()

