
class Bond:
    def __init__(self, coeff1, coeff2, coeff3):
        self.coeff1 = coeff1
        self.coeff2 = coeff2
        self.coeff3 = coeff3
        self.bonded_atoms = set()

    def add_bonded_atoms(self, atom_id1, atom_id2):
        if atom_id1 > atom_id2:
            atom_id1, atom_id2 = atom_id2, atom_id1
        self.bondef_atoms.add((atom_id1, atom_id2))
