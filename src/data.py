# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# data tool

from os import popen
from copy import deepcopy
oneline = "Read, write, manipulate LAMMPS data files"

docstr = """
d = data("data.poly")            read a LAMMPS data file, can be gzipped
d = data()			 create an empty data file

d.map(1,"id",3,"x")              assign names to atom columns (1-N)

coeffs = d.get("Pair Coeffs")    extract info from data file section
q = d.get("Atoms",4)

  1 arg = all columns returned as 2d array of floats
  2 args = Nth column returned as vector of floats

d.reorder("Atoms",1,3,2,4,5)     reorder columns (1-N) in a data file section

  1,3,2,4,5 = new order of previous columns, can delete columns this way

d.title = "My LAMMPS data file"	 set title of the data file
d.headers["atoms"] = 1500        set a header value
d.sections["Bonds"] = lines      set a section to list of lines (with newlines)
d.delete("bonds")		 delete a keyword or section of data file
d.delete("Bonds")
d.replace("Atoms",5,vec)      	 replace Nth column of section with vector
d.newxyz(dmp,1000)		 replace xyz in Atoms with xyz of snapshot N

  newxyz assumes id,x,y,z are defined in both data and dump files
    also replaces ix,iy,iz if they are defined
  
index,time,flag = d.iterator(0/1)          loop over single data file snapshot
time,box,atoms,bonds,tris,lines = d.viz(index)   return list of viz objects

  iterator() and viz() are compatible with equivalent dump calls
  iterator() called with arg = 0 first time, with arg = 1 on subsequent calls
    index = timestep index within dump object (only 0 for data file)
    time = timestep value (only 0 for data file)
    flag = -1 when iteration is done, 1 otherwise
  viz() returns info for specified timestep index (must be 0)
    time = 0
    box = [xlo,ylo,zlo,xhi,yhi,zhi]
    atoms = id,type,x,y,z for each atom as 2d array
    bonds = id,type,x1,y1,z1,x2,y2,z2,t1,t2 for each bond as 2d array
      NULL if bonds do not exist
    tris = NULL
    lines = NULL
    
d.write("data.new")             write a LAMMPS data file
"""

# History
#   8/05, Steve Plimpton (SNL): original version
#   11/07, added triclinic box support
#   05/21/2019 Porting to python3
# ToDo list

# Variables
#   title = 1st line of data file
#   names = dictionary with atom attributes as keys, col #s as values
#   headers = dictionary with header name as key, value or tuple as values
#   sections = dictionary with section name as key, array of lines as values
#   nselect = 1 = # of snapshots

# Imports and external programs


try:
    tmp = PIZZA_GUNZIP
except:
    PIZZA_GUNZIP = "gunzip"

# Class definition


class data:

    # --------------------------------------------------------------------

    def __init__(self, *list):
        self.nselect = 1

        if len(list) == 0:
            self.title = "LAMMPS data file"
            self.names = {}
            self.headers = {}
            self.sections = {}
            return

        filename = list[0]
        if filename[-3:] == ".gz":
            file = popen("%s -c %s" % (PIZZA_GUNZIP, filename), 'r')
        else:
            file = open(filename)

        self.title = file.readline()
        #print("Title: %s" %self.title)
        self.names = {}
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
            #print("current line %s" %line)
            if line.isspace():
                line = file.readline()
                continue
            for pair in skeywords:
                keyword, length = pair[0], pair[1]
                if keyword in line:
                    found = True
                    if length not in headers.keys():
                        raise Exception("data section %s has no matching header value" % line)
                    file.readline()
                    list = []
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


        file.close()
        self.headers = headers
        self.sections = sections

    def append(self, other, sep_z):
        """Append contents of other data file to this file"""
        #other.sections["Velocities"], self.sections["Velocities"] = None, None
        #del other.sections["Velocities"]
        #del self.sections["Velocities"]
        z_idx = 5
        type_idx = 2
        delta_z = 0
        print("Check index for z coordinates of atoms currenctly using %d" %z_idx)
        old_headers = deepcopy(self.headers)
        for keyword in hkeywords:
            if  keyword in self.headers.keys() or keyword in other.headers.keys():
                if keyword == "xlo xhi" or keyword == "ylo yhi":
                    pair = self.headers[keyword]
                    o_pair = other.headers[keyword]
                    #print(pair[0], o_pair[0])
                    p0 = min(pair[0], o_pair[0])
                    p1 = max(pair[1], o_pair[1])
                    self.headers[keyword] = (p0, p1)
                elif keyword == "zlo zhi":
                    pair = self.headers[keyword]
                    o_pair = other.headers[keyword]
                    delta_z = pair[1] - o_pair[0]
                    self.headers[keyword] = (pair[0], o_pair[1] + delta_z + sep_z)
                elif keyword == "xy xz yz":
                    triple = self.headers[keyword]
                else:
                    if keyword not in other.headers.keys():
                        other.headers[keyword] = 0
                    if keyword not in self.headers.keys():
                        self.headers[keyword] = 0
                    self.headers[keyword] += other.headers[keyword]
        atom_count       = old_headers["atoms"]
        atom_type_count  = old_headers["atom types"]
        bond_type_count  = old_headers["bond types"]
        angle_type_count = 0 if "angle types" not in old_headers else old_headers["angle types"] 
        for pair in skeywords:
            keyword = pair[0]
            if keyword in other.sections.keys():
                start_id = old_headers[pair[1]]
                for line in other.sections[keyword]:
                    words = line.split()
                    id = int(words[0])
                    id += start_id
                    words[0] = str(id)
                    if keyword == "Atoms":
                        words[z_idx] = str( float(words[z_idx]) + delta_z + sep_z )
                        words[type_idx] = str(int(words[type_idx]) + atom_type_count)
                    if keyword == "bonds" or keyword == "Bonds":
                        words[1] = str( int(words[1]) +  bond_type_count)
                        words[2] = str( int(words[2]) +  atom_count)
                        words[3] = str( int(words[3]) +  atom_count)
                    if keyword == "angles" or keyword == "Angles":
                        words[1] = str( int(words[1]) +  angle_type_count)
                        words[2] = str( int(words[2]) +  atom_count)
                        words[3] = str( int(words[3]) +  atom_count)
                        words[4] = str( int(words[4]) +  atom_count)
                    newline = ' '.join(words)
                    self.sections[keyword].append(newline)
                    self.sections[keyword].append('\n')


    # --------------------------------------------------------------------
    # assign names to atom columns

    def map(self, *pairs):
        if len(pairs) % 2 != 0:
            raise StandardError( "data map() requires pairs of mappings")
        for i in range(0, len(pairs), 2):
            j = i + 1
            self.names[pairs[j]] = pairs[i]-1

    # --------------------------------------------------------------------
    # extract info from data file fields

    def get(self, *list):
        if len(list) == 1:
            field = list[0]
            array = []
            lines = self.sections[field]
            for line in lines:
                words = line.split()
                values = map(float, words)
                array.append(values)
            return array
        elif len(list) == 2:
            field = list[0]
            n = list[1] - 1
            vec = []
            lines = self.sections[field]
            for line in lines:
                words = line.split()
                vec.append(float(words[n]))
            return vec
        else:
            raise StandardError( "invalid arguments for data.get()")

    # --------------------------------------------------------------------
    # reorder columns in a data file field

    def reorder(self, name, *order):
        n = len(order)
        natoms = len(self.sections[name])
        oldlines = self.sections[name]
        newlines = natoms*[""]
        for index in order:
            for i in xrange(len(newlines)):
                words = oldlines[i].split()
                newlines[i] += words[index-1] + " "
        for i in xrange(len(newlines)):
            newlines[i] += "\n"
        self.sections[name] = newlines

    # --------------------------------------------------------------------
    # replace a column of named section with vector of values

    def replace(self, name, icol, vector):
        lines = self.sections[name]
        newlines = []
        j = icol - 1
        for i in xrange(len(lines)):
            line = lines[i]
            words = line.split()
            words[j] = str(vector[i])
            newline = ' '.join(words) + '\n'
            newlines.append(newline)
        self.sections[name] = newlines

    # --------------------------------------------------------------------
    # replace x,y,z in Atoms with x,y,z values from snapshot ntime of dump object
    # assumes id,x,y,z are defined in both data and dump files
    # also replaces ix,iy,iz if they are defined

    def newxyz(self, dm, ntime):
        nsnap = dm.findtime(ntime)

        dm.sort(ntime)
        x, y, z = dm.vecs(ntime, "x", "y", "z")

        self.replace("Atoms", self.names['x']+1, x)
        self.replace("Atoms", self.names['y']+1, y)
        self.replace("Atoms", self.names['z']+1, z)

        if dm.names.has_key("ix") and self.names.has_key("ix"):
            ix, iy, iz = dm.vecs(ntime, "ix", "iy", "iz")
            self.replace("Atoms", self.names['ix']+1, ix)
            self.replace("Atoms", self.names['iy']+1, iy)
            self.replace("Atoms", self.names['iz']+1, iz)

    # --------------------------------------------------------------------
    # delete header value or section from data file

    def delete(self, keyword):

        if self.headers.has_key(keyword):
            del self.headers[keyword]
        elif self.sections.has_key(keyword):
            del self.sections[keyword]
        else:
            raise StandardError("keyword not found in data object")

    # --------------------------------------------------------------------
    # write out a LAMMPS data file

    def write(self, file):
        f = open(file, "w")
        print(self.title, end="\n\n", file=f)
        for keyword in hkeywords:
            if  keyword in self.headers.keys():
                if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                        keyword == "zlo zhi":
                    pair = self.headers[keyword]
                    print(pair[0], pair[1], keyword, file=f)
                elif keyword == "xy xz yz":
                    triple = self.headers[keyword]
                    print(triple[0], triple[1], triple[2], keyword, file=f)
                else:
                    print(self.headers[keyword], keyword, file = f)
        for pair in skeywords:
            keyword = pair[0]
            if keyword in self.sections.keys():
                print("\n%s\n" % keyword, file=f)
                for line in self.sections[keyword]:
                    print(line, end = "", file=f)
        f.close()

    # --------------------------------------------------------------------
    # iterator called from other tools

    def iterator(self, flag):
        if flag == 0:
            return 0, 0, 1
        return 0, 0, -1

    # --------------------------------------------------------------------
    # time query from other tools

    def findtime(self, n):
        if n == 0:
            return 0
        raise StandardError("no step %d exists" % (n))

    # --------------------------------------------------------------------
    # return list of atoms and bonds to viz for data object

    def viz(self, isnap):
        if isnap:
            raise StandardError("cannot call data.viz() with isnap != 0")

        id = self.names["id"]
        type = self.names["type"]
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]

        xlohi = self.headers["xlo xhi"]
        ylohi = self.headers["ylo yhi"]
        zlohi = self.headers["zlo zhi"]
        box = [xlohi[0], ylohi[0], zlohi[0], xlohi[1], ylohi[1], zlohi[1]]

        # create atom list needed by viz from id,type,x,y,z

        atoms = []
        atomlines = self.sections["Atoms"]
        for line in atomlines:
            words = line.split()
            atoms.append([int(words[id]), int(words[type]),
                          float(words[x]), float(words[y]), float(words[z])])

        # create list of current bond coords from list of bonds
        # assumes atoms are sorted so can lookup up the 2 atoms in each bond

        bonds = []
        if self.sections.has_key("Bonds"):
            bondlines = self.sections["Bonds"]
            for line in bondlines:
                words = line.split()
                bid, btype = int(words[0]), int(words[1])
                atom1, atom2 = int(words[2]), int(words[3])
                atom1words = atomlines[atom1-1].split()
                atom2words = atomlines[atom2-1].split()
                bonds.append([bid, btype,
                              float(atom1words[x]), float(atom1words[y]),
                              float(atom1words[z]),
                              float(atom2words[x]), float(atom2words[y]),
                              float(atom2words[z]),
                              float(atom1words[type]), float(atom2words[type])])

        tris = []
        lines = []
        return 0, box, atoms, bonds, tris, lines

    # --------------------------------------------------------------------
    # return box size

    def maxbox(self):
        xlohi = self.headers["xlo xhi"]
        ylohi = self.headers["ylo yhi"]
        zlohi = self.headers["zlo zhi"]
        return [xlohi[0], ylohi[0], zlohi[0], xlohi[1], ylohi[1], zlohi[1]]

    # --------------------------------------------------------------------
    # return number of atom types

    def maxtype(self):
        return self.headers["atom types"]

# --------------------------------------------------------------------
# data file keywords, both header and main sections


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



def add_sphere_tip(M, N, T, r, Dx, sep_z):
    sep_z = 2**(1/6) + 0.5
    d = data("../lammpsinput/clean_quenched_stiff_M%d_N%d_T%g_smooth.data" %(M, N, T))
#    d2 = data("../lammpsinput/flattip_Dx%d.dat" %(Dx))
    d2 = data("../lammpsinput/sphere_r%d_Dx%d.dat" %(r, Dx))
    d.append(d2, sep_z)
    d.write("../lammpsinput/spheretip_data_quenched_stiff_M%d_N%d_T%g_sphR%d_smooth" %(M, N, T, r))
    print("sepz: %g" %sep_z)
    

def main():
    M = 2000
    N = 256
    T = 0.2
    r = 25#30
    cone_ang = 0#30
    Dx = 91#85
    sep_z = 2**(1/6) + 0.5
    add_sphere_tip(M, N, T, r, Dx, sep_z)
    return
    d = data("../lammpsinput/clean_quenched_stiff_M%d_N%d_T%g.data" %(M, N, T))
#    d2 = data("../lammpsinput/flattip_Dx%d.dat" %(Dx))
    d2 = data("../lammpsinput/tip_r%d_Dx%d_cang%d.dat" %(r, Dx, cone_ang))
    d.append(d2, sep_z)
    d.write("../lammpsinput/wtip_data_quenched_stiff_M%d_N%d_T%g_sphR%d_cang%d" %(M, N, T, r, cone_ang))
    print("sepz: %g" %sep_z)


if __name__ == "__main__":
    main()
