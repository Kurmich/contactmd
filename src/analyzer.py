from fcc import *
import secrets
from scipy import stats
from scipy.spatial import ConvexHull
from scipy.spatial import KDTree
forces = [500.0]
poisson = 0.5
G_shear_mod = 16.0
E_star = 1#4 * G_shear_mod
R = 1000.0
sigma = 1.0 #2**(1/6)
d = 2.0**(1.0/6.0) * sigma
atom_N = [1, 2, 3, 4]
epsilon = 0.000001

#idx_id, idx_mol, idx_type, idx_x, idx_y, idx_z, idx_fx, idx_fy, idx_fz, _ = line.split(' ')


class SimulationSettings:
    def __init__(self, force, poisson, G_shear_mod, R, sigma, d):
        self.force = force
        self.poisson = poisson
        self.G_shear_mod = G_shear_mod
        self.R = R
        self.d = d
    
class AtomicForces:
    def __init__(self, a_id, type, x, y, z, fx, fy, fz):
        self.id = a_id
        self.x, self.y, self.z = x, y, z
        self.fx, self.fy, self.fz = fx, fy, fz
        self.radius = self.get_radius(x, y, 0)
        self.attributes = {}
        self.neighbors = [] 

    def get_radius(self, x, y, z):
        return math.sqrt(x**2 + y**2 + z**2)
    def __lt__(self, other):
        return self.radius < other.radius
    def __eq__(self, other):
        return self.radius == other.radius


def contruct_neighbors(atom_forces, rc):
    '''rc - cutoff radius'''
    N = len(atom_forces)
    point2idx = {}
    data = np.zeros([N, 3])
    for i in range(N):
        af = atom_forces[i]
        x, y, z = af.x, af.y, af.z
        data[i, 0], data[i, 1], data[i, 2] = x, y, z
        point2idx[(x, y, z)] = i
    tree = KDTree(data)
    #NEED TO WRAP PBC
    for i in range(N):
        point = data[i, :]
        nbrs = tree.query_ball_point(point, rc)
        atom_forces[i].neighbors = nbrs
        #for nbr in nbrs:
            #atom_forces[i].neighbors.append(atom_forces[nbr])
            #print(data[nbr, :])
        #break
        #print(nbrs)


def add_neighbors(all_pair_ids, all_res, atype):
    N = len(all_pair_ids)
    N1 = len(all_res)
    assert N == N1
    for i in range(N):
        pair_count = 0
        pair_ids = all_pair_ids[i]
        atomic_forces = all_res[i]
        M = len(pair_ids)
        for (id1, id2) in pair_ids:
            if id1 > M or id2 > M: continue
            af1 = atomic_forces[atype][id1-1] ## ASSUMING ATOMS ARE SORTED ACCORDING TO THEIR IDS
            af2 = atomic_forces[atype][id2-1]
            print(af1.id, id1, af2.id, id2)
            assert af1.id == id1 and af2.id == id2
            af1.neighbors.append(af2)
            af2.neighbors.append(af1)
            pair_count += 1
        print("i: %d number of pairs: %d" %(i, pair_count))
def binary_search_up(atom_forces, max_r):
    """Returns the index of atom located at radius closest to max_r (i.e. atom.radius <= max_r)"""
    lo = 0
    hi = len(atom_forces)
    while hi - lo > 1:
        mid = (hi + lo)//2
        if atom_forces[mid].radius <= max_r:
            lo = mid
        else:
            hi = mid
    return lo


def binary_search_low(atom_forces, min_r):
    """Returns the index of atom located at radius closest to max_r (i.e. atom.radius <= max_r)"""
    lo = -1
    hi = len(atom_forces)-1
    while hi - lo > 1:
        mid = (hi + lo)//2
        if atom_forces[mid].radius >= min_r:
            hi = mid
        else:
            lo = mid
    return hi


def get_pair_interactions(filename, t_start, t_end):
    assert t_start <= t_end
    r_time = False
    r_entry_count = False
    r_boundary = False
    r_entries = False
    skip = False
    d = 0
    t = 0
    all_pair_ids = []
    res = []
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_entry_count:
                count = int(line)
                res = []
                print("# of pair ids: %d" %count)
                r_entry_count = False
            elif r_boundary:
                d -= 1
                r_boundary = False if d == 0 else True
            elif r_entries:
                if 'ITEM: TIMESTEP' in line:
                    r_entries = False
                else:
                    #print("reading atoms")
                   # print(len(line.split(' ')))
                    if skip:
                        continue
                    line = line.strip()
                    id1, id2  = line.split(' ')
                    id1, id2 = int(id1), int(id2)
                    res.append((id1, id2))
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if len(res) != 0:
                    all_pair_ids.append(res)
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                print("Next is timestep")
            elif 'ITEM: NUMBER OF ENTRIES' in line:
                r_entry_count = True
                print("Next is number of entries")
            elif 'ITEM: BOX BOUNDS pp pp mm' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ENTRIES' in line:
                r_entries = True
                print("Pair id lines are coming")
    return all_pair_ids

def get_interactions(filename, t_start, t_end, types, interacting = False):
    assert t_start <= t_end
    r_time = False
    r_atom_count = False
    r_boundary = False
    r_atoms = False
    skip = False
    d = 0
    t = 0
    max_r = 0
    all_res = []
    res = dict()
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                res = dict()
                print("# of atoms: %d" %count)
                r_atom_count = False
            elif r_boundary:
                d -= 1
                r_boundary = False if d == 0 else True
            elif r_atoms:
                if 'ITEM: TIMESTEP' in line:
                    r_atoms = False
                else:
                    #print("reading atoms")
                   # print(len(line.split(' ')))
                    if skip:
                        continue
                    line = line.strip()
                    #print(line)
                    #print(line.split(' ', 10))
                    '''Return a list of the words in the string, using sep as the delimiter string. 
                    If maxsplit is given, at most maxsplit splits are done (thus, the list will have at most maxsplit+1 elements). 
                    If maxsplit is not specified or -1, then there is no limit on the number of splits (all possible splits are made).'''
                    a_id, mol, atype, x, y, z, fx, fy, fz, arributes  = line.split(' ', 9)
                    a_id, mol, atype, x, y, z, fx, fy, fz = int(a_id), int(mol), int(atype), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    if atype not in types: continue
                    if interacting and abs(fx) < epsilon and abs(fy) < epsilon and abs(fz) < epsilon: continue
                    #choose interacting atoms
                    radius = math.sqrt(x**2 + y**2)
                    if atype not in res: res[atype] = []
                    res[atype].append(AtomicForces(a_id, atype, x, y, z, fx, fy, fz))
                    max_r = max(max_r, radius)
                    
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if len(res) != 0:
                    for type, atom_forces in res.items():
                        l = sorted(atom_forces, key = lambda af: af.id)
                        res[type] = l
                    all_res.append(res)
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                print("Next is timestep")
            elif 'ITEM: NUMBER OF ATOMS' in line:
                r_atom_count = True
                print("Next is number of atoms")
            elif 'ITEM: BOX BOUNDS pp pp mm' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ATOMS' in line:
                r_atoms = True
                atr_line = line[line.index('S')+1:]
                atr_words = atr_line.split(' ')
                print("Atom coordinate lines are coming")
    return all_res



def get_displ_pbr(x_next, x_prev, L):
        #get displacement vector pointing form x_prev to x_next for periodic boundary condition
        #periodicity is L
        sign = np.sign(x_next - x_prev)
        dx = abs(x_next - x_prev)
        if dx > L/2:
            dx = -(L - dx)
        return sign * dx

def append_bondlens(filename, types, chain_count, chain_len):
    N_atoms = chain_len * chain_count
    newfilename =  filename
    Lx, Ly, Lz = 0 , 0 , 0
    new_file = open("wblen.txt", "w+")
    epsilon = 0.0000000000001
    r_time = False
    r_atom_count = False
    r_boundary = False
    r_atoms = False
    skip = False
    d = 0
    t = 0
    max_r = 0
    all_res = []
    cur_frame = None
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if not r_atom_count and not r_atoms:
                if 'ITEM: ATOMS' in line:
                    words = line.split()
                    newline = ""
                    for word in words:
                        newline += (word + " ")
                    newline += "bondlen\n"
                    line = newline
                new_file.write(line)
            
                
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                cur_frame = np.zeros([N_atoms, 10])
                new_file.write(str(N_atoms) + "\n")
                print("# of atoms: %d" %count)
                r_atom_count = False
            elif r_boundary:
                words = line.split()
                L = float(words[1]) - float(words[0])
                if d == 3: 
                    Lx = L
                elif d == 2:
                    Ly = L
                elif d == 1:
                    Lz = L
                print(Lx, Ly, Lz)
                d -= 1
                r_boundary = False if d == 0 else True
            elif r_atoms:
                if 'ITEM: TIMESTEP' in line:
                    r_atoms = False
                else:
                    #print("reading atoms")
                   # print(len(line.split(' ')))
                    id, mol, type, x, y, z, fx, fy, fz, _  = line.split(' ')
                    id, mol, type, x, y, z, fx, fy, fz = int(id), int(mol), int(type), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    if fx < epsilon: fx = 0
                    if fy < epsilon: fy = 0
                    if fz < epsilon: fz = 0
                    if type not in types: continue
                    cur_frame[id-1, 0], cur_frame[id-1, 1], cur_frame[id-1, 2] = id, mol, type
                    cur_frame[id-1, 3], cur_frame[id-1, 4], cur_frame[id-1, 5] = x, y, z
                    cur_frame[id-1, 6], cur_frame[id-1, 7], cur_frame[id-1, 8] = fx, fy, fz
                    
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if cur_frame is not None:
                    #mod_atoms = ""
                    for i in range(chain_count):
                        s = i * chain_len
                        e = s + chain_len
                        chain_frame = cur_frame[s:e, :]
                        #drs = np.sqrt(np.sum(np.square(np.diff(chain_frame[:, 3:6], axis = 0)), axis = 1))
                        #chain_frame[0:chain_len-1, 9] = drs
                        for j in range(chain_len - 1):
                            dx = get_displ_pbr(chain_frame[j, 3], chain_frame[j+1, 3], Lx)
                            dy = get_displ_pbr(chain_frame[j, 4], chain_frame[j+1, 4], Ly)
                            #dz assumes that there shouldn't be problems with prc since it isn't periodic in z direction + tip is added
                            dz = get_displ_pbr(chain_frame[j, 5], chain_frame[j+1, 5], Lz) 
                            dr = math.sqrt(dx*dx + dy*dy + dz*dz)
                            chain_frame[j, 9] = dr
                        chain_frame[chain_len - 1, 9] = chain_frame[chain_len - 2, 9]
                        #chainstring = ""
                        for j in range(chain_len):
                            cur_l = ""
                            for k in range(3):
                                cur_l += str(int(chain_frame[j, k]))
                                cur_l += " "
                            for k in range(3, 10):
                                cur_l += str(chain_frame[j, k])
                                cur_l += " "
                            new_file.write(cur_l + "\n")
                    new_file.write(line)
                            #chainstring += cur_l + "\n"
                        #mod_atoms += chainstring
                    #print(mod_atoms)
                r_time = True
                t += 1
                '''if t == 3:
                    new_file.flush()
                    return
                '''
            elif 'ITEM: NUMBER OF ATOMS' in line:
                r_atom_count = True
                print("Next is number of atoms")
            elif 'ITEM: BOX BOUNDS' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ATOMS' in line:
                r_atoms = True
                print("Atom coordinate lines are coming")
        new_file.flush()
        new_file.close()



def get_frames(filename, t_start, t_end, types = None):
    assert t_start <= t_end
    r_time = False
    r_atom_count = False
    r_boundary = False
    r_atoms = False
    skip = False
    d = 0
    t = 0
    max_r = 0
    cur_frame = None
    real_count = 0
    frames = []
    #type_count = {}
    #for type in types:
     #   type_count[type] = 0
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                if not skip:
                    cur_frame = np.zeros([count, 8])
                print("# of atoms: %d" %count)
                r_atom_count = False
            elif r_boundary:
                d -= 1
                r_boundary = False if d == 0 else True
            elif r_atoms:
                if 'ITEM: TIMESTEP' in line:
                    r_atoms = False
                else:
                    #print("reading atoms")
                   # print(len(line.split(' ')))
                    if skip:
                        continue
                        
                    id, mol, type, x, y, z, fx, fy, fz, _  = line.split(' ')
                    id, mol, type, x, y, z, fx, fy, fz = int(id), int(mol), int(type), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    #if types is not None and type in types and time == t_start:
                    #    typecount[type] += 1
                    #if type not in types: continue
                    cur_frame[id-1, 0], cur_frame[id-1, 1] = mol, type
                    cur_frame[id-1, 2], cur_frame[id-1, 3], cur_frame[id-1, 4] = x, y, z
                    cur_frame[id-1, 5], cur_frame[id-1, 6], cur_frame[id-1, 7] = fx, fy, fz
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if cur_frame is not None and not skip:
                    frames.append(cur_frame)
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                print("Next is timestep")
            elif 'ITEM: NUMBER OF ATOMS' in line:
                r_atom_count = True
                print("Next is number of atoms")
            elif 'ITEM: BOX BOUNDS pp pp mm' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ATOMS' in line:
                r_atoms = True
                print("Atom coordinate lines are coming")
    return  frames


def broken_bonds(frames, type, chain_count, chain_len, max_bond_length):
    broken_bond_count = []
    bond_count = chain_len - 1
    for frame in frames:
        frame = frame[frame[:, 1] == type, :]
        print(frame)
        #print(np.diff(frame[:, 2:5], axis = 0))
        print(np.sum(np.square(np.diff(frame[:, 2:5], axis = 0)), axis = 1))
        drs = np.sqrt(np.sum(np.square(np.diff(frame[:, 2:5], axis = 0)), axis = 1))
        print(drs)
        count = 0
        for i in range(chain_count):
            print(i)
            s = i * chain_len
            e = s + chain_len
            #print(s, e)
            drs = np.sqrt(np.sum(np.square(np.diff(frame[s:e, 2:5], axis = 0)), axis = 1))
            c = drs[drs[s:e] > max_bond_length]
            count += len(drs[drs[s:e] > max_bond_length])
        broken_bond_count.append(count)
    ts = [i+1 for i in range(len(broken_bond_count))]
    plt.plot(ts, broken_bond_count)
    plt.show()



def get_avg_pressure(atom_forces, r):
    count = len(atom_forces)
    normal_pressures = np.zeros(count)
    i = 0
    for atom_f in atom_forces:
        normal_pressures[i] = atom_f.fz/(math.pi*((d/2)**2))
        i += 1
    if len(normal_pressures) == 0:
        print("NO ATOMS TO COMPUTE AVG PRESSURE")
        return 0
    total_pressure = np.sum(normal_pressures)
    print("Radius r: %f total pressure: %f number of atoms: %d" %(r, total_pressure, count))
    return abs(total_pressure/count) #IS ABS VALUE ALWAYS VALID?

def plot_avg_pressure(filename, type):
    res, frames = get_interactions(filename, 21, 22)
    atom_forces = res[type]
    max_r = atom_forces[-1].radius
    r_limit = max_r * math.cos(math.pi/4)
    bins = np.arange(sigma/2, r_limit + sigma/4, sigma )
    #calculate hertz predictions
    hertz_max_r = (3*forces[-1]*R/(4*E_star))**(1/3) #adhoc here on choice of forces
    hertz_bins = np.arange(0, hertz_max_r, hertz_max_r/100)
    hertz_bins = np.append(hertz_bins, hertz_max_r)
    hertz_pres = [(2*hertz_max_r/(math.pi*R)) * math.sqrt(1 - (r/hertz_max_r)**2) for r in hertz_bins]
    #calculate average pressures
    avg_pressures = []
    for bin_r in bins:
        low_bound = bin_r - sigma/2
        up_bound = bin_r + sigma/2
        low_idx = binary_search_low(atom_forces, low_bound)
        up_idx  = binary_search_up(atom_forces, up_bound)
        atom_force_bin = atom_forces[low_idx:(up_idx+1)]
        if bin_r == 45.5:
            print(get_avg_pressure(atom_forces[0:low_idx+1], 41.5))
            print(get_avg_pressure(atom_forces[low_idx:(len(atom_forces))], 41.5))
        avg_pressures.append(get_avg_pressure(atom_force_bin, bin_r))
        print(bin_r, low_idx, up_idx)
    new_avg_pressures = [p/E_star for p in avg_pressures]
    print(hertz_bins.size)
    print("total pressure: %f" %(sum(avg_pressures)))
    #plt.plot(hertz_bins, hertz_pres, 'r')
    plt.plot(bins, new_avg_pressures, 'bo')
    plt.suptitle('Normal pressure distribution in contact zone. Normal load: %d' %(forces[-1]), fontsize = 20) #adhoc here on choice of forces
    plt.ylabel(r'$p/E^{\ast}$', fontsize = 16)
    plt.xlabel(r'$r/\sigma$', fontsize = 16)
    plt.show()


def plot_nforce_vs_cont_area():
    max_radii = []
    t1 = np.arange(0.0, 0.06, 0.01)
    for force in forces:
        filename = 'visualize_%d.out' %force
        res = get_interactions(filename)
        max_r = 0
        for type, atom_forces in res.items():
            max_r = max(max_r, atom_forces[-1].radius)
        max_radii.append(max_r)
    new_forces = [(3*force/(4*E_star*R**2))**(1/3) for force in forces]
    new_radii = [r/R for r in max_radii]
    print(t1)
    plt.plot(t1, t1, 'r')
    plt.plot(new_forces, new_radii, 'bo')
    plt.suptitle('Contact radius a vs normal load N', fontsize = 20)
    plt.ylabel(r'$a/R$', fontsize = 16)
    plt.xlabel(r'$(3N/4E^{\ast}R^2)^{1/3}$', fontsize = 16)
    plt.show()

def plot_layer_density(l_type, frames, t):
    frame_l = select_from_frame(frames[t-1], l_type)
    z_idx = 4
    z_vals = frame_l[:, z_idx]
    print(frame_l)
    print(z_vals)
    hist, bin_edges = np.histogram(z_vals, bins = 100)
    bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])
    fig, ax = plt.subplots()
    plt.xlabel("z/sigma")
    plt.title("Number of atoms vs z.")
    ax.plot(bincenters, hist, '-', marker = 'o', fillstyle = 'none', markerfacecolor = 'r')
    #ax.hist(z_vals, bins='auto')
    ax.legend()
    plt.show()

def plot_color_map(data):
    print(data)
    fig, ax = plt.subplots()
    ax.pcolormesh(data)
    #fig = plt.figure()
    #ax = Axes3D(fig)
    #ax.scatter(data[:, 0], data[:, 1], data[:, 2])
    #ax.set_title('pcolormesh')
    plt.ylabel("S")
    plt.xlabel("log(t/tau)")
    plt.show()


def get_displacements(type, frames, t1, t2):
    frame1 = select_from_frame(frames[ t1 - 1 ], type)
    frame2 = select_from_frame(frames[ t2 - 1 ], type)
    delta = frame2[:, 2:5] - frame1[:, 2:5]
    nRows, nCols = delta.shape
    #x = np.arange(0, nRows, 1)
    #xx, yy = np.meshgrid(x, x, sparse=True)
    #print(xx)
    delta_d = np.sqrt( np.sum(np.square(delta) , axis = 1))
    cmap = frame1[:, 2:5]
    cmap[:, 2] = delta_d
    print(delta_d.shape, frame1.shape)
    return cmap
    

def select_from_frame(frame, type):
    '''Get given type of molecules from set of all molecules'''
    idx = frame[:, 1] == type
    return frame[idx, :]


def print_total_load(frame, type):
    frame = select_from_frame(frame, type)
    print(np.sum(frame[:, 7]))



def interacting_particles(atom_forces):
    interacting = []
    for af in atom_forces:
        if abs(af.fx) > epsilon or abs(af.fy) > epsilon or abs(af.fz) > epsilon:
            interacting.append(af)
    N = len(interacting)
    points = np.zeros([N, 2])
    for i in range(N):
        points[i, 0], points[i, 1] = interacting[i].x, interacting[i].y
    return interacting, points


def shoelace_area(xs, ys):
    '''Let 'vertices' be an array of N pairs (x,y), indexed from 0
        Let 'area' = 0.0
        for i = 0 to N-1, do
          Let j = (i+1) mod N
          Let area = area + vertices[i].x * vertices[j].y
          Let area = area - vertices[i].y * vertices[j].x
        end for'''
    N = len(xs)
    area = 0
    for i in range(N):
        j = (i+1) % N
        area += (xs[i] * ys[j] - ys[i] * xs[j])
    return abs(area)/2

def orientation(atom_f1, atom_f2, atom_f3):
    #http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    val = (atom_f2.y - atom_f1.y) * (atom_f3.x - atom_f2.x) - (atom_f3.y - atom_f2.y) * (atom_f2.x - atom_f1.x)
    if val == 0: #collinear
        return 0
    elif val > 0: #clockwise
        return 1
    else:
        return 2 # counterlockwise
def further(atom_f1, atom_f2, atom_f3):
    d1 = (atom_f2.y - atom_f1.y) * (atom_f2.y - atom_f1.y) + (atom_f2.x - atom_f1.x) * (atom_f2.x - atom_f1.x)
    d2 = (atom_f3.y - atom_f1.y) * (atom_f3.y - atom_f1.y) + (atom_f3.x - atom_f1.x) * (atom_f3.x - atom_f1.x)
    if d2 > d1:
        return 1 #return 1 if atom 2 is further than atom 1
    return -1

def jarvis(atom_forces, visualize = False):
    leftmost = atom_forces[0]
    for af in atom_forces:
        if af.x < leftmost.x:
            leftmost = af
    hull = [leftmost]
    max_r = 0
    cur, cand = leftmost, None
    xs, ys = [], []
    remaining = [af for af in atom_forces if af not in hull]
    while True:
        cand = secrets.choice(atom_forces)
        for af in atom_forces:
            if af != cand:
                if orientation(cur, cand, af) == 2:
                    #print("here")
                    cand = af
                elif orientation(cur, cand, af) == 0 and further(cur, cand, af) == 1:
                    cand = af
        if cand == leftmost: break
        cur = cand
        max_r = max(max_r, cur.radius)
        xs.append(cur.x)
        ys.append(cur.y)
        hull.append(cand)

    area = shoelace_area(xs, ys)

    if visualize:
        xs, ys = [], []
        for af in hull:
            xs.append(af.x)
            ys.append(af.y)
        print("vis area %g" %(shoelace_area(xs, ys)))
        xs.append(hull[0].x)
        ys.append(hull[0].y)
        plt.plot(xs, ys, color = 'r')
        xs, ys = [], []
        for af in atom_forces:
            xs.append(af.x)
            ys.append(af.y)
        plt.scatter(xs, ys)
        plt.show()
    return hull, max_r, area
#def get_edge_avg(atom_forces):


def get_total_forces(atomic_forces):
    fx, fy, fz = 0, 0, 0
    for af in atomic_forces:
        fx += af.fx
        fy += af.fy
        fz += af.fz
    return fx, fy, fz


def get_contact_depth(interacting_af):
    zlo, zhi = 100000000000, -100000000000
    for af in interacting_af:
        zlo = min(zlo, af.z)
        zhi = max(zhi, af.z)
    return zhi - zlo

def plot_stresszz_d(all_res, type):
    areas = []
    times = []
    ds  = []
    strs_z = []
    fzs = []
    num_intrs = []
    dt = 0.01 * 500000
    v = 0.0001
    d_start = 0
    t0 = 0
    d0 = 0
    first_contact = True
    ts = 0
    avg_strs, cnt = 0, 0
    fz_prev = 0
    E_modulus = 0
    for t in range(len(all_res)):
        print("New timestep: %d" %t)
        res = all_res[t]
        interacting_af, points = interacting_particles(res[type])
        count = len(interacting_af)
        if count == 0:
            continue
        elif count < 3:
            #continue
            if first_contact:
                t0 = t-1
                print("First contact happens at t0 : %d" %t0)
                first_contact = False
           # continue
            fx, fy, fz = get_total_forces(interacting_af)
            area = count * math.pi * (d/2)**2
            hc = get_contact_depth(interacting_af)
        else:
            #if first_contact:
            #    t0 = t-1
            #    first_contact = False
            #_, _, all_fz = get_total_forces(res[type])
            fx, fy, fz = get_total_forces(interacting_af)
            hull = ConvexHull(points)
            #plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
            #plt.show()
            print("scipy area: %g" %hull.area)
            area = shoelace_area(points[hull.vertices,0], points[hull.vertices,1])
            hc = get_contact_depth(interacting_af)
            #hull, max_r, area = jarvis(interacting_af, visualize=True)
            #area = math.pi * max_r**2
        

        if fz_prev == 0: fz_prev = fz
        if fz_prev != fz:
            del_fz    = fz - fz_prev
            del_h     = v * dt
            stiffness = del_fz/del_h
            fz_prev   = fz
            E_modulus = (stiffness * (math.pi**(1/2))) / (2 * area**(1/2))
        
        del_z = (t - t0) * v * dt
        if del_z < d0:
            continue
        else:
            del_z -= d0
        str_z = fz / area
        #if str_z > 3:
        #    continue
        areas.append(area)
        fzs.append(fz)
        num_intrs.append(count)
        times.append(t)
        ds.append(del_z)
        strs_z.append(str_z)
        if del_z > 8 and del_z < 12:
            avg_strs += str_z
            cnt += 1
        print("Displacement d: %g Contact Depth: %g Num of particles: %d Total Fz: %g Area: %g StressZ: %g E: %g" %(del_z,hc,count,fz, area, str_z, E_modulus))
    fz_polcoeffs = np.polyfit(ds, fzs, 2)
    az_polcoeffs = np.polyfit(ds, areas, 2)
    print(fz_polcoeffs)
    print(az_polcoeffs)
    polfz = np.poly1d(fz_polcoeffs)
    fig, ax = plt.subplots(2, 2)
    E_moduli = [2 * fzs[i] * math.tan(math.radians(45))/ areas[i] for i in range(len(ds)) ]
    ax[0][0].plot(ds, num_intrs, 'g')
    ax[0][0].set_ylabel("N", rotation = 0)
    ax[0][1].plot(ds, areas, 'r')
    ax[0][1].set_ylabel("A", rotation = 0)
    ax[1][0].plot(ds, fzs, 'k')
    ax[1][0].set_ylabel("$F_z$", rotation = 0)
    ax[1][1].plot(ds, strs_z, 'b')
    ax[1][1].set_ylabel("$\sigma_{zz}$", rotation = 0)
    for i in range(2):
        for j in range(2):
            ax[i][j].set_xlabel("$d$")
    fig.suptitle("Shift d = %g" %d0)
    print("Avg stress: %g" %(avg_strs/cnt))
    plt.show()
    '''
    plt.plot(ds, fzs, label = "F_z vs d")
    ds2 = [d**2 for d in ds]
    plt.plot(ds2, fzs)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ds2,fzs)
    print("slope: %g intercept: %g r_val: %g p_val: %g std_err: %g" %(slope, intercept, r_value, p_value, std_err))
    plt.ylabel('$F_z$')
    plt.xlabel('$d^2$')
    plt.plot(ds, polfz(ds))
    plt.legend()
    plt.show()
    plt.plot(ds, strs_z)
    plt.show()
    print(E_moduli)
    '''

def main():
    #plot_nforce_vs_cont_area()
    #substrate_type = 1
    #tip_type = 2
    #oligomer_type = 3
    M, N = 2000, 256
    r = 10
    cang = 45
    tip_type = 2
    glass = 1
    t_start = 1
    t_end = 5
    types = [glass]
    rc= 1.5
    #bond testing
    #filename = '../visfiles/viscomp_M%d_N%d.out' %(M, N)
    #frames = get_frames(filename, t_start, t_end)
    #broken_bonds(frames, glass, M, N, 1.5)
    #return
    filename = '../visfiles/visualize_M%d_N%d_r%d_cang%d.out' %(M, N, r, cang)
    filenameinteractions = '../visfiles/pairids_M%d_N%d_r%d_cang%d.out' %(M, N, r, cang)
    #append_bondlens(filename, types, M, N)
    #return
    all_res = get_interactions(filename, t_start, t_end, types, interacting = True)
    all_inter = get_pair_interactions(filenameinteractions, t_start, t_end)
    add_neighbors(all_inter, all_res, glass)
    return
    #atom_forces = all_res[glass][1]
    #add_neighbors(atom_forces, rc)
#    plot_layer_density(glass, frames, t_start + 1)
 #   plot_layer_density(glass, frames, t_start + 20)
    #print_total_load(frames[0], substrate_type)
    #print_total_load(frames[0], tip_type)
    #data = get_displacements(2, frames, 1, 2)
    #plot_color_map(data)
    #plot_avg_pressure(filename, substrate_type)
    #plot_avg_pressure(filename, tip_type)
    #frame = frames[t_end-1]
    #indices = np.where(frame[:, 1] == tip_type)
    #tipframe = frame[indices, :]
    #interacting_af = interacting_particles(res[tip_type])
    #jarvis(interacting_af, visualize = True)
    plot_stresszz_d(all_res, tip_type)
    '''square = np.zeros([4, 2])
    square[1, 0] = 1
    square[2, 1] = 1
    square[3, 0], square[3, 1] = 1, 1
    hull = ConvexHull(square)
    plt.plot(square[hull.vertices,0], square[hull.vertices,1], 'r--', lw=2)
    plt.show()
    print(hull.area)'''
    print("M: %d N: %d r: %d cang: %d t_start: %d t_end: %d" %(M, N, r, cang, t_start, t_end))

if __name__=="__main__":
    main()

#plt.plot(timestep, forcesx)
#plt.suptitle('Friction force vs timestep', fontsize = 20)
#plt.ylabel('Fx', fontsize = 20)
#plt.xlabel('t', fontsize = 20)
#plt.axhline(0, color = 'black')
#plt.show()
