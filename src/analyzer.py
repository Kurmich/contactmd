from fcc import *
import secrets
from scipy import stats
from scipy.spatial import ConvexHull
from scipy.spatial import KDTree
import os
from thermoparser import heat_stats
from boxcell import *
import gc
forces = [500.0]
poisson = 0.5
G_shear_mod = 16.0
E_star = 1#4 * G_shear_mod
R = 1000.0
sigma = 1.0 #2**(1/6)
d = 2.0**(1.0/6.0) * sigma
atom_N = [1, 2, 3, 4]
epsilon = 0.000001
vis_data_path = '../visfiles/'
out_data_path = '../visfiles/'
#idx_id, idx_mol, idx_type, idx_x, idx_y, idx_z, idx_fx, idx_fy, idx_fz, _ = line.split(' ')
#vis_data_path = '/home/kkurman1/equilibrate/identtip/visfiles/'
#out_data_path =  '/home/kkurman1/equilibrate/identtip/outputfiles/'




class SimulationSettings:
    def __init__(self, force, poisson, G_shear_mod, R, sigma, d):
        self.force = force
        self.poisson = poisson
        self.G_shear_mod = G_shear_mod
        self.R = R
        self.d = d


class ConesimSettings:
    def __init__(self, M, N, T, r, cang, vz, dt):
        self.M = M       #number of chains
        self.N = N       #number of monomers per chain
        self.T = T       #Temperature of the system
        self.r = r       #radius of the spherical tip of the cone
        self.cang = cang #angle that cone makes with substrate plane
        self.vz = vz     #velocity of the cone
        self.dt = dt     #timestep of the simulation
        print("Simulation settings used for analysis")
        print("M : %d, N: %d, T: %g, r: %g, cang: %g, vz: %g, dt: %g" %(M, N, T, r, cang, vz, dt))
    def set_analysisvals(self, t_init, t_final, t_step):
        self.t_init  = t_init  #analysis starts from this dump 
        self.t_final = t_final #analysis continue until this dump
        self.t_step  = t_step  #us this much dumps for analysis at single iteration
        print("Analysis starts from dump time: %g continues until: %g in steps of %g" %(t_init, t_final, t_step))
        
Temp = 0.2
css = ConesimSettings(2000, 256, Temp, 10, 45, 0.0001, 0.01)
css.set_analysisvals(17, 20, 1)

  
def reconstuct_ave_lj_bonds(all_res, atype, all_bounds, percent):
    '''If average positions are dumped this method is to reconstruct'''

    N = len(all_res)
    rc = 1.5 #cut off radius
    r_lim = 2*rc + 0.5  #cutoff after which dfs search for new neighbors ends 
    for i in range(N):
        bounds = all_bounds[i]
        atom_forces = all_res[i][atype]
        new_nbrs = {}
        for af in atom_forces:
            seen_ids = []
            new_nbrs[af.id] = []
            #print("Next atom")
            for child_af in af.neighbors:
                inser_new_nbrs(af, child_af, new_nbrs[af.id], seen_ids, rc, r_lim, bounds.Lx, bounds.Ly, bounds.Lz)
        for af in atom_forces:
            af.neighbors = new_nbrs[af.id]
            
            

def inser_new_nbrs(root_af, child_af, nbrs, seen_ids, rc, r_lim, Lx, Ly, Lz):
    '''create a new neighbor list for root af such that neighbors are within <= rc'''
    if child_af.id in seen_ids: return
    dx, dy, dz = get_displ_pbr(root_af.x, child_af.x, Lx), get_displ_pbr(root_af.y, child_af.y, Ly), get_displ_pbr(root_af.z, child_af.z, Lz)
    r_cur  = math.sqrt( dx**2 + dy**2 + dz**2 )
    seen_ids.append(child_af.id)
    if r_cur > r_lim: return
    if r_cur <= rc:
        nbrs.append(child_af)
    for af in child_af.neighbors:
        inser_new_nbrs(root_af, af, nbrs, seen_ids, rc, r_lim, Lx, Ly, Lz )
        
memory_atom_forces_p = None

def get_lj_bond_stats(all_res, atype, all_bounds, delta_r, step):
    '''ASSUMPTION: that all atom_forces are sorted by their IDs'''
    breaks        = []
    formations    = []
    changes_comp  = []
    changes_ext   = []
    comp_then_ext = []
    ext_then_comp = []
    global memory_atom_forces_p
    N = len(all_res)
    assert N >= step
    for i in range(step, N):
        lj_change_count_c = 0
        lj_change_count_e = 0
        broken_count      = 0 
        formed_count      = 0
        atom_forces       = all_res[i][atype]
        atom_forces_p     = all_res[i-step][atype]
        bounds            = all_bounds[i]
        bounds_p          = all_bounds[i-step]
        M                 = len(atom_forces)
        print("Periodicities bounds ", bounds.pbcX, bounds.pbcY, bounds.pbcZ)
        print("Periodicities bounds_p ", bounds_p.pbcX, bounds_p.pbcY, bounds_p.pbcZ)
        #fig = plt.figure()
        #ax = Axes3D(fig)
        for j in range(M):
            af   = atom_forces[j]
            af_p = atom_forces_p[j]
            assert af.id == af_p.id
            lj_change_comp, lj_change_ext, broken, formed = get_stats(atom_forces, atom_forces_p, af_p, af_p.neighbors, af, af.neighbors, bounds_p, bounds, delta_r)
            lj_change_count_c  += len(lj_change_comp)
            lj_change_count_e  += len(lj_change_ext)
            broken_count       += len(broken) 
            formed_count       += len(formed)
            af_p.atr['comp']   = len(lj_change_comp)
            af_p.atr['ext']    = len(lj_change_ext)
            af_p.atr['broken'] = len(broken)
            af_p.atr['formed'] = len(formed)
            #for (af1, af2) in lj_change:
                #ax.plot([af1.x, af2.x], [af1.y, af2.y], [af1.z, af2.z])
        changes_comp.append(lj_change_count_c/2)
        changes_ext.append(lj_change_count_e/2)
        breaks.append(broken_count/2)
        formations.append(formed_count/2)
        if i == step and memory_atom_forces_p is not None:
            print("Assigning from memory")
            atom_forces_p = memory_atom_forces_p
        elif i == N-1:
            print("Assigning memory")
            memory_atom_forces_p = atom_forces
        comp_back, ext_back = reverse_changes(atom_forces_p, atom_forces)
        comp_then_ext.append(ext_back/2)
        ext_then_comp.append(comp_back/2)
        #plt.show()
        #fig.savefig("t: %d.png" %i)
        print("i: %d change comp: %g change ext: %d broken: %d formed: %d c-e: %d e-c: %d" %(i, lj_change_count_c, lj_change_count_e, broken_count, formed_count, ext_back, comp_back), flush = True)
    return changes_comp, changes_ext, breaks, formations, comp_then_ext, ext_then_comp



def reverse_changes(atom_forces_p, atom_forces):
    comp_back, ext_back = 0, 0
    for i in range(len(atom_forces)):
        af   = atom_forces[i]
        af_p = atom_forces_p[i]
        for id_ext in af.ext_ids:
            if id_ext in af_p.comp_ids:
                ext_back += 1
        for id_comp in af.comp_ids:
            if id_comp in af_p.ext_ids:
                comp_back += 1
    return comp_back, ext_back

def maxchange_criteria(r_prev, r_cur, delta_r):
    if abs(r_prev - r_cur) > delta_r:
        if r_cur > r_prev:
            return 1
        else:
            return -1
    return 0


def fracchange_criteria(r_prev, r_cur, percent):
    if abs(r_prev - r_cur) / r_prev > percent and r_cur > r_prev: 
            return 1
    elif abs(r_prev - r_cur) / r_prev > percent:
            return -1
    return 0
        

def get_stats(atom_forces, atom_forces_p, af_p, nbr_list_p, af, nbr_list, bounds_p, bounds, delta_r):
    '''bounds - current bounds
       bounds_p - previous bounds, here _p stands for previous dump
    '''

    r_cutoff = 1.501 
    lj_change_comp = []
    lj_change_ext = []
    prev_ids = {}
    cur_ids = {}
    
    for i in range(len(nbr_list_p)): prev_ids[nbr_list_p[i].id] = i 
    for i in range(len(nbr_list)):   cur_ids[nbr_list[i].id] = i 
    for nbr_prev in nbr_list_p:
        #for each neighbor in neighbor list of previous timestep
        nbr_id = nbr_prev.id
        #af_nbr_p = nbr_list_p[prev_ids[nbr_id]]
        af_nbr_p = nbr_prev
        if nbr_id in cur_ids:
            af_nbr   = nbr_list[cur_ids[nbr_id]] #if the nbr_prev is still a neighbor
        else:
            af_nbr = atom_forces[nbr_id - 1]     #if it isn't neighbor anymore
        '''ASSUMING NON-PBC in z direction'''
        dx, dy, dz = get_dxdydz(af_p, af_nbr_p, bounds_p)
        r_prev = math.sqrt( dx**2 + dy**2 + dz**2 )
        
        dx, dy, dz = get_dxdydz(af, af_nbr, bounds)
        r_cur  = math.sqrt( dx**2 + dy**2 + dz**2 )
        
        #assert af in af_nbr.neighbors
        #assert af_p in af_nbr_p.neighbors
        assert af_nbr_p.id == af_nbr.id
#        assert r_prev < r_cutoff and r_cur < r_cutoff, "Interaction distance is > rc = %g r_prev: %g r_cur: %g" %(r_cutoff, r_prev, r_cur)
        #assert r_prev < r_cutoff, "Interaction distance is > rc = %g r_prev: %g" %(r_cutoff, r_prev)
        flag = maxchange_criteria(r_prev, r_cur, delta_r)
        if flag == 1:
            af.ext_ids.add(af_nbr.id)
            af_nbr.ext_ids.add(af.id)
            lj_change_ext.append((af.id, af_nbr.id))
        elif flag == -1:
            lj_change_comp.append((af.id, af_nbr.id))
            af.comp_ids.add(af_nbr.id)
            af_nbr.comp_ids.add(af.id)
            
    rbond = 1.2
    
    broken = []
    for nbr_prev in nbr_list_p:
        if nbr_prev.id not in cur_ids:
            dx, dy, dz = get_dxdydz(af_p, nbr_prev, bounds_p)
            r_prev = math.sqrt( dx**2 + dy**2 + dz**2)
            if r_prev < rbond:
                broken.append((af_p.id, nbr_prev.id))
    
    formed = []
    for nbr_cur in nbr_list:
        if nbr_cur.id not in prev_ids:
            
            dx, dy, dz = get_dxdydz(af, nbr_cur, bounds)
            r_cur  = math.sqrt( dx**2 + dy**2 + dz**2 )
            
            if r_cur < rbond:
                formed.append((af, nbr_cur))
            nbr_id = nbr_cur.id
            af_nbr_p = atom_forces_p[nbr_id - 1] 
            assert af_nbr_p.id == nbr_cur.id
            af_nbr = nbr_cur
            
            dx, dy, dz = get_dxdydz(af_p, af_nbr_p, bounds_p)
            r_prev = math.sqrt( dx**2 + dy**2 + dz**2 )
            
            dx, dy, dz = get_dxdydz(af, af_nbr, bounds)
            r_cur  = math.sqrt( dx**2 + dy**2 + dz**2 )
            flag = maxchange_criteria(r_prev, r_cur, delta_r)
            if flag == 1: 
                print("Error bond can't be extended")
                lj_change_ext.append((af.id, af_nbr.id))
            elif flag == -1:
                lj_change_comp.append((af.id, af_nbr.id))
                af.comp_ids.add(af_nbr.id)
                af_nbr.comp_ids.add(af.id)
        
    return lj_change_comp, lj_change_ext, broken, formed
            

#def get_stat(af_p, nbr_list_p, af, nbr_list, bounds_p, bounds, percent):
           
def is_neighbor(a_id, neighbor_list):
    for nbr in neighbor_list:
        if a_id == nbr.id:
            return True
    return False

def get_separation(af1, af2):
    return math.sqrt((af1.x - af2.x)**2 + (af1.y - af2.y)**2 + (af1.z - af2.z)**2)


def add_neighbors(all_pair_ids, all_res, atype):
    N = len(all_pair_ids)
    N1 = len(all_res)
    assert N == N1
    contact_idx = 1000000
    for i in range(N):
        pair_count = 0
        pair_ids = all_pair_ids[i]
        atom_forces = all_res[i][atype]
        M = len(atom_forces)
        for (id1, id2) in pair_ids:
            if min(id1, id2) <= M and max(id1, id2) > M: contact_idx = min(contact_idx, i)
            if id1 > M or id2 > M: continue
            af1 = atom_forces[id1-1] ## ASSUMING ATOMS ARE SORTED ACCORDING TO THEIR IDS
            af2 = atom_forces[id2-1]
            #print(af1.id, id1, af2.id, id2)
            assert af1.id == id1 and af2.id == id2
            af1.neighbors.append(af2)
            af2.neighbors.append(af1)
            pair_count += 1
        print("i: %d number of pairs: %d" %(i, pair_count))
    return contact_idx


def construct_neighbors(all_res, all_bounds, atype, rc):
    N = len(all_res)
    pair_counts = []
    for i in range(N):
        atom_forces = all_res[i][atype]
        bounds      = all_bounds[i]
        pair_count  = create_neighbor_lists(atom_forces, bounds, rc)
        pair_counts.append(pair_count)
    return pair_counts
        
def visualize_neighbors(atom_forces, bounds):
    '''Assumption: atom_forces are sorted by id'''
    fig = plt.figure()
    ax = Axes3D(fig)
    for af in atom_forces:
        for nbr in af.neighbors:
            dx, dy, dz = get_dxdydz(nbr, af, bounds)
            x_next = af.x + dx
            y_next = af.y + dy
            z_next = af.z + dz
            if af.id - nbr.id == -1 or af.id - nbr.id == 1:
                ax.plot([af.x, x_next], [af.y, y_next], [af.z, z_next], color = 'black')
            else:
                ax.plot([af.x, x_next], [af.y, y_next], [af.z, z_next], color =  'cyan')
            #ax.scatter(nbr.x, nbr.y, nbr.z, c = 'green')
        ax.scatter(af.x, af.y, af.z, c = 'black')
    plt.show()
    
def plot_neighbor_changes(times,  changes_comp, changes_ext, breaks, formations, pair_counts, vz, dt, d0):
    ccfrac = [changes_comp[i]/pair_counts[i]    for i in range(len(changes_comp))]
    cefrac = [changes_ext[i]/pair_counts[i]    for i in range(len(changes_ext))]
    bfrac = [breaks[i]/pair_counts[i]     for i in range(len(breaks))]
    ffrac = [formations[i]/pair_counts[i] for i in range(len(formations))]
    ds = [vz*dt*t - d0 for t in times]
    plt.plot(ds[:-1], ccfrac, label = "Compressed")
    plt.plot(ds[:-1], cefrac, label = "Extended")
    plt.plot(ds[:-1], bfrac, label = "Broke")
    plt.plot(ds[:-1], ffrac, label = "Formed")
    plt.legend()
    plt.show()

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
    pair_counts = []
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
                pair_counts.append(count)
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
            elif 'ITEM: BOX BOUNDS' in line:
                r_boundary = True
                d = 3
                print("Next 3 lines are bondaries")
            elif 'ITEM: ENTRIES' in line:
                r_entries = True
                print("Pair id lines are coming")
    return all_pair_ids, pair_counts

def get_interactions(filename, t_start, t_end, types, interacting = False):
    assert t_start <= t_end
    r_time = False
    r_atom_count = False
    r_boundary = False
    r_atoms = False
    skip = False
    dim = 0
    t = 0
    max_r = 0
    pbcX, pbcY, pbcZ = False, False, False #boundary periodicity flags
    all_res = []
    all_bounds = []
    times = []
    cur_b = []
    res = dict()
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                times.append(time)
                #print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                res = dict()
                #print("# of atoms: %d" %count)
                r_atom_count = False
            elif r_boundary:
                if dim > 0:
                    words = line.strip().split()
                    lo, hi = map(float, words)
                    if   dim == 3: bounds.set_x_bounds(lo, hi)
                    elif dim == 2: bounds.set_y_bounds(lo, hi)
                    elif dim == 1: bounds.set_z_bounds(lo, hi)
                dim -= 1
                if dim == 0:
                    all_bounds.append(bounds)
                    r_boundary = False
            elif r_atoms:
                if 'ITEM: TIMESTEP' in line:
                    r_atoms = False
                    r_time = True
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
                    words = line.split(' ', 9)
                    if len(words) > 9:
                        a_id, mol, atype, x, y, z, fx, fy, fz, arributes  = words
                    else:
                        a_id, mol, atype, x, y, z, fx, fy, fz  = words
                    a_id, mol, atype, x, y, z, fx, fy, fz = int(a_id), int(mol), int(atype), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    if atype not in types: continue
                    if interacting and abs(fx) < epsilon and abs(fy) < epsilon and abs(fz) < epsilon: continue
                    #choose interacting atoms
                    radius = math.sqrt(x**2 + y**2)
                    if atype not in res: res[atype] = []
                    res[atype].append(AtomicForces(a_id, mol, atype, x, y, z, fx, fy, fz))
                    max_r = max(max_r, radius)
                    
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if len(res) != 0:
                    for type, atom_forces in res.items():
                        l = sorted(atom_forces, key = lambda af: af.id)
                        res[type] = l
                    all_res.append(res)
                elif len(times) > 0:
                    times.pop()
                r_time = True
                t += 1
                if t > t_end:
                    break
                if t < t_start:
                    skip = True
                else:
                    skip = False
                #print("Next is timestep")
            elif 'ITEM: NUMBER OF ATOMS' in line:
                r_atom_count = True
                #print("Next is number of atoms")
            elif 'ITEM: BOX BOUNDS' in line:
                words = line.strip().split()
                pbcX = True if words[-3] == "pp" else False
                pbcY = True if words[-2] == "pp" else False
                pbcZ = True if words[-1] == "pp" else False
                r_boundary = True
                bounds = Boundary(line)
                dim = 3
                #print("Next 3 lines are bondaries")
            elif 'ITEM: ATOMS' in line:
                r_atoms = True
                atr_line = line[line.index('S')+1:]
                atr_words = atr_line.split(' ')
                #print("Atom coordinate lines are coming")
    return all_res, all_bounds, times



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

def plot_stresszz_d(all_res, times, v, type):
    areas = []
    #times = []
    ds  = []
    strs_z = []
    fzs = []
    num_intrs = []
    dt = 0.01
    d_start = 0
    t0 = 0
    d0 = 0
    del_z0 = 0
    first_contact = True
    ts = 0
    avg_strs, cnt = 0, 0
    fz_prev = 0
    E_modulus = 0
    print(len(all_res), len(times))
    for i in range(len(all_res)):
        t = times[i]
        print("New timestep: %d" %t)
        res = all_res[i]
        interacting_af, points = interacting_particles(res[type])
        count = len(interacting_af)
        if count == 0:
            continue
        elif count < 3:
            #continue
            if first_contact:
                t0 = times[i]
                d0 = t0 * v * dt
                print("First contact happens at t0 : %d d0: %g" %(t0, d0))
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
        
        del_z = t * v * dt - d0
        if del_z < del_z0:
            continue
        else:
            del_z -= del_z0
        str_z = fz / area
        #if str_z > 3:
        #    continue
        areas.append(area)
        fzs.append(fz)
        num_intrs.append(count)
        #times.append(t)
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


def plot_changes(ds, ccfrac, cefrac, contactd, delta_r):
    cfrac = [ccfrac[i] + cefrac[i] for i in range(len(ccfrac))]
    fig = plt.figure()
    plt.title("Changes by >= %g vs displacement of the tip" %(delta_r))
    plt.plot(ds, ccfrac, label = "Compressed")
    plt.plot(ds, cefrac, label = "Extended")
    plt.plot(ds, cfrac,  label = "Changed" )
    plt.axvline(x=contactd, color = 'red', label = "Contact point: %g" %contactd)
    plt.xlabel("d")
    plt.ylabel("Fraction")
    plt.legend()
    #plt.show()

def plot_breaks(ds, bfrac, ffrac, contactd):
    plt.title("Broken and formed bond fractions  vs displacement of the tip")
    plt.plot(ds, bfrac, label = "Broke")
    plt.plot(ds, ffrac, label = "Formed")
    plt.xlabel("d")
    plt.ylabel("Fraction")
    plt.axvline(x=contactd, color = 'red', label = "Contact point: %g" %contactd)
    plt.legend()
    #plt.show()


def append_vis_file(all_res, all_bounds, times, filename, step):
    N = len(times)
    with open(filename, 'a') as f:
        for i in range(N-step):
            bounds = all_bounds[i]
            res = all_res[i]
            #print time
            time = times[i]
            f.write("ITEM: TIMESTEP\n%d\n" %time)
            #print atom count
            atom_count = 0
            for key, val in res.items():
                atom_count += len(val)
            f.write("ITEM: NUMBER OF ATOMS\n%d\n" %atom_count)
            #print bounds
            f.write(bounds.get_bounds_text())
            f.write("ITEM: ATOMS id mol type x y z fx fy fz ext comp chg brk frm\n")
            for key, atom_forces in res.items():
                for af in atom_forces:
                    f.write("%d %d %d %f %f %f %f %f %f " %(af.id, af.mol, af.type, af.x, af.y, af.z, af.fx, af.fy, af.fz))
                    f.write("%f %f %f %f %f\n" %(af.atr['ext'], af.atr['comp'],af.atr['ext'] + af.atr['comp'], af.atr['broken'], af.atr['formed']))
            

def save_bond_change_stats(data, filename):
    rows, cols = data.shape
    with open(filename, 'w') as f:
        f.write("time compressions extentions breaks formation c-e e-c pair-counts\n")
        for i in range(rows):
            text = ""
            for j in range(cols):
                val = int(data[i, j])
                text += "%d " %val
            text = text.strip()
            f.write(text)
            f.write("\n")
def visualize_lj_bond_stats(css):
    #M, N = 2000, 256
    #T = 0.0001
    #r = 10
    dt = css.dt
    #cang = 45
    tip_type = 2
    glass = 1
    types = [glass]
    atype = glass
    #rc = 1.5
    vz = css.vz
    d0 = 0 #2.2 
    t_init, t_final = css.t_init, css.t_final
    t_step = css.t_step
    step = 1 #step for calculations
    ccfrac, cefrac, bfrac, ffrac  = [], [], [], []
    comp_ext_frac, ext_comp_frac = [], []
    data_count = t_final-t_init
    print("data count", data_count)
    data = np.zeros( (data_count, 8) )
    ds = []
    delta_r = 0.3
    contactd = 100000000000
    filename = vis_data_path + 'visualize_M%d_N%d_T%g_r%d_cang%d.out' %(css.M, css.N, css.T, css.r, css.cang)
    filenameinteractions = vis_data_path + 'pairids_M%d_N%d_T%g_r%d_cang%d.out' %(css.M, css.N, css.T, css.r, css.cang)
    filename_heat = out_data_path + "conetip_M%d_N%d_T%g_sphR%d_cang%d_nve.txt" %(css.M, css.N, css.T, css.r, css.cang)
    changes_filename = "changes_M%d_N%d_T%g_r%d_cang%d_p%g_stp%d.png" %(css.M, css.N, css.T, css.r, css.cang, delta_r, step)
    breaks_filename = "breaks_M%d_N%d_T%g_r%d_cang%d_p%g_stp%d.png" %(css.M, css.N, css.T, css.r, css.cang,  delta_r, step)
    vischanges_filename = "visualizechanges_M%d_N%d_T%g_r%d_cang%d_p%g.out" %(css.M, css.N, css.T, css.r, css.cang, delta_r)
    #remove files if they already exist
    remove_file(changes_filename)
    remove_file(breaks_filename)
    remove_file(vischanges_filename)
    rc = 1.5 + 0.001
    idx = 0 #CHANGE THIS
    contactd = 3
    j = 0
    for t_start in range(t_init, t_final, t_step):
        t_end = t_start + t_step + step - 1
        all_res, all_bounds, times = get_interactions(filename, t_start, t_end, types, interacting = False)
        pair_counts = construct_neighbors(all_res, all_bounds, atype, rc)
        #all_inter, pair_counts = get_pair_interactions(filenameinteractions, t_start, t_end)
        #for i in range(len(all_res)):
        #    test_neighbor_lists(all_res[i][glass], all_inter[i], glass, all_bounds[i], rc)
        #
        #idx = add_neighbors(all_inter, all_res, glass)
        #test_neighbor_lists(all_res[0][glass], all_inter[0], glass, all_bounds[0], rc)
        #if css.T > 0.01:
        #    reconstuct_ave_lj_bonds(all_res, atype, all_bounds, percent)
        #if idx < t_final:
        #    contactd = min(contactd, times[idx]*vz*dt)
        changes_comp, changes_ext, breaks, formations, comp_then_ext, ext_then_comp = get_lj_bond_stats(all_res, glass, all_bounds, delta_r, step)
        for i in range(len(changes_comp)):
            if j >= data_count: break
            data[j, 0] = times[i]
            data[j, 1] = changes_comp[i]
            data[j, 2] = changes_ext[i]
            data[j, 3] = breaks[i]
            data[j, 4] = formations[i]
            data[j, 5] = comp_then_ext[i]
            data[j, 6] = ext_then_comp[i]
            data[j, 7] = pair_counts[i]
            j          += 1

        ccfrac.extend( [changes_comp[i]/pair_counts[i]    for i in range(len(changes_comp))] )
        cefrac.extend( [changes_ext[i]/pair_counts[i]     for i in range(len(changes_ext))] )
        bfrac.extend(  [breaks[i]/pair_counts[i]          for i in range(len(breaks))] )
        ffrac.extend(  [formations[i]/pair_counts[i]      for i in range(len(formations))] )
        ds.extend(     [vz*dt*times[i] - d0               for i in range(len(times)-step)] )
        comp_ext_frac.extend([comp_then_ext[i]/pair_counts[i]    for i in range(len(comp_then_ext))])
        ext_comp_frac.extend([ext_then_comp[i]/pair_counts[i]    for i in range(len(ext_then_comp))])
        append_vis_file(all_res, all_bounds, times, vischanges_filename, step)
        print("t start: %d\n" %t_start, flush=True)
        gc.collect() #
    print(len(ds), len(ffrac))
    save_bond_change_stats(data, "stats_M%d_N%d_T%g_r%d_cang%d_p%g.txt" %(css.M, css.N, css.T, css.r, css.cang, delta_r))
    plot_changes(ds, ccfrac, cefrac, contactd, delta_r)
    plt.plot(ds, comp_ext_frac, label = "CE")
    plt.plot(ds, ext_comp_frac, label = "EC")
    plt.legend()
    #Write rate of rearragements to a file
    deld = ds[step] - ds[0]
    rate = [(ccfrac[i] + cefrac[i])/deld for i in range(len(ccfrac))]
    with open("rate_T%g_deld%g.txt" %(css.T, deld), 'w') as f:
        for item in rate:
            f.write("%s\n" %item)
    
    if css.T < 0.01:
        scale =  0.8214 * 0.00000416
        heat_stats(filename_heat, scale, step)
    plt.savefig(changes_filename)
    plt.close()
    plot_breaks(ds, bfrac, ffrac, contactd)
    plt.savefig(breaks_filename)
    plt.close()
    


def visualize_particles(css):
    dt = css.dt
    #cang = 45
    tip_type = 2
    glass = 1
    types = [glass]
    atype = glass
    #rc = 1.5
    vz = css.vz
    d0 = 0 #2.2
    t_init, t_final = css.t_init, css.t_final
    t_step = css.t_step
    filename = vis_data_path + 'particles_M%d_N%d_T%g_r%d_cang%d.out' %(css.M, css.N, css.T, css.r, css.cang)
    ts = []
    pos_dict = {}
    labels = ['x', 'y', 'z', 'r' ]
    for t_start in range(t_init, t_final, t_step):
        t_end = t_start + t_step
        all_res, all_bounds, times = get_interactions(filename, t_start, t_end, types, interacting = False)
        ts.extend(times[:-1])
        for i in range(t_step):            
            atom_forces = all_res[i][atype]
            if not pos_dict:
                for af in atom_forces:
                    pos_dict[af.id] = [ [], [], [], [] ]
            for af in atom_forces:
                pos_dict[af.id][0].append(af.x)
                pos_dict[af.id][1].append(af.y)
                pos_dict[af.id][2].append(af.z)
                pos_dict[af.id][3].append(math.sqrt(af.x**2 + af.y**2 + af.z**2))
                
    for key, val in pos_dict.items():
        for i in range(len(labels)):
            fig = plt.figure()
            plt.title("Positions of atom %d vs timestep" %(key))
            plt.plot(ts, val[i], label = labels[i])
            plt.xlabel("t")
            plt.ylabel(labels[i])
            plt.legend()
            plt.savefig("positions_%d_%s.png" %(key, labels[i]))
            plt.close()
        ##plt.plot(ts, val[1], label = "y")
        #plt.plot(ts, val[2], label = "z" )
        #plt.plot(ts, val[3], label = "r" )
       
    
                
            
def remove_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    
def main():
    #plot_nforce_vs_cont_area()
    #substrate_type = 1
    #tip_type = 2
    #oligomer_type = 3
    #visualize_particles(css)
    visualize_lj_bond_stats(css)
    return
    M, N = 2000, 256
    r = 10
    dt = 0.01
    cang = 60
    tip_type = 2
    glass = 1
    t_start = 1
    t_end = 100
    types = [tip_type]
    rc= 1.5
    vz = 0.0001
    d0 = 0 #2.2
    #bond testing
    #filename = '../visfiles/viscomp_M%d_N%d.out' %(M, N)
    #frames = get_frames(filename, t_start, t_end)
    #broken_bonds(frames, glass, M, N, 1.5)
    #return
    filename = '../visfiles/visualize_M%d_N%d_r%d_cang%d.out' %(M, N, r, cang)
    filenameinteractions = '../visfiles/pairids_M%d_N%d_r%d_cang%d.out' %(M, N, r, cang)
    #append_bondlens(filename, types, M, N)
    #return
    all_res, bounds, times = get_interactions(filename, t_start, t_end, types, interacting = True)
    #print(times)
    #return
    #all_inter, pair_counts = get_pair_interactions(filenameinteractions, t_start, t_end)
    #add_neighbors(all_inter, all_res, glass)
    #changes_comp, changes_ext, breaks, formations = get_lj_bond_stats(all_res, glass, bounds, 0.2)
    #plot_neighbor_changes(times, changes_comp, changes_ext, breaks, formations, pair_counts, vz, dt, d0)
    #af_2_visualize = all_res[2][glass][0:20]
    #visualize_neighbors(af_2_visualize, bounds)
    #return
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
    plot_stresszz_d(all_res, times, vz, tip_type)
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
