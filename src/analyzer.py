from fcc import *
import secrets
from scipy import stats
from scipy.spatial import ConvexHull
from scipy.spatial import KDTree
import os
import argparse
from thermoparser import heat_stats
from boxcell import *
from dumpparser import *
import gc
epsilon = 0.000001
vis_data_path = '../visfiles/'
out_data_path = '../outputfiles/'
#idx_id, idx_mol, idx_type, idx_x, idx_y, idx_z, idx_fx, idx_fy, idx_fz, _ = line.split(' ')
#vis_data_path = '/home/kkurman1/equilibrate/identtip/visfiles/'
#out_data_path =  '/home/kkurman1/equilibrate/identtip/outputfiles/'
#vis_data_path = '/home/kkurman1/GlassExperiments/stiff/indentation/visfiles/'
#out_data_path =  '/home/kkurman1/GlassExperiments/stiff/indentation//outputfiles/'


class SimulationSettings:
    def __init__(self, force, poisson, G_shear_mod, R, sigma, d):
        self.force = force
        self.poisson = poisson
        self.G_shear_mod = G_shear_mod
        self.R = R
        self.d = d


class FileNames:
    def __init__(self, M, N, T, r, cang, stiff, is_cone):
        if is_cone:
            self.make_conetip_names(M, N, T, r, cang, stiff)
        else:
            self.make_spheretip_names(M, N, T, r, stiff)
    def make_conetip_names(self, M, N, T, r, cang, stiff):
        name  = 'visualize_stiff_M%d_N%d_T%g_r%d_cang%d.out'    %(M, N, T, r, cang) if stiff else 'visualize_M%d_N%d_T%g_r%d_cang%d.out'   %(M, N, T, r, cang)
        self.vis   = vis_data_path + name
        name  = 'pairids_stiff_M%d_N%d_T%g_r%d_cang%d.out'      %(M, N, T, r, cang) if stiff else 'pairids_M%d_N%d_T%g_r%d_cang%d.out'     %(M, N, T, r, cang)
        self.inter = vis_data_path + name
        name  = 'conetip_stiff_M%d_N%d_T%g_r%d_cang%d_nve.txt'  %(M, N, T, r, cang) if stiff else 'conetip_M%d_N%d_T%g_r%d_cang%d_nve.txt' %(M, N, T, r, cang)
        self.heat  = out_data_path + name
        name  = 'stats_stiff_M%d_N%d_T%g_r%d_cang%d.out'        %(M, N, T, r, cang) if stiff else 'stats_M%d_N%d_T%g_r%d_cang%d.out'       %(M, N, T, r, cang)
        self.stats  = out_data_path + name
        #changes = "changes_M%d_N%d_T%g_r%d_cang%d_p%g_stp%d.png" %(css.M, css.N, css.T, css.r, css.cang, delta_r, step)
        #breaks = "breaks_M%d_N%d_T%g_r%d_cang%d_p%g_stp%d.png" %(css.M, css.N, css.T, css.r, css.cang,  delta_r, step)
        #vis_changes = "visualizechanges_M%d_N%d_T%g_r%d_cang%d_p%g.out" %(css.M, css.N, css.T, css.r, css.cang, delta_r)
    def make_spheretip_names(self, M, N, T, r, stiff):
        name  = 'vis_sphere_stiff_M%d_N%d_T%g_r%d.out'    %(M, N, T, r) if stiff else 'vis_sphere_M%d_N%d_T%g_r%d.out'   %(M, N, T, r)
        self.vis   = vis_data_path + name
        name  = 'pairids_stiff_M%d_N%d_T%g_r%d.out'       %(M, N, T, r) if stiff else 'pairids_M%d_N%d_T%g_r%d.out'      %(M, N, T, r)
        self.inter = vis_data_path + name
        name  = 'conetip_stiff_M%d_N%d_T%g_r%d_nve.txt'   %(M, N, T, r) if stiff else 'conetip_M%d_N%d_T%g_r%d_nve.txt'  %(M, N, T, r)
        self.heat  = out_data_path + name
        name  = 'stats_stiff_M%d_N%d_T%g_r%d.out'         %(M, N, T, r) if stiff else 'stats_M%d_N%d_T%g_r%d.out'        %(M, N, T, r)
        self.stats  = out_data_path + name
        #changes = "changes_M%d_N%d_T%g_r%d_cang%d_p%g_stp%d.png" %(css.M, css.N, css.T, css.r, css.cang, delta_r, step)
        #breaks = "breaks_M%d_N%d_T%g_r%d_cang%d_p%g_stp%d.png" %(css.M, css.N, css.T, css.r, css.cang,  delta_r, step)
        #vis_changes = "visualizechanges_M%d_N%d_T%g_r%d_cang%d_p%g.out" %(css.M, css.N, css.T, css.r, css.cang, delta_r)

class ConesimSettings:
    def __init__(self, M, N, T, r, cang, vz, dt):
        self.M    = M       #number of chains
        self.N    = N       #number of monomers per chain
        self.T    = T       #Temperature of the system
        self.r    = r       #radius of the spherical tip of the cone
        self.cang = cang    #angle that cone makes with substrate plane
        self.vz   = vz      #velocity of the cone
        self.dt   = dt      #timestep of the simulation
        print("Simulation settings used for analysis")
        print("M : %d, N: %d, T: %g, r: %g, cang: %g, vz: %g, dt: %g" %(M, N, T, r, cang, vz, dt), flush = True)
    def set_analysisvals(self, t_init, t_final, t_step):
        self.t_init  = t_init  #analysis starts from this dump 
        self.t_final = t_final #analysis continue until this dump
        self.t_step  = t_step  #us this much dumps for analysis at single iteration
        print("Analysis starts from dump time: %g continues until: %g in steps of %g" %(t_init, t_final, t_step))

'''      
Temp = 0.0001
is_stiff = True
css = ConesimSettings(2000, 256, Temp, 0, 0, 0.0001, 0.01)
css.set_analysisvals(1, 30, 1)
filenames = FileNames(2000, 256, Temp, 0, 0, is_stiff)
'''
filenames, css = None, None

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
            d = 1.12
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
    
    plt.show()
    fig, ax = plt.subplots(1, 1)
    plt.plot(ds, strs_z, 'b')
    plt.ylabel("$\sigma_{zz}$", rotation = 0)
    plt.xlabel("$d$")
    plt.legend()
    plt.plot()
    plt.show()
    print("Avg stress: %g" %(avg_strs/cnt))
    
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
    filename             = filenames.vis
    filenameinteractions = filenames.inter
    filename_heat        = filenames.heat
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
    save_bond_change_stats(data, filenames.stats)
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
    
    #if css.T < 0.3:
    #    scale =  0.8214 * 0.00000416
    #    heat_stats(filename_heat, scale, step)
    plt.savefig(changes_filename)
    plt.close()
    plot_breaks(ds, bfrac, ffrac, contactd)
    plt.savefig(breaks_filename)
    plt.close()


def get_average_stretches(atomic_forces, M, N, bounds):
    '''Assumption atom_forces are sorted by id 
    and M chains are of length N first polymer (chain) is at idx location 0:N-1'''
    rc = 1.5
    rcsq = rc**2
    broken_count  = 0
    Rtot = 0.0
    Rs = np.zeros([M, 1])
    for i in range(0, M*N, N):
        Rx, Ry, Rz = 0, 0, 0
        for j in range(i, i+N-1):
            af_cur  = atomic_forces[j]
            af_next = atomic_forces[j+1]
            dx, dy, dz = get_dxdydz(af_next, af_cur, bounds)
            Rx += dx
            Ry += dy
            Rz += dz
            rsq = dx*dx + dy*dy + dz*dz
            if rcsq < rsq:
                broken_count += 1
        idx = int(i/N)
        Rs[idx]    = ( Rx*Rx + Ry*Ry + Rz*Rz )**(1/2)
    
    R_ave = np.mean(Rs)
    #print(R_ave)
    R_std_err = np.std(Rs) / len(Rs)**(1/2)
    return R_ave, R_std_err, broken_count

def visualize_stretches(css):
    print("Visualizing stretches")
    glass = 1
    atype = glass
    types = [atype]
    t_init, t_final = css.t_init, css.t_final
    t_step = css.t_step
    vz = css.vz
    d0 = 0
    R_avgs, R_stderrs = [], []
    chain_breaks = []
    ds = []
    filename = filenames.vis
    for t_start in range(t_init, t_final, t_step):
        t_end = t_start + t_step - 1
        all_res, all_bounds, times = get_interactions(filename, t_start, t_end, types, interacting = False)
        for i in range(len(all_res)):
            res                      = all_res[i]
            bounds                   = all_bounds[i]
            time                     = times[i]
            atomic_forces            = res[atype]
            Rave, R_std_err, broken  = get_average_stretches(atomic_forces, css.M, css.N, bounds)
            R_avgs.append(Rave)
            R_stderrs.append(R_std_err)
            chain_breaks.append(broken)
        ds.extend([vz * css.dt * times[i] - d0   for i in range(len(times))])
    with open("stretches_stiff_M%d_N%d_T%g_r%d_cang%d.txt" %(css.M, css.N, css.T, css.r, css.cang),'w') as f:
        for i in range(len(ds)):
            f.write("%g %g %g %g\n" %(ds[i], R_avgs[i], R_stderrs[i], chain_breaks[i]))
    plt.plot(ds, R_avgs)
    plt.show()

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
       

def vis_hardness(css):
    print("Visualizing hardness")
    dt = css.dt
    #cang = 45
    tip_type = 2
    glass = 1
    types = [tip_type]
    atype = glass
    #rc = 1.5
    vz = css.vz
    d0 = 0 #2.2 
    t_init, t_final = css.t_init, css.t_final
    t_step = css.t_step
    filename = filenames.vis
    all_res, bounds, times = get_interactions(filename, t_init, t_final, types, interacting = True)
    plot_stresszz_d(all_res, times, vz, tip_type)
             
            
def remove_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

    
def main():
    parser = argparse.ArgumentParser(description = "Contact analysis")
    parser.add_argument('--M',    type=int,   default = 2000,   help = '# of chains in a melt')
    parser.add_argument('--N',    type=int,   default = 256,    help = '# of monomers per chain')
    parser.add_argument('--T',    type=float, default = 0.0001, help = 'Temperature of the system')
    parser.add_argument('--vz',   type=float, default = 0.0001, help = 'Indenting tip velocity')
    parser.add_argument('--dt',   type=float, default = 0.01,   help = 'Integration time step of the simluation')
    parser.add_argument('--r',    type=float, default = 10,     help = 'Radius of the spherical tip')
    parser.add_argument('--cang', type=float, default = 45,     help = 'Angle the cone surface makes with the horizontal plane')
    parser.add_argument('--stiff',     action = 'store_true', default = False, help = 'True if polymer stiff (i.e. there is an angle style defined)')
    parser.add_argument('--hardness',  action = 'store_true', default = False, help = 'True if analysis of hardness is needed')
    parser.add_argument('--bondstats', action = 'store_true', default = False, help = 'True if analysis of hardness is LJ bond lengths is needed')
    parser.add_argument('--stretches', action = 'store_true', default = False, help = 'True if analysis of end-end polymer lengths is LJ bond lengths is needed')
    parser.add_argument('--conetip',   action = 'store_true', default = False, help = 'True if tip is of spherical shape')
    args = parser.parse_args()
    Temp = 0.0001
    is_stiff = True
    global css, filenames
    print(args.M, args.N, args.T, args.r, args.cang, args.stiff, args.conetip)
    css = ConesimSettings(args.M, args.N, args.T, args.r, args.cang, args.vz, args.dt)
    css.set_analysisvals(20, 22, 1)
    filenames = FileNames(args.M, args.N, args.T, args.r, args.cang, args.stiff, args.conetip)

    print(args.stiff, args.M)
    #plot_nforce_vs_cont_area()
    #substrate_type = 1
    #tip_type = 2
    #oligomer_type = 3
    #sargs.hardness = True
    if args.hardness:
        vis_hardness(css)
    if args.bondstats:
        visualize_lj_bond_stats(css)
    if args.stretches:
        visualize_stretches(css)
    #vis_hardness(css)
    #visualize_particles(css)
    #visualize_stretches(css)
    #vis_hardness(css)
    #visualize_lj_bond_stats(css)
    return

if __name__=="__main__":
    main()

#plt.plot(timestep, forcesx)
#plt.suptitle('Friction force vs timestep', fontsize = 20)
#plt.ylabel('Fx', fontsize = 20)
#plt.xlabel('t', fontsize = 20)
#plt.axhline(0, color = 'black')
#plt.show()
