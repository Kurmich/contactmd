from fcc import *
import secrets
from scipy.spatial import ConvexHull
forces = [500.0]
poisson = 0.5
G_shear_mod = 16.0
E_star = 1#4 * G_shear_mod
R = 1000.0
sigma = 1.0 #2**(1/6)
d = 2.0**(1.0/6.0) * sigma
atom_N = [1, 2, 3, 4]


#idx_id, idx_mol, idx_type, idx_x, idx_y, idx_z, idx_fx, idx_fy, idx_fz, _ = line.split(' ')


class SimulationSettings:
    def __init__(self, force, poisson, G_shear_mod, R, sigma, d):
        self.force = force
        self.poisson = poisson
        self.G_shear_mod = G_shear_mod
        self.R = R
        self.d = d
    
class AtomicForces:
    def __init__(self, type, x, y, z, fx, fy, fz):
        self.x, self.y, self.z = x, y, z
        self.fx, self.fy, self.fz = fx, fy, fz
        self.radius = self.get_radius(x, y, 0)

    def get_radius(self, x, y, z):
        return math.sqrt(x**2 + y**2 + z**2)
    def __lt__(self, other):
        return self.radius < other.radius
    def __eq__(self, other):
        return self.radius == other.radius


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


def get_interactions(filename, t_start, t_end):
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
    cur_frame = None
    frames = []
    with open(filename) as file:
        for line in file:
            #check what kind of data to expect in this line
            if r_time:
                time = int(line)
                print("Time step: %d" %time)
                r_time = False
            elif r_atom_count:
                count = int(line)
                cur_frame = np.zeros([count, 8])
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
                    id, mol, type, x, y, z, fx, fy, fz, _  = line.split(' ')
                    id, mol, type, x, y, z, fx, fy, fz = int(id), int(mol), int(type), float(x), float(y), float(z), float(fx), float(fy), float(fz)
                    cur_frame[id-1, 0], cur_frame[id-1, 1] = mol, type
                    cur_frame[id-1, 2], cur_frame[id-1, 3], cur_frame[id-1, 4] = x, y, z
                    cur_frame[id-1, 5], cur_frame[id-1, 6], cur_frame[id-1, 7] = fx, fy, fz
                    #choose interacting atoms
                    radius = math.sqrt(x**2 + y**2)
                    if type not in res: res[type] = []
                    res[type].append(AtomicForces(type, x, y, z, fx, fy, fz))
                    max_r = max(max_r, radius)
                    
                    

            #set what kind of data to expect in next lines
            if 'ITEM: TIMESTEP' in line:
                if cur_frame is not None:
                    frames.append(cur_frame)
                if len(res) != 0:
                    for type, atom_forces in res.items():
                        l = sorted(atom_forces)
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
                print("Atom coordinate lines are coming")
    return all_res, frames


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
    epsilon = 0.000001
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


def plot_stresszz_d(all_res, type):
    areas = []
    times = []
    ds  = []
    strs_z = []
    dt = 0.005 * 500000
    v = 0.0001
    d_start = 10 * dt * v
    for t in range(len(all_res)):
        print("New timestep: %d" %t)
        res = all_res[t]
        interacting_af, points = interacting_particles(res[type])
        #_, _, all_fz = get_total_forces(res[type])
        fx, fy, fz = get_total_forces(interacting_af)
        hull = ConvexHull(points)
        #plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
        #plt.show()
        print("scipy area: %g" %hull.area)
        area = shoelace_area(points[hull.vertices,0], points[hull.vertices,1])
        #hull, max_r, area = jarvis(interacting_af, visualize=True)
#        area = math.pi * max_r**2
        areas.append(area)
        times.append(t)
        ds.append(d_start + t * v * dt)
        strs_z.append(fz/area)
        print(fz, area)
    plt.plot(ds, strs_z)
    plt.show()

def main():
    #plot_nforce_vs_cont_area()
    #substrate_type = 1
    #tip_type = 2
    #oligomer_type = 3
    tip_type = 3
    glass = 2
    t_start = 17
    t_end = 70
    #filename = '../visualize_f500_r1000.out'
    filename = '../visualize_M500_N500_r5.out'
    all_res, frames = get_interactions(filename, t_start, t_end)
    #plot_layer_density(glass, frames, 20)
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

if __name__=="__main__":
    main()

#plt.plot(timestep, forcesx)
#plt.suptitle('Friction force vs timestep', fontsize = 20)
#plt.ylabel('Fx', fontsize = 20)
#plt.xlabel('t', fontsize = 20)
#plt.axhline(0, color = 'black')
#plt.show()
