from fcc import *

stepsize = 0.004
i = 0
sigma = 1.0
d = 2.0**(1.0/6.0)
normal_load = 1000
tip_radius = 1000

def get_force_vs_dist(filename):
    forces_x =  []
    forces_y =  [] 
    dists = []
    i = 0
    max_fx = -1
    max_fy = -1
    with open(filename) as file:
        for line in file:
            if i == 0:
                i += 1
                continue
            t,_,_,_,fx,fy,_ = line.split(' ')
            fx, fy = float(fx), float(fy)
            max_fx = min(max_fx, fx) #taking minimum because force is in negative direction
            max_fy = min(max_fy, fy)
            forces_x.append(fx)
            forces_y.append(fy)
            dists.append( (stepsize * sigma) * float(t) / d)
    return forces_x, forces_y, dists

def plot_force_vs_dist(load_radius_pairs):
    fig, ax = plt.subplots()
    for normal_load, tip_radius in load_radius_pairs:
        filename = '../friction_data_f%d_r%d.dat' %(normal_load, tip_radius)
        forces_x, forces_y, dists = get_force_vs_dist(filename)
        norm_forces_x = [fx/normal_load for fx in forces_x]
        norm_forces_y = [fy/normal_load for fy in forces_y]
        ax.plot(dists, norm_forces_x, label = 'Load: %d' %normal_load)
    ax.legend()
    plt.ylabel('Fx/N', fontsize = 20)
    plt.xlabel(r'$x/d$', fontsize = 20)
    plt.title("Friction force vs displacement.")
    plt.show()


def main():
    load_radius_pairs = [(601, 1000),(600, 1000)]#, (500, 1000), (350, 1000)]#, (300, 1000), (200, 1000), (100, 1000)]
    plot_force_vs_dist(load_radius_pairs)

if __name__ == "__main__":
    main()


'''
norm_forcesx = [fx/max_fx for fx in forcesx]
plt.plot(timestep, norm_forcesx)
plt.suptitle('Friction force vs displacement. Load: %d, Radius: %d'  %(normal_load, tip_radius), fontsize = 20)
plt.axhline(0, color = 'black')
plt.show()
'''
