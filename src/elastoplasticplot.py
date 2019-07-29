
import matplotlib.pyplot as plt
import math

theta  = 20
theta_rad = math.radians(theta)
r = 10

d_boundary = r * ( 1 - math.cos(theta_rad) )
a_boundary = r * math.sin(theta_rad)

print(a_boundary, d_boundary)

filename = "../push_displ_r%d.dat" %r
z0 = 0 
z = []
fzs = []
sigmaz = []
areas = []
i = 0
with open(filename, 'r') as file:
    for line in file:
        i+=1
        if i < 3: continue
        words = line.split()
        fz = float(words[5])
        if i == 3:
            z0 = float(words[2])
        delta_z = abs(float(words[2])-z0)
        r_eff = 0
        if delta_z <= 0 or delta_z > 15: continue
        
        if delta_z <= d_boundary:
            r_eff = r * math.sin(math.acos(1 - delta_z/r))
            area = math.pi * (r_eff**2)
        else:
            r_eff =  a_boundary + (delta_z - d_boundary) / math.tan(theta_rad)
            area = math.pi * (r_eff**2)
        #if fz / (delta_z**(2)) > 100: continue
        
        z.append(delta_z)
        fzs.append(fz)
        areas.append(area)
        sigmaz.append(fz/area)

plt.title("Tip R = %d, Cone angle (w.r.t. radial direction) = %d" %(r,theta))
plt.ylabel('$\sigma$', fontsize = 20)
plt.xlabel(r'$d$', fontsize = 20) 
plt.plot(z, sigmaz)
plt.legend()
plt.show()
