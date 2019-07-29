'''
boundary p p p
units lj
variable angle equal 60*PI/180
variable r   equal 10
variable dz  equal $r-$r*cos(${angle})
variable rlo equal $r*sin(${angle})
variable h   equal 30
variable rhi equal ${rlo}+$h*tan(${angle})
variable maxh equal ${dz}+$h
atom_style atomic
region mysphere sphere 0 0 $r $r
region mycone cone z 0 0 ${rlo} ${rhi} ${dz} ${maxh}
region spherecone union 2 mycone mysphere
create_box 1 spherecone
mass 1 1
lattice fcc 7
create_atoms 1 region spherecone
run 0
write_data cone nocoeff
'''

import matplotlib.pyplot as plt
import math
import numpy as np

angle = 70
angle_rad = angle * math.pi / 180
r = 1
dz = r - r * math.cos(angle_rad)
rlo = r * math.sin(angle_rad)
h = 25
rhi = rlo + h /  math.tan( angle_rad )
maxh = dz + h
circle = plt.Circle((0,r),r, fill = False)

fig, ax = plt.subplots()

ax.add_patch(circle)
xs_right = [rlo, rhi]
ys = [dz, maxh]
xs_left = [-rlo, -rhi]
ax.plot(xs_right, ys)
ax.plot(xs_left, ys)
ax.plot([-rhi, rhi],[maxh, maxh])
#fig.savefig('plotcircle.png')
plt.show()
