import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# make lbm importable
sys.path.append(os.path.abspath(os.path.dirname(__file__) + "/.."))
from lbm.streamcollide import set_equi
from lbm.util import CT, v, v_inv, t
from lbm.time_evolve import  evolve
from lbm.init import init


# set the domain variables
nx = 100
ny = 100
omega = np.array([1.0])       # relaxation parameter
gravity = np.array([0, 0])
total_timesteps = 100
# get the variables
fin, fout, equi, fdist, inlet, u, rho, mass, cell_type = init(nx, ny)

# make a fluid square
radius1 = 10
radius2 = 10
x1 = 50
y1 = 50

x2 = 65
y2 = 50


for ix in range(mass.shape[0]):
    for iy in range(mass.shape[1]):
        distsq = (ix - x1)**2 + (iy - y1)**2
        if (distsq < radius1**2):
            cell_type[ix, iy] = CT.FLUID
            mass[ix, iy] = 1
        elif (distsq >= radius1**2) and (distsq < (radius1+1)**2):
            cell_type[ix, iy] = CT.INTERFACE
            mass[ix,iy] = 1 - (distsq**(1/2) - radius1)      
'''
for ix in range(mass.shape[0]):
    for iy in range(mass.shape[1]):
        distsq = (ix - x2)**2 + (iy - y2)**2
        if (distsq < radius2**2):
            cell_type[ix, iy] = CT.FLUID
            mass[ix, iy] = 1
        elif (distsq >= radius2**2) and (distsq < (radius2+1)**2):
            cell_type[ix, iy] = CT.INTERFACE
            mass[ix,iy] = 1 - (distsq**(1/2) - radius2)   
'''

# initialize fin

for ix in range(rho.shape[0]):
    for iy in range(rho.shape[1]):
        set_equi(ix, iy, equi, rho[ix, iy], u[:, ix, iy], v, t)
fin = equi.copy()

# evolve the dam break
u, cell_types = evolve(total_timesteps, fin, fout, equi, fdist, inlet, u, rho, mass, cell_type, omega, v, v_inv, t, gravity)


# animate
def animate(i):
    im.set_array((cell_types[i] & (CT.FLUID | CT.INTERFACE)) != 0)
    return [im]


fig = plt.figure()
im = plt.imshow((cell_types[-1] & (CT.INTERFACE)) != 0, cmap=plt.cm.Blues)
plt.savefig('test_cell.png')
plt.close()
im = plt.imshow(mass, cmap=plt.cm.Blues)
plt.savefig('test_mass.png')
plt.close()

