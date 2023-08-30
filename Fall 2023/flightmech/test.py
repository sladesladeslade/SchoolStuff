import numpy as np
from animation import animation
from sliders import sliders
import simparams as SIM


# make cube
w = 5
verts = np.array([[w/2, -w/2, -w/2],
        [-w/2, -w/2, -w/2],
        [-w/2, w/2, -w/2],
        [w/2, w/2, -w/2],
        [w/2, -w/2, w/2],
        [-w/2, -w/2, w/2],
        [-w/2, w/2, w/2],
        [w/2, w/2, w/2]])
def obj(v):
    return np.array([[v[0], v[1], v[2], v[3]],
                     [v[4], v[5], v[6], v[7]],
                     [v[2], v[3], v[7], v[6]],
                     [v[1], v[0], v[4], v[5]],
                     [v[0], v[3], v[7], v[4]],
                     [v[2], v[6], v[5], v[1]]])
faces = ['g', 'r', 'r', 'r', 'r','r']

state = np.array([[0], [0], [-1], [0], [0], [0], [0], [0], [0], [0], [0], [0]])
cube = animation()
slider = sliders()

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Q to exit...")
n=state[0,0]
e=state[1,0]
d=state[2,0]
phi=state[6,0]
theta=state[7,0]
psi=state[8,0]

# loop and sim
while sim_time < SIM.end_time:
    # reading from sliders
    phi = slider.roll_slider.val
    theta = slider.pitch_slider.val
    psi = slider.yaw_slider.val
    
    # updating from those vals
    cube.update(verts, obj, n, e, d, phi, theta, psi, facecolors=faces)
    
    # -------increment time-------------
    sim_time += SIM.ts_simulation