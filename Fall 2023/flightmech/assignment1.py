# AEEM4012 Flight Mechanics Assignment #1
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
from animation import animation
from sliders import sliders
import simparams as SIM
import keyboard
from hangar import SampleUAV_verts, SampleUAV_obj


# create vehicle
verts = SampleUAV_verts
obj = SampleUAV_obj
faces = ["g"]

# init animation class
plane = animation()

# ICs
state = np.array([[0], [0], [-1], [0], [0], [0], [0], [0], [0], [0], [0], [0]])

# init sliders
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

while sim_time < SIM.end_time:
    # reading from sliders
    phi = slider.roll_slider.val
    theta = slider.pitch_slider.val
    psi = slider.yaw_slider.val
    
    # updating from those vals
    plane.update(verts, obj, n, e, d, phi, theta, psi, facecolors=faces)
    
    # -------increment time-------------
    sim_time += SIM.ts_simulation
    
    # stop on q
    if keyboard.is_pressed("q"):
        break