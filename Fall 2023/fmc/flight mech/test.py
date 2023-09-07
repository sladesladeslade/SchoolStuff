# AEEM4012 Flight Mechanics Assignment #1
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\")
import numpy as np
from lib.animation import animation
from lib.sliders import sliders
import lib.simparams as SIM
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj


# create vehicle
verts = sampleUAV_verts
verts2 = sampleUAV_verts
obj = sampleUAV_obj
obj2 = sampleUAV_obj
faces = ["b"]
faces2 = ["g"]

# init animation class
plane = animation(limits=10, alpha=0.5)

# ICs
state = np.array([0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
state2 = np.array([0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# init sliders
slider = sliders()

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Q to exit...")
n = state[0]
d = state[2]
phi = state[6]
theta = state[7]
psi = state[8]

while sim_time < SIM.end_time:
    # reading from sliders
    phi = slider.roll_slider.val
    theta = slider.pitch_slider.val
    psi = slider.yaw_slider.val
    
    # updating from those vals
    plane.update(verts, n, -5, d, phi, theta, psi, obj, facecolors=faces)
    plane.update(verts2, n, 5, d, phi, theta, psi, obj2, facecolors=faces2)
    
    # increment time
    sim_time += SIM.ts_simulation
    
    # stop on q
    if keyboard.is_pressed("q"): break