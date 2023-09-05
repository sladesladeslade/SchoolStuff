# AEEM4012 Flight Mechanics Assignment #1
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\")
import numpy as np
from animation import animation
from sliders import sliders
import simparams as SIM
import keyboard
from hangar import sampleUAV_verts, sampleUAV_obj, visionJet_verts, x59_verts
from signalGenerator import signalGenerator


# create vehicle
# verts = sampleUAV_verts
# obj = sampleUAV_obj
verts = x59_verts
obj = None
faces = ["b"]

# init animation class
plane = animation(limits=10, alpha=0.25)

# ICs
state = np.array([0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# init sliders
slider = sliders()

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Q to exit...")
n = state[0]
e = state[1]
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
    plane.update(verts, n, e, d, phi, theta, psi, obj, facecolors=faces)
    
    # increment time
    sim_time += SIM.ts_simulation
    
    # stop on q
    if keyboard.is_pressed("q"): break