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
from hangar import SampleUAV_verts, SampleUAV_obj, visionJet_verts, visionJet_obj
from signalGenerator import signalGenerator


# create vehicle
verts = SampleUAV_verts
obj = SampleUAV_obj
# verts = visionJet_verts
# obj = visionJet_obj
faces = ["b"]

# init animation class
plane = animation(limits=5)

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

# define signal
sig = signalGenerator(np.deg2rad(90), 0.1)

while sim_time < SIM.end_time:
    # reading from sliders
    phi = slider.roll_slider.val
    theta = slider.pitch_slider.val
    psi = slider.yaw_slider.val
    
    # updating from those vals
    theta = sig.sin(sim_time)
    plane.update(verts, obj, n, e, d, phi, theta, psi, facecolors=faces)
    
    # increment time
    sim_time += SIM.ts_simulation
    
    # stop on q
    if keyboard.is_pressed("q"): break