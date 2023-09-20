# AEEM4012 Flight Mechanics Assignment #3
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\flight mech\\lib\\")
import matplotlib.pyplot as plt
import lib.simparams as SIM
import lib.animation as animation
import lib.UAVdynamics as dynamics
import lib.UAVaero as aero
import lib.wind as wind
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj
from lib.signalGenerator import signalGenerator
import numpy as np


# define stuff
uav = dynamics.UAVdynamics()
anim = animation.animation(limits=10, alpha=0.5, flag=False)
aero = aero.UAVaero()
Vs = np.array([[0],[0],[0]])        # steady wind m/s
Vg = np.array([[0],[0],[0]])        # max gust m/s
wind = wind.wind(Vs)

# create vehicle
verts = sampleUAV_verts
obj = sampleUAV_obj
faces = ["b"]

# initial state
state = np.array([[0.],[0.],[-10.],[10.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]])
deltaa = 0
deltae = 0
deltar = 0
deltat = 1
Va = np.sqrt(state[3][0]**2 + state[4][0]**2 + state[5][0]**2)

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Q to exit...")
while sim_time < SIM.end_time:
    # keyboard inputs
    if keyboard.is_pressed("down arrow"): deltae -= np.deg2rad(0.5)
    if keyboard.is_pressed("up arrow"): deltae += np.deg2rad(0.5)
    if keyboard.is_pressed("right arrow"): deltaa += np.deg2rad(0.5)
    if keyboard.is_pressed("left arrow"): deltaa -= np.deg2rad(0.5)
    
    # update everything
    Va, alpha, beta = wind.windout(state, Va)
    fx, fy, fz = aero.forces(state, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
    l, m, n = aero.moments(state, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
    state = uav.update(fx, fy, fz, l, m, n)
    pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = state.flatten()
    anim.update(verts, pn, pe, pd, phi, theta, psi, obj, faces)
    
    # increment time
    sim_time += SIM.ts_simulation
    
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()