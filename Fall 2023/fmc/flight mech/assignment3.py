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
import lib.UAVparams as P
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj, f18_verts
import numpy as np


# define stuff
uav = dynamics.UAVdynamics()
anim = animation.animation(limits=10, alpha=0.35, flag=False)
aero = aero.UAVaero()
Vs = np.array([[5],[2],[0]])        # steady wind m/s
wind = wind.wind(Vs)

# create vehicle
verts = sampleUAV_verts
verts = f18_verts
obj = sampleUAV_obj
obj = None
faces = ["b"]

# initial state
state = P.states0
deltaa = 0
deltae = 0
deltar = 0
deltat = 1
Va = np.sqrt(state[3][0]**2 + state[4][0]**2 + state[5][0]**2)

# add subplots
throttle = plt.figure(1).add_subplot(1, 10, 1)

# initialize the simulation time
sim_time = SIM.start_time
t_next_plot = sim_time + SIM.ts_plotting

# main simulation loop
print("Press Q to exit...")
while sim_time < SIM.end_time:
    
    # keyboard inputs
    if keyboard.is_pressed("down arrow"): deltae -= np.deg2rad(1)
    if keyboard.is_pressed("up arrow"): deltae += np.deg2rad(1)
    if keyboard.is_pressed("right arrow"): deltaa += np.deg2rad(0.5); deltar -= np.deg2rad(0.25)
    if keyboard.is_pressed("left arrow"): deltaa -= np.deg2rad(0.5); deltar += np.deg2rad(0.25)
    if keyboard.is_pressed("space"): deltae = 0; deltaa = 0; deltar = 0
    if keyboard.is_pressed("shift"):
        if deltat < 1: deltat += 0.05
    if keyboard.is_pressed("left control"):
        if deltat > 0: deltat -= 0.05
    
    # update everything
    Va, alpha, beta = wind.windout(state, Va)
    fx, fy, fz = aero.forces(state, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
    l, m, n = aero.moments(state, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
    state = uav.update(fx, fy, fz, l, m, n)
    pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = state.flatten()
    anim.update(verts, pn, pe, pd, phi, theta, psi, obj, faces)
    
    # plotting
    if sim_time > t_next_plot:
        throttle.clear(); throttle.bar(0, deltat); throttle.set_ylim(0, 1); throttle.set_xticklabels(""); throttle.set_xticks([])
        t_next_plot = sim_time + SIM.ts_plotting
    
    # kill you
    if pd >= 0:
        plt.figure(1).text(x=0.25, y=0.5, s="You Crashed", fontsize=100)
        plt.pause(0.01)
        sim_time = SIM.end_time
    
    # increment time
    sim_time += SIM.ts_simulation

    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()