# AEEM4012 Flight Mechanics Assignment #1
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\flight mech\\lib\\")
import numpy as np
from lib.animation import animation
from lib.sliders import sliders
import lib.simparams as SIM
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj, visionJet_verts, x59_verts
from lib.signalGenerator import signalGenerator


# create vehicle
# verts = sampleUAV_verts
# obj = sampleUAV_obj
verts = x59_verts
obj = None
faces = ["b"]

# init animation class
plane = animation(limits=10, alpha=0.25, flag=True)

# ICs
state = np.array([0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# set up signal generator
angles = signalGenerator(np.deg2rad(45), 1)
translations = signalGenerator(5, 1)

# initialize the simulation time
sim_time = SIM.start_time

# add subplots
simtimes = []
phis = []
thetas = []
psis = []
ns = []
es = []
ds = []
anglesp = plane.fig.add_subplot(2, 2, 2)
transp = plane.fig.add_subplot(2, 2, 4)
anglesp.set_title("Rotation")
transp.set_title("Translation")

# main simulation loop
print("Press Q to exit...")
n = state[0]
e = state[1]
d = state[2]
phi = state[6]
theta = state[7]
psi = state[8]

while sim_time < SIM.end_time:
    # read in positions and rotations
    if sim_time <= 1:
        phi = 0
    elif sim_time <= 2:
        phi = angles.sin(sim_time)
    elif sim_time <= 3:
        theta = angles.sin(sim_time)
    elif sim_time <= 4:
        psi = angles.sin(sim_time)
    elif sim_time <= 5:
        n = translations.sin(sim_time)
    elif sim_time <= 6:
        e = translations.sin(sim_time)
    elif sim_time <= 7:
        d = translations.sin(sim_time)
        
    # store data
    simtimes.append(sim_time)
    phis.append(phi)
    thetas.append(theta)
    psis.append(psi)
    ns.append(n)
    es.append(e)
    ds.append(d)

    # plot it all
    anglesp.clear()
    transp.clear()
    anglesp.plot(simtimes, np.rad2deg(phis), "r-", label="$\phi$")
    anglesp.plot(simtimes, np.rad2deg(thetas), "g-", label="$\\theta$")
    anglesp.plot(simtimes, np.rad2deg(psis), "b-", label="$\psi$")
    transp.plot(simtimes, ns, "r-", label="North")
    transp.plot(simtimes, es, "g-", label="East")
    transp.plot(simtimes, ds, "b-", label="Height")
    anglesp.legend(loc="upper right")
    anglesp.grid()
    anglesp.set_ylim(bottom=-60, top=60)
    transp.legend(loc="upper right")
    transp.grid()
    transp.set_ylim(bottom=-6, top=6)
    
    # updating from those vals
    plane.update(verts, n, e, d, phi, theta, psi, obj, facecolors=faces)
    
    # increment time
    sim_time += SIM.ts_simulation
    
    # stop on q
    if keyboard.is_pressed("q"): break