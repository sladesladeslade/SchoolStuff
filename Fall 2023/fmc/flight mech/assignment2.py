# AEEM4012 Flight Mechanics Assignment #2
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\flight mech\\lib\\")
import matplotlib.pyplot as plt
import lib.simparams as SIM
import lib.animation as animation
import lib.UAVdynamics as dynamics
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj


# define stuff
uav = dynamics.UAVdynamics()
anim = animation.animation(limits=10, alpha=0.5, flag=False)

# create vehicle
verts = sampleUAV_verts
obj = sampleUAV_obj
faces = ["b"]

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Q to exit...")
while sim_time < SIM.end_time:
    fx = 0.
    fy = 0.
    fz = 0.
    l = 0.
    m = 0.
    n = 0.
    y = uav.update(fx, fy, fz, l, m, n)
    anim.update(verts, y[0][0], y[1][0], y[2][0], y[6][0], y[7][0], y[8][0], obj, faces)
    plt.pause(0.01)
    
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()