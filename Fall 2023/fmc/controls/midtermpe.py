# AEEM4042 Controls Midterm Part E
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.bbsimparams as P
import lib.ballBeamAnim as animation
import lib.ballBeamDynamics as dynamics
import numpy as np
import keyboard
import lib.signalGenerator as sig
import lib.bbPD as ctr


# define system
anim = animation.ballBeamAnim(0.6, True)
bb = dynamics.ballBeamDynamics()
ctr = ctr.controller()
tdev = sig.signalGenerator(0.15, 0.01)

# set up plots
forcep = anim.fig.add_subplot(522)
zp = anim.fig.add_subplot(524); zdotp = anim.fig.add_subplot(526)
tp = anim.fig.add_subplot(528); tdotp = anim.fig.add_subplot(5, 2, 10)
forces = []
zs = []
zdots = []
thetas = []
tdots = []
ts = []
targets = []

# define initial states
state = np.array([[P.z0], [P.theta0], [P.zdot0], [P.thetadot0]])

# initial force
u = 13.23
target = 0.25

# sim loop
t = P.t_start
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        # update target
        target = 0.25 + tdev.square(t)
        
        # add force
        u = ctr.update(target, bb.state)
        
        # update dynamics
        bb.update(u)
        
        # update animation
        anim.update(bb.state)
        
        # increment time
        t += P.Ts
        
    # do plotting
    # get new vals
    z, th, zd, td = bb.state.flatten()
    forces.append(u); zs.append(z); zdots.append(zd); targets.append(target)
    thetas.append(np.rad2deg(th)); tdots.append(np.rad2deg(td)); ts.append(t)
    # plot and format
    forcep.clear(); zp.clear(); zdotp.clear(); tp.clear(); tdotp.clear()
    zp.plot(ts, zs, label="Actual"); zp.plot(ts, targets, "r-", label="Target"); zp.legend(loc="lower left")
    zp.set_xlim((0, ts[-1])); zp.set_ylim((0, 0.5)); zp.set_ylabel("z (m)"); zp.grid()
    forcep.plot(ts, forces); forcep.set_xlim((0, ts[-1])); forcep.set_ylim((0, 15))
    forcep.set_ylabel("Force Input (N)"); forcep.grid()
    zdotp.plot(ts, zdots); zdotp.set_xlim((0, ts[-1])); zdotp.set_ylabel("$\dot{z}$ (m/s)"); zdotp.grid()
    tp.plot(ts, thetas); tp.set_xlim((0, ts[-1])); tp.set_ylabel("$\\theta$ (deg)"); tp.grid()
    tdotp.plot(ts, tdots); tdotp.set_xlim((0, ts[-1])); tdotp.set_ylabel("$\dot{\\theta}$ (deg/s)"); tdotp.grid()
    
    # draw stuff and increment time
    plt.pause(0.01)
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()