# AEEM4042 Controls Midterm Part 1
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.bbsimparams as P
import lib.ballBeamAnim as animation
import numpy as np
import keyboard
import lib.signalGenerator as sig


# define system
anim = animation.ballBeamAnim(0.6, True)
ball = sig.signalGenerator(0.5, 0.1)
beam = sig.signalGenerator(np.deg2rad(45), 0.01)

# set up plots
zp = anim.fig.add_subplot(222)
tp = anim.fig.add_subplot(224)
buffer = 100
zs = np.zeros(buffer)
thetas = np.zeros(buffer)
ts = np.zeros(buffer)

# define initial states
state = np.array([[P.z0], [P.theta0], [P.zdot0], [P.thetadot0]])

# sim loop
t = P.t_start
while t < P.t_end:
    # setting positions
    state[0][0] = ball.sin(t)
    if state[0][0] < 0:
        state[0][0] *= -1
    state[1][0] = beam.sin(t)
    
    # store data
    zs = np.concatenate((zs[1:], [state[0][0]]))
    thetas = np.concatenate((thetas[1:], np.rad2deg([state[1][0]])))
    ts = np.concatenate((ts[1:], [t]))
    
    # update animation
    anim.update(state)
    
    # update plots
    zp.clear(); tp.clear()
    zp.plot(ts, zs); tp.plot(ts, thetas)
    zp.set_xlabel("Time (s)"); zp.set_ylabel("Z (m)"); zp.set_ylim((0, 0.5))
    tp.set_xlabel("Time (s)"); tp.set_ylabel("$\\theta$ (deg)"); tp.set_ylim((-60, 60))
    zp.set_xlim((ts[0], ts[-1])); tp.set_xlim((ts[0], ts[-1]))
    zp.grid(); tp.grid()
    
    # draw stuff and increment time
    plt.pause(0.001)
    t += P.Ts
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()