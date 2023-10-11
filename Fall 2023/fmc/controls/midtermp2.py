# AEEM4042 Controls Midterm Part 2
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


# define system
anim = animation.ballBeamAnim(0.6, False)
bb = dynamics.ballBeamDynamics()
uf = sig.signalGenerator(0.01, 0.1)

# define initial states
state = np.array([[P.z0], [P.theta0], [P.zdot0], [P.thetadot0]])

# initial force
u = 13.23

# sim loop
t = P.t_start
while t < P.t_end:
    # add force
    u += uf.sin(t)
    
    # update dynamics
    bb.update(u)
    
    # update animation
    anim.update(bb.state)
    
    # draw stuff and increment time
    plt.pause(0.1)
    t += P.Ts
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()