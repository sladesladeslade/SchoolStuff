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
anim = animation.ballBeamAnim(0.6, False)
ball = sig.signalGenerator(0.5, 0.1)
beam = sig.signalGenerator(np.deg2rad(45), 0.01)

# define initial states
state = np.array([[P.z0], [P.theta0], [P.zdot0], [P.thetadot0]])

# sim loop
t = P.t_start
while t < P.t_end:
    state[0][0] = ball.sin(t)
    if state[0][0] < 0:
        state[0][0] *= -1
    state[1][0] = beam.sin(t)
        
    anim.update(state)
    
    plt.pause(0.01)
    t += P.Ts
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()