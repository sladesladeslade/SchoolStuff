# AEEM4042 Controls Assignment #1
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.simparams as P
import lib.massSpringAnim as animation
import lib.massSpringDynamics as dynamics
import keyboard
import lib.dataPlotter as dp
from lib.signalGenerator import signalGenerator

# define system
msd = dynamics.massSpringDynamics()
anim = animation.massSpringAnim(limits=10)
reference = signalGenerator(1, 1)
dataplot = dp.dataPlotter()

# simulation loop
t = P.t_start
while t < P.t_end:
    # do dynamics
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        u = 1
        y = msd.update(u)
        t += P.Ts
    
    # update animation
    anim.update(msd.state)
    dataplot.update(t, r, msd.state, u)
    plt.pause(0.001)
    
    if keyboard.is_pressed("q"): break