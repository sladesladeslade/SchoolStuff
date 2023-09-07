# AEEM4042 Controls Assignment #1 Part 1
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.mssimparams as P
import lib.massSpringAnim as animation
import lib.massSpringDynamics as dynamics
import keyboard


# define system
msd = dynamics.massSpringDynamics()
anim = animation.massSpringAnim(limits=1, flag=True)

# add subplots
zes = anim.fig.add_subplot(2, 2, 2)
fes = anim.fig.add_subplot(2, 2, 4)
zes.set_ylabel("z (m)")
fes.set_ylabel("Force (N)")
simtimes = []
zs = []
fs = []

# simulation loop
t = P.t_start
while t < P.t_end:
    # do dynamics
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        u = 1
        y = msd.update(u)
        t += P.Ts
        
    # update daterp
    simtimes.append(t)
    zs.append(y)
    fs.append(u)
    
    # update animation and plots
    anim.update(msd.state)
    zes.clear()
    fes.clear()
    zes.plot(simtimes, zs)
    fes.plot(simtimes, fs)
    zes.grid()
    fes.grid()
    zes.set_ylabel("z (m)")
    fes.set_ylabel("Force (N)")
    plt.pause(0.001)
    
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()