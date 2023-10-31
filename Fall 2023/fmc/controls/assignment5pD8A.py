# AEEM4042 Controls Assignment #5 Part D8A
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.mssimparams as P
import lib.massSpringAnim as animation
import lib.massSpringDynamics as dynamics
import keyboard
import lib.msdPD as ctr
import numpy as np


# define system
msd = dynamics.massSpringDynamics()
anim = animation.massSpringAnim(limits=2, flag=True)
ctr = ctr.controller(Fmax=np.inf)
ctr.kD = 7.2
ctr.kP = 3.05

# add subplots
zes = anim.fig.add_subplot(2, 2, 2)
fes = anim.fig.add_subplot(2, 2, 4)
zes.set_ylabel("z (m)")
fes.set_ylabel("Force (N)")
simtimes = []
zs = []
fs = []
targets = []

# initials
z = 0
u = 0
target = 0

# simulation loop
t = P.t_start
while t < P.t_end:
    if t <= 2:
        target = 0
    elif t <= 20:
        target = 1
    elif t <= 30:
        target = 1.5
    elif t <= 40:
        target = 0.5
    else:
        target = 0
        
    u = ctr.update(target, msd.state)
    z = msd.update(u)
    
    # update daterp
    simtimes.append(t)
    zs.append(z)
    fs.append(u)
    targets.append(target)
    
    # update animation and plots
    anim.update(msd.state)
    zes.clear()
    fes.clear()
    zes.plot(simtimes, zs, label="State")
    fes.plot(simtimes, fs)
    zes.plot(simtimes, targets, label="Target")
    zes.legend(loc="lower right")
    zes.grid()
    fes.grid()
    zes.set_ylabel("z (m)")
    fes.set_ylabel("Force (N)")
    plt.pause(0.01)
    
    t += P.Ts
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()