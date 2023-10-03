# AEEM4042 Controls Assignment #1 Part 2
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.VTOLsimparams as P
import lib.VTOLAnim as animation
import lib.VTOLDynamics as dynamics
import keyboard
import lib.vtolPD as ctr


# define system
vtol = dynamics.VTOLDynamics()
anim = animation.VTOLAnim(limits=10, flag=True)
ctr = ctr.controller()

# add subplots
zes = anim.fig.add_subplot(2, 2, 2)
fes = anim.fig.add_subplot(2, 2, 4)
zes.set_ylabel("h (m)")
fes.set_ylabel("Force (N)")
simtimes = []
hs = []
fs = []
targets = []
target = 0.25

# sim loop
t = P.t_start
while t < P.t_end:
    if t <= 2:
        target = 0.25
    elif t <= 30:
        target = 5
    elif t <= 50:
        target = 8
    elif t <= 65:
        target = 0.5
    else:
        target = 0.25
    
    F = ctr.update(target, vtol.state)
    y = vtol.update(F/2., F/2.)
    
    # update daterp
    simtimes.append(t)
    hs.append(y[1][0])
    fs.append(F)
    targets.append(target)
    
    # update anim
    anim.update(vtol.state)
    zes.clear()
    fes.clear()
    zes.plot(simtimes, hs, label="State")
    fes.plot(simtimes, fs)
    zes.plot(simtimes, targets, label="Target")
    zes.legend(loc="lower right")
    zes.grid()
    fes.grid()
    zes.set_ylabel("z (m)")
    fes.set_ylabel("Total Force (N)")
    plt.pause(0.01)
    
    t += P.Ts
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()