# AEEM4042 Controls Assignment #5 Part F10
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.ass5.VTOLsimparams as P
import lib.VTOLAnim as animation
import lib.ass5.VTOLDynamics as dynamics
import keyboard
import lib.ass5.VTOLPID as ctr
import numpy as np


# define system
vtol = dynamics.VTOLDynamics()
anim = animation.VTOLAnim(limits=10, flag=True)
ctr = ctr.controller(10, 0.05, False)

# set controller
ctr.kpt = 0.3721
ctr.kdt = 0.1913
ctr.kpz = -0.0077095
ctr.kdz = -0.032858
ctr.kph = 0.11345
ctr.kdh = 0.5835
ctr.kih = 0.
ctr.kit = 0.
ctr.kiz = 0.

# add subplots
zes = anim.fig.add_subplot(322)
fes = anim.fig.add_subplot(326)
tes = anim.fig.add_subplot(324)
zes.set_ylabel("Location (m)")
fes.set_ylabel("Force (N)")
tes.set_ylabel(r"$\theta$")
hs = [P.h0]
zs = [P.z0]
frs = [0]
fls = [0]
zt = 0.
ht = 0.25
hts = [ht]
zts = [zt]
ts = [0]

# sim loop
t = P.t_start
simtimes = [t]
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        if t > 2:
            ht = 0.25
            zt = 0.

        zee, ach, thet = vtol.h().flatten()
        fr, fl = ctr.update(ht, zt, zee, ach, thet)
        y = vtol.update(fr, fl)
        t += P.Ts
    
    # update daterp
    simtimes.append(t)
    hs.append(y[1][0])
    zs.append(y[0][0])
    ts.append(np.rad2deg(y[2][0]))
    frs.append(fr)
    fls.append(fl)
    zts.append(zt)
    hts.append(ht)
    
    # update anim
    anim.update(vtol.state)
    zes.clear()
    fes.clear()
    tes.clear()
    tes.plot(simtimes, ts)
    zes.plot(simtimes, hs, color="tab:blue", label="Height")
    zes.plot(simtimes, zs, color="tab:orange", label="Z")
    fes.plot(simtimes, frs, label="Fr")
    fes.plot(simtimes, fls, label="Fl")
    zes.plot(simtimes, hts, color="tab:blue", linestyle="--", label="Target Height")
    zes.plot(simtimes, zts, color="tab:orange", linestyle="--", label="Target Z")
    zes.legend(loc="lower right")
    fes.legend(loc="lower right")
    zes.grid()
    fes.grid()
    tes.grid()
    zes.set_ylabel("h (m)")
    fes.set_ylabel("Total Force (N)")
    tes.set_ylabel(r"$\theta$ (deg)")
    
    plt.pause(0.01)
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()