# AEEM4042 Controls Assignment #5 Part F8E
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.VTOLsimparams2 as P
import lib.VTOLAnim as animation
import lib.VTOLDynamics2 as dynamics
import keyboard
import lib.vtolPD2 as ctr
import lib.signalGenerator as sig
import numpy as np


# define system
vtol = dynamics.VTOLDynamics()
anim = animation.VTOLAnim(limits=10, flag=True)
ctr = ctr.controller()
zsig = sig.signalGenerator(0.25, 0.08)

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
zt = 3.
ht = 0.5
hts = [ht]
zts = [zt]
ts = [0]

# sim loop
t = P.t_start
simtimes = [t]
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        if t > 12.5:
            ht = 0.5
            zt = 3 + zsig.square(t)
    
        fr, fl = ctr.update(zt, ht, vtol.state)
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