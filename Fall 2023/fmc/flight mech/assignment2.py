# AEEM4012 Flight Mechanics Assignment #2
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\flight mech\\lib\\")
import matplotlib.pyplot as plt
import lib.simparams as SIM
import lib.animation as animation
import lib.UAVdynamics as dynamics
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj
from lib.signalGenerator import signalGenerator
import numpy as np


# define stuff
uav = dynamics.UAVdynamics()
anim = animation.animation(limits=10, alpha=0.5, flag=False)
forces = signalGenerator(1500., 2)
moments = signalGenerator(10., 2)

# add subplots
forcesp = anim.fig.add_subplot(231)
fpp = forcesp.get_position(); fpp.x0 -= 0.05; fpp.x1 -= 0.05; forcesp.set_position(fpp)
momentsp = anim.fig.add_subplot(234)
mpp = momentsp.get_position(); mpp.x0 -= 0.05; mpp.x1 -= 0.05; momentsp.set_position(mpp)
posp = anim.fig.add_subplot(433)
ppp = posp.get_position(); ppp.x0 += 0.1; ppp.x1 += 0.1; posp.set_position(ppp)
velp = anim.fig.add_subplot(436)
vpp = velp.get_position(); vpp.x0 += 0.1; vpp.x1 += 0.1; velp.set_position(vpp)
angp = anim.fig.add_subplot(439)
app = angp.get_position(); app.x0 += 0.1; app.x1 += 0.1; angp.set_position(app)
ratep = anim.fig.add_subplot(4, 3, 12)
rpp = ratep.get_position(); rpp.x0 += 0.1; rpp.x1 += 0.1; ratep.set_position(rpp)

# create vehicle
verts = sampleUAV_verts
obj = sampleUAV_obj
faces = ["b"]

# initialize the simulation time
sim_time = SIM.start_time

# inputs
fxs = [0.]
fys = [0.]
fzs = [0.]
ls = [0.]
ms = [0.]
ns = [0.]
ts = [sim_time]
pns = [0.]
pes = [0.]
pds = [0.]
us = [0.]
vs = [0.]
ws = [0.]
phis = [0.]
thetas = [0.]
psis = [0.]
ps = [0.]
qs = [0.]
rs = [0.]

# main simulation loop
print("Press Q to exit...")
while sim_time < SIM.end_time:
    # set forces
    if sim_time <= 0.5:
        fx = forces.sin(sim_time)
        fy = 0
        fz = 0
        l = 0
        m = 0
        n = 0
    elif sim_time <= 1:
        fx = 0
        fy = forces.sin(sim_time)
    elif sim_time <= 1.5:
        fy = 0
        fz = -forces.sin(sim_time)
    elif sim_time <= 2:
        fz = 0
        l = moments.sin(sim_time)
    elif sim_time <= 2.5:
        l = 0
        m = moments.sin(sim_time)
    elif sim_time <= 3:
        m = 0
        n = moments.sin(sim_time)
    else:
        n = 0
    
    # store daterp
    fxs.append(fx)
    fys.append(fy)
    fzs.append(fz)
    ls.append(l)
    ms.append(m)
    ns.append(n)
    ts.append(sim_time)
    
    # update everything
    y = uav.update(fx, fy, fz, l, m, n)
    pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = y.flatten()
    anim.update(verts, pn, pe, pd, phi, theta, psi, obj, faces)
    
    # store more
    pns.append(pn)
    pes.append(pe)
    pds.append(-pd)
    us.append(u)
    vs.append(v)
    ws.append(w)
    phis.append(phi)
    thetas.append(theta)
    psis.append(psi)
    ps.append(p)
    qs.append(q)
    rs.append(r)
    
    # plot
    forcesp.clear(); momentsp.clear(); posp.clear(); velp.clear(); angp.clear(); ratep.clear()
    forcesp.plot(ts, fxs, label="$F_x$"); forcesp.plot(ts, fys, label="$F_y$"); forcesp.plot(ts, fzs, label="$F_z$")
    momentsp.plot(ts, ls, label="$l$"); momentsp.plot(ts, ms, label="$m$"); momentsp.plot(ts, ns, label="$n$")
    forcesp.grid(); forcesp.legend(loc="upper right"); forcesp.set_title("Force Input"); forcesp.set_ylabel("Force (N)")
    momentsp.grid(); momentsp.legend(loc="upper left"); momentsp.set_title("Moment Input"); momentsp.set_ylabel("Moment (Nm)")
    posp.plot(ts, pns, label="North"); posp.plot(ts, pes, label="East"); posp.plot(ts, pds, label="Height")
    posp.set_ylabel("Position (m)"); posp.grid(); posp.legend(loc="lower right")
    velp.plot(ts, us, label="$U$"); velp.plot(ts, vs, label="$V$"); velp.plot(ts, ws, label="$W$")
    velp.set_ylabel("Velocity (m/s)"); velp.grid(); velp.legend(loc="upper right")
    angp.plot(ts, np.degrees(phis), label="$\phi$"); angp.plot(ts, np.degrees(thetas), label="$\\theta$")
    angp.plot(ts, np.degrees(psis), label="$\psi$")
    angp.set_ylabel("Angle (deg)"); angp.grid(); angp.legend(loc="upper left")
    ratep.plot(ts, np.degrees(ps), label="p"); ratep.plot(ts, np.degrees(qs), label="q"); ratep.plot(ts, np.degrees(rs), label="r")
    ratep.set_ylabel("Rate (deg/s)"); ratep.grid(); ratep.legend(loc="upper left")
    
    # increment time
    plt.pause(0.01)
    sim_time += SIM.ts_simulation
    
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()