# AEEM4012 Flight Mechanics Assignment #3
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\flight mech\\lib\\")
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import lib.simparams as SIM
import lib.animation as animation
import lib.UAVdynamics as dynamics
import lib.UAVlinaero as aero
import lib.wind as wind
import lib.UAVparams as P
import lib.compute_trim as comp
import keyboard
from lib.hangar import sampleUAV_verts, sampleUAV_obj, f18_verts
import numpy as np
import warnings; warnings.filterwarnings("ignore", category=UserWarning, module="control")


# define stuff
uav = dynamics.UAVdynamics()
anim = animation.animation(limits=10, alpha=0.35, flag=False)
aero = aero.UAVaero()
Vs = np.array([[0],[0],[0]])        # steady wind m/s
wind = wind.wind(Vs)
trim = comp.ComputeTrim()

# create vehicle
verts = sampleUAV_verts
# verts = f18_verts
obj = sampleUAV_obj
# obj = None
faces = ["b"]

# add subplots
throttle = plt.figure(1).add_subplot(1, 20, 1); tpp = throttle.get_position(); tpp.x0-=0.1; tpp.x1-=0.1
throttle.set_position(tpp)
deflections = anim.fig.add_subplot(331)
dpp = deflections.get_position(); dpp.x0 -= 0.015; dpp.x1 -= 0.05; deflections.set_position(dpp)
force = anim.fig.add_subplot(334)
fpp = force.get_position(); fpp.x0 -= 0.015; fpp.x1 -= 0.05; force.set_position(fpp)
moment = anim.fig.add_subplot(337)
mpp = moment.get_position(); mpp.x0 -= 0.015; mpp.x1 -= 0.05; moment.set_position(mpp)
posp = anim.fig.add_subplot(433)
ppp = posp.get_position(); ppp.x0 += 0.125; ppp.x1 += 0.09; posp.set_position(ppp)
velp = anim.fig.add_subplot(436)
vpp = velp.get_position(); vpp.x0 += 0.125; vpp.x1 += 0.09; velp.set_position(vpp)
angp = anim.fig.add_subplot(439)
app = angp.get_position(); app.x0 += 0.125; app.x1 += 0.09; angp.set_position(app)
ratep = anim.fig.add_subplot(4, 3, 12)
rpp = ratep.get_position(); rpp.x0 += 0.125; rpp.x1 += 0.09; ratep.set_position(rpp)

# set up data arrays
buffer = int(1/SIM.ts_plotting)
simtimes = np.zeros(buffer)
defsa = np.zeros(buffer)
defse = np.zeros(buffer)
defsr = np.zeros(buffer)
forces0 = np.zeros(buffer); forces1 = np.zeros(buffer); forces2 = np.zeros(buffer)
moments0 = np.zeros(buffer); moments1 = np.zeros(buffer); moments2 = np.zeros(buffer)
poss0 = np.zeros(buffer); poss1 = np.zeros(buffer); poss2 = np.zeros(buffer)
vels0 = np.zeros(buffer); vels1 = np.zeros(buffer); vels2 = np.zeros(buffer)
angs0 = np.zeros(buffer); angs1 = np.zeros(buffer); angs2 = np.zeros(buffer)
rates0 = np.zeros(buffer); rates1 = np.zeros(buffer); rates2 = np.zeros(buffer)

# initialize the simulation time
sim_time = SIM.start_time

# targets
Vat = 35.
Y = np.deg2rad(0.1)
R = np.inf
xtrim, utrim = trim.compute_trim(Vat, Y, R)
deltae, deltat, deltaa, deltar = utrim.flatten()
# deltae = -deltae
print("--- Trim Conditions ---")
print(f"E: {np.rad2deg(deltae):.2f} deg")
print(f"T: {deltat*100:.2f} %")
print(f"A: {np.rad2deg(deltaa):.2f} deg")
print(f"R: {np.rad2deg(deltar):.2f} deg")

# initial state
Va = Vat
pn = 0.
pe = 0.
pd = -100.
u = xtrim.item(3)
v = xtrim.item(4)
w = xtrim.item(5)
phi = xtrim.item(6)
theta = xtrim.item(7)
psi = xtrim.item(8)
p = xtrim.item(9)
q = xtrim.item(10)
r = xtrim.item(11)
state0 = np.array([[pn], [pe], [pd], [u], [v], [w], [phi], [theta], [psi], [p], [q], [r]])
uav.state = np.ndarray.copy(state0)
state = uav.state

# main simulation loop
print("Press Q to exit...")
while sim_time < SIM.end_time:
    t_next_plot = sim_time + SIM.ts_plotting
    while sim_time < t_next_plot:
        # update everything
        Va, alpha, beta = wind.windout(state, Va, sim_time)
        fx, fy, fz = aero.forces(state, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
        l, m, n = aero.moments(state, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
        state = uav.update(fx, fy, fz, l, m, n)
        # update anim
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = state.flatten()
        anim.update(verts, pn, pe, pd, phi, theta, psi, obj, faces)
        # increment time
        sim_time += SIM.ts_simulation
    
    # write data
    simtimes = np.concatenate((simtimes[1:], [sim_time]))
    defsa = np.concatenate((defsa[1:], [np.degrees(deltaa)]))
    defse = np.concatenate((defse[1:], [np.degrees(deltae)]))
    defsr = np.concatenate((defsr[1:], [np.degrees(deltar)]))
    forces0 = np.concatenate((forces0[1:], [fx]))
    forces1 = np.concatenate((forces1[1:], [fy]))
    forces2 = np.concatenate((forces2[1:], [fz]))
    moments0 = np.concatenate((moments0[1:], [l]))
    moments1 = np.concatenate((moments1[1:], [m]))
    moments2 = np.concatenate((moments2[1:], [n]))
    poss0 = np.concatenate((poss0[1:], [pn]))
    poss1 = np.concatenate((poss1[1:], [pe]))
    poss2 = np.concatenate((poss2[1:], [-pd]))
    vels0 = np.concatenate((vels0[1:], [u]))
    vels1 = np.concatenate((vels1[1:], [v]))
    vels2 = np.concatenate((vels2[1:], [w]))
    angs0 = np.concatenate((angs0[1:], [np.degrees(phi)]))
    angs1 = np.concatenate((angs1[1:], [np.degrees(theta)]))
    angs2 = np.concatenate((angs2[1:], [np.degrees(psi)]))
    rates0 = np.concatenate((rates0[1:], [np.degrees(p)]))
    rates1 = np.concatenate((rates1[1:], [np.degrees(q)]))
    rates2 = np.concatenate((rates2[1:], [np.degrees(r)]))
    
    # plot plots
    x0 = simtimes[0]
    x1 = simtimes[-1]
    throttle.clear(); throttle.bar(0, deltat); throttle.set_ylim(0, 1); throttle.set_xticklabels("")
    throttle.set_xticks([]); throttle.set_title("Throttle")
    deflections.clear(); deflections.plot(simtimes, defsa, label="a"); deflections.plot(simtimes, defse, label="e")
    deflections.plot(simtimes, defsr, label="r"); deflections.set_ylabel("Deflection (deg)")
    deflections.legend(loc="upper right"); deflections.grid(); deflections.set_xlim((x0, x1)); deflections.set_ylim((-100, 100))
    force.clear(); force.plot(simtimes, forces0, label="fx"); force.plot(simtimes, forces1, label="fy")
    force.plot(simtimes, forces2, label="fz"); force.set_ylabel("Force (N)"); force.legend(loc="upper right"); force.grid()
    force.set_xlim((x0, x1)); force.set_ylim((-1000, 1000))
    moment.clear(); moment.plot(simtimes, moments0, label="l"); moment.plot(simtimes, moments1, label="m")
    moment.plot(simtimes, moments2, label="n"); moment.set_ylabel("Moment (Nm)"); moment.legend(loc="upper right")
    moment.grid(); moment.set_xlim((x0, x1)); moment.set_ylim((-40, 40))
    posp.clear(); posp.plot(simtimes, poss0, label="n"); posp.plot(simtimes, poss1, label="e")
    posp.plot(simtimes, poss2, label="h"); posp.set_ylabel("Position (m)"); posp.legend(loc="upper right"); posp.grid()
    posp.set_xlim((x0, x1))
    velp.clear(); velp.plot(simtimes, vels0, label="u"); velp.plot(simtimes, vels1, label="v")
    velp.plot(simtimes, vels2, label="w"); velp.set_ylabel("Velocity (m/s)"); velp.legend(loc="upper right"); velp.grid()
    velp.set_xlim((x0, x1)); velp.set_ylim((-25, 100))
    angp.clear(); angp.plot(simtimes, angs0, label="$\phi$"); angp.plot(simtimes, angs1, label="$\\theta$")
    angp.plot(simtimes, angs2, label="$\psi$"); angp.set_ylabel("Angle (deg)"); angp.legend(loc="upper right"); angp.grid()
    angp.set_xlim((x0, x1)); angp.set_ylim((-180, 180))
    ratep.clear(); ratep.plot(simtimes, rates0, label="p"); ratep.plot(simtimes, rates1, label="q")
    ratep.plot(simtimes, rates2, label="r"); ratep.set_ylabel("Rates (deg/s)"); ratep.legend(loc="upper right"); ratep.grid()
    ratep.set_xlim((x0, x1)); ratep.set_ylim((-360, 360))
    
    # pause baby
    plt.pause(0.1)
    if keyboard.is_pressed("q"): break
    
# wait to close
plt.waitforbuttonpress()