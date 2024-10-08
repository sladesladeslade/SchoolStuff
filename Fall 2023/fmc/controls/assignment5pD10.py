# AEEM4042 Controls Assignment #5 Part D8A
# Slade Brooks
# brooksl@mail.uc.edu

import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\lib\\")
import matplotlib.pyplot as plt
import lib.ass5.mssimparams as P
import lib.massSpringAnim as animation
import lib.ass5.massSpringDynamics as dynamics
import keyboard
import lib.ass5.msdPID as ctr
import lib.signalGenerator as sig


# define system
msd = dynamics.massSpringDynamics()
anim = animation.massSpringAnim(limits=2, flag=True)
ctr = ctr.controller(6., 0.05, False)
sig = sig.signalGenerator(1, 0.05)
print(f"M: {P.m1:.2f}")
print(f"K: {P.k:.2f}")
print(f"B: {P.b:.2f}")

# set controller
ctr.kd = 7.2
ctr.kp = 3.05
ctr.ki = 0.4

# add subplots
zes = anim.fig.add_subplot(2, 2, 2)
fes = anim.fig.add_subplot(2, 2, 4)
zes.set_ylabel("z (m)")
fes.set_ylabel("Force (N)")

# initials
z = 0
u = 0
target = 1
zs = [z]
fs = [u]
targets = [target]

# simulation loop
t = P.t_start
simtimes = [t]
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        target = sig.square(t)
        u = ctr.update(target, msd.h()[0])
        z = msd.update(u)
        t += P.Ts
    
    # update daterp
    simtimes.append(t)
    zs.append(z[0])
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
    
    plt.pause(0.1)
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()