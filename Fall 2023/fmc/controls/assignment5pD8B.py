# AEEM4042 Controls Assignment #5 Part D8B
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


# define system
msd = dynamics.massSpringDynamics()
anim = animation.massSpringAnim(limits=2, flag=True)
ctr = ctr.controller(Fmax=6)

# tune controller
trise = 1.6375
wn = 2.2/trise
damping = 0.7
a1 = 2*damping*wn
a0 = wn**2
ctr.kD = P.m1*(a1 - P.b/P.m1)
ctr.kP = P.m1*(a0 - P.k/P.m1)

# add subplots
zes = anim.fig.add_subplot(222)
fes = anim.fig.add_subplot(224)
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
    else:
        target = 1
        
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