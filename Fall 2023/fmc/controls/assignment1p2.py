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


# define system
vtol = dynamics.VTOLDynamics()
anim = animation.VTOLAnim(limits=1, flag=False)

# sim loop
t = P.t_start
while t < P.t_end:
    # dod ynamics
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        fr = 1.5*9.81
        fl = 1.5*9.81
        y = vtol.update(fr, fl)
        t += P.Ts
        
    # update anim
    anim.update(vtol.state)
    plt.pause(0.001)
    
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()