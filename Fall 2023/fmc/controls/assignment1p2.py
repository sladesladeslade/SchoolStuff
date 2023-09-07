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
anim = animation.VTOLAnim(limits=10, flag=False)

# weight
w = 1.5*9.81/2

# sim loop
t = P.t_start
while t < P.t_end:
    # dod ynamics
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        print(f"{t:.2f}")
        if t < 4:
            fr = w*1.05
            fl = fr
        elif t < 4.1:
            fr = 0
            fl = fr
        elif t < 10:
            fr = w*0.98
            fl = fr
        elif t < 10.25:
            fr = w*1.05
            fl = fr/1.05
        elif t < 10.55:
            fr = w/1.05
            fl = fr*1.05
        elif t < 17:
            fr = w*1.01
            fl = fr
        elif t < 17.25:
            fr = w*1.05
            fl = w/1.05
        else:
            fr = w
            fl = fr
            
        y = vtol.update(fr, fl)
        t += P.Ts
        
    # update anim
    anim.update(vtol.state)
    plt.pause(0.1)
    
    if keyboard.is_pressed("q"): break

# wait to close
plt.waitforbuttonpress()