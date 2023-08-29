
import sys
sys.path.append('.')# one directory up
import numpy as np
from math import cos, sin, tan
import scipy.linalg as linalg
import parameters.simulation_parameters as SIM
from viewers.cube_animation import cube_animation
from tools.signalGenerator import signalGenerator
from tools.sliders import sliders

state=np.array([[0], [0], [-1], [0], [0], [0], [0], [0], [0], [0], [0], [0]])
cube_anim=cube_animation(state, scale=5)
my_slider=sliders()
   
temp = signalGenerator(amplitude=0.5, frequency=0.1)

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Command-Q to exit...")
pn=state[0,0]
pe=state[1,0]
pd=state[2,0]
phi=state[6,0]
theta=state[7,0]
psi=state[8,0]
while sim_time < SIM.end_time:
    # reading from sliders
    phi=my_slider.roll_slider.val
    theta=my_slider.pitch_slider.val
    psi=my_slider.yaw_slider.val
    #phi=temp.sin(sim_time)
    #print(phi)
    cube_anim.update(pn, pe, pd, phi, theta, psi) # -pd for height


    # -------increment time-------------
    sim_time += SIM.ts_simulation