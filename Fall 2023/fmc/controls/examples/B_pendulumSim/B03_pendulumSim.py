import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\examples\\")
import matplotlib.pyplot as plt
import parameters.pendulumParam as P
from tools.signalGenerator import signalGenerator
from viewer.pendulumAnimation import pendulumAnimation
from viewer.dataPlotter import dataPlotter
from dynamics.pendulumDynamics import pendulumDynamics

# instantiate pendulum, controller, and reference classes
pendulum = pendulumDynamics(alpha=0.0)
reference = signalGenerator(amplitude=0.5, frequency=0.02)
force = signalGenerator(amplitude=1, frequency=1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = pendulumAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        u = 0.5*force.sin(t)
        y = pendulum.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts
    # update animation and data plots at rate t_plot
    animation.update(pendulum.state)
    dataPlot.update(t, r, pendulum.state, u)
    plt.pause(0.001)  # allows time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
