import numpy as np

m1 = np.random.uniform(0.8*5, 1.2*5)         # mass of mass (kg)
k = np.random.uniform(0.8*3, 1.2*3)          # spring constant (N/m)
b = np.random.uniform(0.8*0.5, 1.2*0.5)      # damping coefficient (Ns/m)

z0 = 0.         # initial z position (m)
zdot0 = 0.      # initial zdot (m/s)

t_start = 0.    # sim start time
t_end = 60.     # sim end time
Ts = 0.1       # sim sample time
t_plot = 0.1    # plot/anim update rate

F_max = 5.      # max force (N)

w = 0.5         # width of block (m)