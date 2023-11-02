import numpy.random as r

mr = 0.25       # rotor mass (kg)
mc = r.uniform(0.8*1, 1.2*1)         # center mass (kg)
Jc = r.uniform(0.8*0.0042, 1.2*0.0042)     # jesus (kgm^2)
d = r.uniform(0.8*0.3, 1.2*0.3)         # arm diameter (m)
u = r.uniform(0.8*0.1, 1.2*0.1)         # somethin (kg/s)
g = 9.81        # dont change (m/s^2)

z0 = 0.         # initial z position (m)
zdot0 = 0.      # initial zdot (m/s)
h0 = 0.25
hdot0 = 0.
theta0 = 0.
thetadot0 = 0.

t_start = 0.    # sim start time
t_end = 80.     # sim end time
Ts = 0.1       # sim sample time
t_plot = 0.1    # plot/anim update rate

F_max = 50.      # max force (N)

w = 0.5         # width of block (m)
h = 0.5        # height of block (m)