mr = 0.25       # rotor mass (kg)
mc = 1.         # center mass (kg)
Jc = 0.0042     # jesus (kgm^2)
d = 0.3         # arm diameter (m)
u = 0.1         # somethin (kg/s)
g = 9.81        # dont change (m/s^2)

z0 = 3.         # initial z position (m)
zdot0 = 0.      # initial zdot (m/s)
h0 = 0.5
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