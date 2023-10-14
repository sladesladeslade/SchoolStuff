m1 = 0.35           # ball mass (kg)
m2 = 2.             # beam mass (kg)
l = 0.5             # beam length (m)
g = 9.8             # gravity (m/s^2)
r = 0.05            # radius of ball (m)

z0 = 0.             # initial position of ball on beam (m)
zdot0 = 0.          # initial speed of ball on beam (m/s)
theta0 = 0.         # initial angle of beam (rad)
thetadot0 = 0.      # initial theta dot of beam (rad/s)

t_start = 0.        # sim start time
t_end = 600.        # sim end time
Ts = 0.001          # sim sample time
t_plot = 0.1        # plot/anim update rate