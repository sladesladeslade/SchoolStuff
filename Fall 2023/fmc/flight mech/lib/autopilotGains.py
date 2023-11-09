# Compute Autopilot Gains Class
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import UAVparams as P


# set xtrim and utrim
x_trim = np.array([0., 0., 0., 34.9957342, -0.178955066,  0.516297664, -0.0123004054, 0.0147520897, 0., 0., 0., 0.])
u_trim = np.array([-0.027192652340385412, 0.46236388117478694, -0.0038422387096093266, -0.004796900558456178])

# do the tf stuff to get the answers and woo
Va_trim = np.sqrt(x_trim[3]**2 + x_trim[4]**2 + x_trim[5]**2)
alpha_trim = np.arctan(x_trim[5]/x_trim[3])
beta_trim = np.arctan(x_trim[4]/Va_trim)
theta_trim = x_trim[7]
gamma = P.jx * P.jz - P.jxz ** 2
gamma3 = P.jz / gamma
gamma4 = P.jxz / gamma
Cpp = gamma3 * P.C_ell_p + gamma4 * P.C_n_p
Cpdelta_a = gamma3 * P.C_ell_delta_a + gamma4 * P.C_n_delta_a
a_phi1 = -.5 * P.rho * Va_trim**2 * P.S_wing * P.b * Cpp * (P.b / (2 * Va_trim))
a_phi2 = .5 * P.rho * Va_trim**2 * P.S_wing * P.b * Cpdelta_a
a_theta1 = -(P.rho * Va_trim**2 * P.c * P.S_wing)/(2 * P.jy) * P.C_m_q * (P.c/(2 * Va_trim))
a_theta2 = -(P.rho * Va_trim**2 * P.c * P.S_wing)/(2 * P.jy) * P.C_m_alpha
a_theta3 = (P.rho * Va_trim**2 * P.c * P.S_wing)/(2 * P.jy) * P.C_m_delta_e
a_V1 = ((P.rho * Va_trim * P.S_wing)/P.M) * (P.C_D_0 + (P.C_D_alpha * alpha_trim) \
        + (P.C_D_delta_e * u_trim[0])) + (P.rho * P.S_prop)/(P.M) * P.C_prop * Va_trim
a_V2 = (P.rho * P.S_prop)/(P.M) * P.C_prop * P.k_motor**2 * u_trim[3]
a_V3 = P.gravity * np.cos(theta_trim - alpha_trim)

# zeta
zeta = .707

# roll
tr_roll = 0.5
wn_roll = 2.2/tr_roll
kp_roll = (wn_roll**2)/a_phi2
kd_roll = (2*zeta*wn_roll - a_phi1)/a_phi2
ki_roll = 0

# course hold
tr_course = 5.
wn_course = 2.2/tr_course
kp_course = (2*zeta*wn_course*Va_trim)/P.gravity
kd_course = 0
ki_course = ((wn_course**2)*Va_trim)/P.gravity

# pitch attitude hold
zeta_pitch = 0.1
tr_pitch = 0.1
wn_pitch = 2.2/tr_pitch
kp_pitch = ((wn_pitch**2) - a_theta2)/a_theta3
kd_pitch = (2*zeta_pitch*wn_pitch - a_theta1)/a_theta3
ki_pitch = 0
ktheta_DC = (kp_pitch*a_theta3)/(a_theta2 + kp_pitch*a_theta3)

# altitude from pitch gain
tr_altitude = 1.
wn_altitude = 2.2/tr_altitude
kp_altitude = (2*zeta*wn_altitude)/(ktheta_DC*Va_trim)
kd_altitude = 0
ki_altitude = (wn_altitude**2)/(ktheta_DC*Va_trim)

# airspeed from pitch
tr_airspeed = 0.01
wn_airspeed = 2.2/tr_airspeed
kp_airspeed = (a_V1 - 2*zeta*wn_airspeed)/ktheta_DC
kd_airspeed = 0
ki_airspeed = (wn_airspeed**2)/(ktheta_DC*P.gravity)

# airspeed from throttle
tr_throttle = 0.1
wn_throttle = 2.2/tr_throttle
kp_throttle = (2*zeta*wn_throttle - a_V1)/a_V2
kd_throttle = 0
ki_throttle = (wn_throttle**2)/a_V2