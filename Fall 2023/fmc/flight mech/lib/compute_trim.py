import numpy as np
# load message types
#from message_types.msg_state import MsgState
import lib.UAVparams as P
# from library.rotations import Quaternion2Euler, Quaternion2Rotation, Euler2Rotation
from numpy import cos, sin, tan
# import control
# from control.matlab import *
from lib.UAVdynamics import UAVdynamics
from lib.UAVlinaero import UAVaero
from scipy.optimize import minimize


class ComputeTrim:

    def __init__(self):
        self.P = P
        self.Ts = P.Ts
        self.forces_mom = UAVaero()
        self.mav = UAVdynamics()
        
    def compute_trim(self, Va, Y, R):
        x0 = np.array([0, 0, 0])
        res = minimize(lambda x: self.compute_trim_cost(x, Va, Y, R), x0, method='nelder-mead', options={'xatol': 1e-8, 'disp': False})
        x_trim, u_trim = self.compute_trim_states_input(res.x, Va, Y, R)
        return x_trim, u_trim

    def compute_trim_states_input(self, x, Va, Y, R):
        # Inertial parameters
        jx = self.P.jx
        jy = self.P.jy
        jz = self.P.jz
        jxz = self.P.jxz
        gravity = self.P.g
        mass = self.P.m
        gamma = jx*jz - jxz**2
        gamma1 = (jxz*(jx - jy + jz))/gamma
        gamma2 = (jz*(jz - jy) + jxz**2)/gamma
        gamma3 = jz/gamma
        gamma4 = jxz/gamma
        gamma5 = (jz - jx)/jy
        gamma6 = jxz/jy
        gamma7 = ((jx - jy)*jx + jxz**2)/gamma
        gamma8 = jx/gamma
        #Va0=self.P.Va0

        ## aerodynamic parameters
        S_wing        = self.P.S
        b             = self.P.b
        c             = self.P.c
        S_prop        = self.P.S_prop
        rho           = self.P.rho
        e             = self.P.e
        AR            = self.P.AR
        C_L_0         = self.P.C_L_0
        C_D_0         = self.P.C_D_0
        C_m_0         = self.P.C_m_0
        C_L_alpha     = self.P.C_L_alpha
        C_D_alpha     = self.P.C_D_alpha
        C_m_alpha     = self.P.C_m_alpha
        C_L_q         = self.P.C_L_q
        C_D_q         = self.P.C_D_q
        C_m_q         = self.P.C_m_q
        C_L_delta_e   = self.P.C_L_delta_e
        C_D_delta_e   = self.P.C_D_delta_e
        C_m_delta_e   = self.P.C_m_delta_e
        M             = self.P.M 
        alpha0        = self.P.alpha0
        epsilon       = self.P.epsilon
        C_D_p         = self.P.C_D_p
        C_Y_0         = self.P.C_Y_0
        C_ell_0       = self.P.C_ell_0
        C_n_0         = self.P.C_n_0
        C_Y_beta      = self.P.C_Y_beta
        C_ell_beta    = self.P.C_ell_beta
        C_n_beta      = self.P.C_n_beta 
        C_Y_p         = self.P.C_Y_p
        C_ell_p       = self.P.C_ell_p
        C_n_p         = self.P.C_n_p
        C_Y_r         = self.P.C_Y_r
        C_ell_r       = self.P.C_ell_r
        C_n_r         = self.P.C_n_r
        C_Y_delta_a   = self.P.C_Y_delta_a
        C_ell_delta_a = self.P.C_ell_delta_a
        C_n_delta_a   = self.P.C_n_delta_a
        C_Y_delta_r   = self.P.C_Y_delta_r
        C_ell_delta_r = self.P.C_ell_delta_r
        C_n_delta_r   = self.P.C_n_delta_r
        C_prop        = self.P.C_prop
        k_motor       = self.P.k_motor

        alpha = x[0]
        beta = x[1]
        phi = x[2]
        u = Va*cos(alpha)*cos(beta)
        v = Va*sin(beta)
        w = Va*sin(alpha)*cos(beta)
        theta = alpha + Y
        p = (-Va/R)*sin(theta)
        q = (Va/R)*sin(phi)*cos(theta)
        r = (Va/R)*cos(phi)*cos(theta)
        
        x_trim = np.array([[0.],
                           [0.],
                           [0.],
                           [u],
                           [v],
                           [w],
                           [phi],
                           [theta],
                           [0.],
                           [p],
                           [q],
                           [r]], dtype=float)

        C_L = C_L_0 + C_L_alpha*alpha
        C_D = C_D_0 + C_D_alpha*alpha
        
        C_X = -C_D*cos(alpha) + C_L*sin(alpha)
        C_X_q = -C_D_q*cos(alpha) + C_L_q*sin(alpha)
        C_X_delta_e = -C_D_delta_e*cos(alpha) + C_L_delta_e*sin(alpha)
        
        C_Z = -C_D*sin(alpha) - C_L*cos(alpha)
        C_Z_q = -C_D_q*sin(alpha) - C_L_q*cos(alpha)
        C_Z_delta_e = -C_D_delta_e*sin(alpha) - C_L_delta_e*cos(alpha)
        d_e = (((jxz*(p**2 - r**2) + (jx - jz)*p*r)/(0.5*rho*(Va**2)*c*S_wing)) - C_m_0 - C_m_alpha*alpha - C_m_q*((c*q)/(2*Va)))/C_m_delta_e
        d_t = np.sqrt(((2*mass*(-r*v + q*w + gravity*sin(theta)) - rho*(Va**2)*S_wing*(C_X + C_X_q*((c*q)/(2*Va)) + \
                        C_X_delta_e*d_e))/(rho*S_prop*C_prop*k_motor**2)) + ((Va**2)/(k_motor**2)))
        temp_1 = np.linalg.inv(np.array([[C_ell_delta_a, C_ell_delta_r],
                                         [C_n_delta_a, C_n_delta_r]]))
        temp_2 = np.array([[((-gamma1*p*q+gamma2*q*r)/(0.5*rho*(Va**2)*S_wing*b)) - C_ell_0 - C_ell_beta*beta - C_ell_p*((b*p)/(2*Va)) - C_ell_r*((b*r)/(2*Va))],
                        [((-gamma7*p*q+gamma1*q*r)/(0.5*rho*(Va**2)*S_wing*b)) - C_n_0 - C_n_beta*beta - C_n_p*((b*p)/(2*Va)) - C_n_r*((b*r)/(2*Va))]])
        
        temp_3=np.matmul(temp_1, temp_2)
        
        d_a = temp_3[0][0]
        d_r = temp_3[1][0]
        
        u_trim = np.array([[d_e],
                           [d_t],
                           [d_a],
                           [d_r]], dtype=float)

        return (x_trim, u_trim)
    
    def compute_trim_cost(self, x, Va, Y, R):

        alpha = x[0]
        beta = x[1]

        #compute X_dot_star
        x_dot_star = np.array([[0.],
                        [0.],
                        [-Va*sin(Y)], # I am using Pd_dot not hdot..that is why there is a sign change
                        [0.],
                        [0.],
                        [0.],
                        [0.],
                        [0.],
                        [Va/R],
                        [0.],
                        [0.],
                        [0.]], dtype=float)
        
        #compute trimmed states
        x_trim, u_trim = self.compute_trim_states_input(x, Va, Y, R)
        d_e = u_trim[0][0]
        d_t = u_trim[1][0]
        d_a = u_trim[2][0]
        d_r = u_trim[3][0]
        
        fx, fy, fz = self.forces_mom.forces(x_trim, alpha, beta, d_a, d_e, d_r, d_t, Va).flatten()
        l, m, n = self.forces_mom.moments(x_trim, alpha, beta, d_a, d_e, d_r, d_t, Va).flatten()
        
        states_dot = self.mav.f(x_trim, fx, fy, fz, l, m, n) # 

        J = np.linalg.norm(x_dot_star - states_dot)**2

        return J