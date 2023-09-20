# UAV Aerodynamics Class
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import UAVparams as P


class UAVaero():
    """
    TODO
    """
    
    def __init__(self):
        self.S = P.S_wing
        self.b = P.b
        self.c = P.c
        self.S_prop = P.S_prop
        self.rho = P.rho
        self.e = P.e
        self.AR = P.AR
        self.C_L_0 = P.C_L_0
        self.C_D_0 = P.C_D_0
        self.C_m_0 = P.C_m_0
        self.C_L_alpha = P.C_L_alpha
        self.C_D_alpha = P.C_D_alpha
        self.C_m_alpha = P.C_m_alpha
        self.C_L_q = P.C_L_q
        self.C_D_q = P.C_D_q
        self.C_m_q = P.C_m_q
        self.C_L_delta_e = P.C_L_delta_e
        self.C_D_delta_e = P.C_D_delta_e
        self.C_m_delta_e = P.C_m_delta_e
        self.M = P.M
        self.alpha0 = P.alpha0
        self.epsilon = P.epsilon
        self.C_D_p = P.C_D_p
        self.C_Y_0 = P.C_Y_0
        self.C_ell_0 = P.C_ell_0
        self.C_n_0 = P.C_n_0
        self.C_Y_beta = P.C_Y_beta
        self.C_ell_beta = P.C_ell_beta
        self.C_n_beta = P.C_n_beta
        self.C_Y_p = P.C_Y_p
        self.C_ell_p = P.C_ell_p
        self.C_n_p = P.C_n_p
        self.C_Y_r = P.C_Y_r
        self.C_ell_r = P.C_ell_r
        self.C_n_r = P.C_n_r
        self.C_Y_delta_a = P.C_Y_delta_a
        self.C_ell_delta_a = P.C_ell_delta_a
        self.C_n_delta_a = P.C_n_delta_a
        self.C_Y_delta_r = P.C_Y_delta_r
        self.C_ell_delta_r = P.C_ell_delta_r
        self.C_n_delta_r = P.C_n_delta_r
        self.C_prop = P.C_prop
        self.k_motor = P.k_motor
        self.m = P.mass
        self.g = P.g
        
    
    def sigma(self, alpha):
        C1 = np.exp(-self.M*(alpha - self.alpha0))
        C2 = np.exp(self.M*(alpha - self.alpha0))
        return (1 + C1 + C2)/((1 + C1)*(1 + C2))
    
    
    def CDalpha(self, alpha): return self.C_D_p + ((self.C_L_0 + self.C_L_alpha*alpha)**2)/(np.pi*self.e*self.AR)
    
    
    def CLalpha(self, alpha): return (1 - self.sigma(alpha))*(self.C_L_0 + self.C_L_alpha*alpha) +\
            self.sigma(alpha)*(2*np.sign(alpha)*(np.sin(alpha)**2)*np.cos(alpha))
    
    
    def Cx(self, alpha): return -self.CDalpha(alpha)*np.cos(alpha) + self.CLalpha(alpha)*np.sin(alpha)
    
    
    def Cxq(self, alpha): return -self.C_D_q*np.cos(alpha) + self.C_L_q*np.sin(alpha)
    
    
    def Cxdeltae(self, alpha): return -self.C_D_delta_e*np.cos(alpha) + self.C_L_delta_e*np.sin(alpha)
    
    
    def Cz(self, alpha): return -self.CDalpha(alpha)*np.sin(alpha) - self.CLalpha(alpha)*np.cos(alpha)
    
    
    def Czq(self, alpha): return -self.C_D_q*np.sin(alpha) - self.C_L_q*np.cos(alpha)
    
    
    def Czdeltae(self, alpha): return -self.C_D_delta_e*np.sin(alpha) - self.C_L_delta_e*np.cos(alpha)
    
    
    def forces(self, states, alpha, beta, deltaa, deltae, deltar, deltat, Va):
        """
        """
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        grav = np.array([[-self.m*self.g*np.sin(theta)],
                       [self.m*self.g*np.cos(theta)*np.sin(phi)],
                       [self.m*self.g*np.cos(theta)*np.cos(phi)]])
        aero = 0.5*self.rho*(Va**2)*self.S*np.array([[self.Cx(alpha) + self.Cxq(alpha)*(self.c/(2*Va))*q +\
                                                self.Cxdeltae(alpha)*deltae],
                                                   [self.C_Y_0 + self.C_Y_beta*beta + self.C_Y_p*(self.b/(2*Va))*p +\
                                                self.C_Y_r*(self.b/(2*Va))*r + self.C_Y_delta_a*deltaa +\
                                                self.C_Y_delta_r*deltar],
                                                   [self.Cz(alpha) + self.Czq(alpha)*(self.c/(2*Va))*q +\
                                                self.Czdeltae(alpha)*deltae]])
        prop = 0.5*self.rho*self.S_prop*self.C_prop*np.array([[(self.k_motor*deltat)**2 - Va**2], [0], [0]])
        forces = grav + aero + prop
        return forces