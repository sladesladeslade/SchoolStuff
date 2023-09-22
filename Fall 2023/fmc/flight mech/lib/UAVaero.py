# UAV Aerodynamics Class
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
from UAVparams import *


class UAVaero():
    """
    TODO
    """   
    
    @staticmethod
    def sigma(alpha):
        C1 = np.exp(-M*(alpha - alpha0))
        C2 = np.exp(M*(alpha - alpha0))
        return (1 + C1 + C2)/((1 + C1)*(1 + C2))
    
    
    @staticmethod
    def CDalpha(alpha): return C_D_p + ((C_L_0 + C_L_alpha*alpha)**2)/(np.pi*e*AR)
    
    
    def CLalpha(self, alpha): return (1 - self.sigma(alpha))*(C_L_0 + C_L_alpha*alpha) +\
            self.sigma(alpha)*(2*np.sign(alpha)*(np.sin(alpha)**2)*np.cos(alpha))
    
    
    def Cx(self, alpha): return -self.CDalpha(alpha)*np.cos(alpha) + self.CLalpha(alpha)*np.sin(alpha)
    
    
    @staticmethod
    def Cxq(alpha): return -C_D_q*np.cos(alpha) + C_L_q*np.sin(alpha)
    
    
    @staticmethod
    def Cxdeltae(alpha): return -C_D_delta_e*np.cos(alpha) + C_L_delta_e*np.sin(alpha)
    
    
    def Cz(self, alpha): return -self.CDalpha(alpha)*np.sin(alpha) - self.CLalpha(alpha)*np.cos(alpha)
    
    
    @staticmethod
    def Czq(alpha): return -C_D_q*np.sin(alpha) - C_L_q*np.cos(alpha)
    
    
    @staticmethod
    def Czdeltae(alpha): return -C_D_delta_e*np.sin(alpha) - C_L_delta_e*np.cos(alpha)
    
    
    def forces(self, states, alpha, beta, deltaa, deltae, deltar, deltat, Va):
        """
        """
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        grav = np.array([[-m*g*np.sin(theta)],
                       [m*g*np.cos(theta)*np.sin(phi)],
                       [m*g*np.cos(theta)*np.cos(phi)]])
        aero = 0.5*rho*(Va**2)*S*np.array([[self.Cx(alpha) + self.Cxq(alpha)*(c/(2*Va))*q + self.Cxdeltae(alpha)*deltae],
                                           [C_Y_0 + C_Y_beta*beta + C_Y_p*(b/(2*Va))*p + C_Y_r*(b/(2*Va))*r +\
                                            C_Y_delta_a*deltaa + C_Y_delta_r*deltar],
                                           [self.Cz(alpha) + self.Czq(alpha)*(c/(2*Va))*q + self.Czdeltae(alpha)*deltae]])
        prop = 0.5*rho*S_prop*C_prop*np.array([[(k_motor*deltat)**2 - Va**2], [0], [0]])
        forces = grav + aero + prop
        return forces
    
    
    def moments(self, states, alpha, beta, deltaa, deltae, deltar, deltat, Va):
        """
        """
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        moments = 0.5*rho*(Va**2)*S*np.array([[b*(C_ell_0 + C_ell_beta*beta + C_ell_p*(b/(2*Va))*p +\
                                                C_ell_r*(b/(2*Va))*r + C_ell_delta_a*deltaa + C_ell_delta_r*deltar)],
                                              [c*(C_m_0 + C_m_alpha*alpha + C_m_q*(c/(2*Va))*q + C_m_delta_e*deltae)],
                                              [b*(C_n_0 + C_n_beta*beta + C_n_p*(b/(2*Va))*p + C_n_r*(b/(2*Va))*r +\
                                                C_n_delta_a*deltaa + C_n_delta_r*deltar)]]) +\
                                            np.array([[-k_t_p*(k_omega*deltat)**2],[0],[0]])
        return moments