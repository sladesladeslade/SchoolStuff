import numpy as np 
from numpy import sin as s
from numpy import cos as c
from numpy import tan as t
import UAVparams as P

class mavAero():
    def __init__(self):
        # Initial state conditions
        self.state = P.states0
        # Setting parameters
        self.S_wing = P.S_wing       
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
        self.mass = P.mass
        self.g = P.gravity
        self.k_T_p = P.k_t_p
        self.k_omega = P.k_omega

    # Blending Function
    def sigma(self, alpha):
        C1 = np.exp(-self.M*(alpha - self.alpha0))
        C2 = np.exp(self.M*(alpha - self.alpha0))
        return (1 + C1 + C2)/((1 + C1)*(1 + C2))
    
    # Define all of the different constants
    def CDalpha(self, alpha):
        # return self.C_D_p + ((self.C_L_0 + self.C_L_alpha*alpha)**2)/(np.pi*self.e*self.AR)
        return (self.C_D_0 + self.C_D_alpha*alpha)
    
    def CLaplpha(self, alpha):
        # return (1 - self.sigma(alpha))*(self.C_L_0 + self.C_L_alpha*alpha) + \
        #     self.sigma(alpha)*(2*np.sign(alpha)*(s(alpha)**2)*c(alpha))
        return (self.C_L_0 + self.C_L_alpha*alpha)

    def Cx(self, alpha):
        return -self.CDalpha(alpha)*c(alpha) + self.CLaplpha(alpha)*s(alpha)
    
    def Cxq(self, alpha):
        return -self.C_D_q*c(alpha) + self.C_L_q*s(alpha)
    
    def Cxdeltae(self, alpha):
        return -self.C_D_delta_e*c(alpha) + self.C_L_delta_e*s(alpha)
    
    def Cz(self, alpha):
        return -self.CDalpha(alpha)*s(alpha) - self.CLaplpha(alpha)*c(alpha)

    def Czq(self, alpha):
        return -self.C_D_q*s(alpha) - self.C_L_q*c(alpha)
    
    def Czdeltae(self, alpha):
        return -self.C_D_delta_e*s(alpha) - self.C_L_delta_e*c(alpha)
    
    # Define external forces
    def forces(self, state, alpha, beta, deltaa, deltae, deltar, deltat, Va):
        # Grab states
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = state.flatten()

        # Create matrices for gravity, aerodynamics, and propulsion
        grav = np.array([[-self.mass*self.g*s(theta)],
                         [self.mass*self.g*c(theta)*s(phi)],
                         [self.mass*self.g*c(theta)*c(phi)]])
        aero = 0.5*self.rho*(Va**2)*self.S_wing*np.array([[self.Cx(alpha) + self.Cxq(alpha)*(self.c/(2*Va))*q + self.Cxdeltae(alpha)*deltae],
                                                          [self.C_Y_0 + self.C_Y_beta*beta + self.C_Y_p*(self.b/(2*Va))*p + \
                                                                self.C_Y_r*(self.b/(2*Va)*r) + self.C_Y_delta_a*deltaa + self.C_Y_delta_r*deltar],
                                                          [self.Cz(alpha) + self.Czq(alpha)*(self.c/(2*Va))*q + self.Czdeltae(alpha)*deltae]])
        prop = 0.5*self.rho*self.S_prop*self.C_prop*np.array([[(self.k_motor*deltat)**2 - Va**2],
                                                              [0],
                                                              [0]])
    
        # Create forces matrix
        forces = grav + aero + prop
        if forces.shape == (3, 3, 1):     # some weird thing is happening with array shapes at the start so this is my fix
            forces = forces[0]
            forces[0][0] = forces[0][0][0]
        else:
            forces = forces
        fx, fy, fz = forces.flatten()

        return fx, fy, fz

    def moments(self, state, alpha, beta, deltaa, deltae, deltar, deltat, Va):
        # Grab states
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = state.flatten()

        # Create matrices for aerodynamics and propulsion
        aero = 0.5*self.rho*(Va**2)*self.S_wing*np.array([[self.b*(self.C_ell_0 + self.C_ell_beta*beta + self.C_ell_p*(self.b/(2*Va))*p + \
                                                                  self.C_ell_r*(self.b/(2*Va))*r + self.C_ell_delta_a*deltaa + self.C_ell_delta_r*deltar)],
                                                          [self.c*(self.C_m_0 + self.C_m_alpha*alpha + self.C_m_q*(self.c/(2*Va))*q + self.C_m_delta_e*deltae)],
                                                          [self.b*(self.C_n_0 + self.C_n_beta*beta + self.C_n_p*(self.b/(2*Va))*p + \
                                                                  self.C_n_r*(self.b/(2*Va))*r + self.C_n_delta_a*deltaa + self.C_n_delta_r*deltar)]])
        prop = np.array([[-self.k_T_p*(self.k_omega*deltat)**2],
                         [0],
                         [0]])
        
        # Create moments matrix 
        moments = aero + prop
        if moments.shape == (3, 3, 1):    # some weird thing is happening with array shapes at the start so this is my fix
            moments[0][0] = moments[0][0][0]
            moments = moments[0]
        else:
            moments = moments
        l, m, n = moments.flatten()

        return l, m, n