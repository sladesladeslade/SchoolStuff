# Ball-Beam Dynamics Class
# Slade Brooks
# brooksl@mail.uc.edu
# bruh 2

import numpy as np
import bbsimparams as P


class ballBeamDynamics():
    """
    Mass spring damper system dynamics class.
    
    Methods
    -------
    # TODO
    """
    
    def __init__(self):
        self.state = np.array([[P.z0], [P.theta0],
                               [P.zdot0], [P.thetadot0]])       # initial states
        self.Ts = P.Ts      # sim time step
        self.m1 = P.m1      # mass of ball
        self.m2 = P.m2      # mass of beam
        self.l = P.l        # length of beam
        self.g = P.g        # gravity
        self.r = P.r        # radius of ball
        
    
    def update(self, u):
        """"""
        self.rk4_step(u)
        y = self.h()
        return y
    
    
    def f(self, state, u):
        """"""
        # get states
        z, theta, zdot, thetadot = state.flatten()
        F = u
        
        # get zddot from EOM
        zddot = z*thetadot**2 - self.g*np.sin(theta)
        
        # get thetaddot from EOM
        thetaddot = (1/((self.m2*self.l**2)/3 + self.m1*z**2))*(F*self.l*np.cos(theta) - 2*self.m1*z*zdot*thetadot \
                    - self.m1*self.g*z*np.cos(theta) - ((self.m2*self.g*self.l)/2)*np.cos(theta))
        
        # build xdot
        xdot = np.array([[zdot], [thetadot], [zddot], [thetaddot]])
        return xdot
    
    
    def h(self):
        """"""
        z, theta, _, _ = self.state.flatten()
        y = np.array([[z], [theta]])
        return y
    
    
    def rk4_step(self, u):
        """"""
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts/2*F1, u)
        F3 = self.f(self.state + self.Ts/2*F2, u)
        F4 = self.f(self.state + self.Ts*F3, u)
        self.state += self.Ts/6*(F1 + 2*F2 + 2*F3 + F4)