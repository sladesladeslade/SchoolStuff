# VTOL Dynamics Class
# Slade Brooks
# brooksl@mail.uc.edu
# bruh 2

import numpy as np
import VTOLsimparams2 as P


class VTOLDynamics():
    """
    TODO
    """
    
    def __init__(self):
        self.state = np.array([[P.z0], [P.h0], [P.theta0], [P.zdot0], [P.hdot0], [P.thetadot0]])
        self.Ts = P.Ts
        self.mc = P.mc
        self.mr = P.mr
        self.Jc = P.Jc
        self.u = P.u
        self.d = P.d
        self.g = P.g
        self.force_limit = P.F_max
    
    
    def update(self, fr, fl):
        """
        TODO
        """
        self.rk4_step(fr, fl)
        y = self.h()
        return y
    
    
    def f(self, state, fr, fl):
        """
        TODO
        """
        # get states
        z = state[0][0]
        zdot = state[3][0]
        h = state[1][0]
        hdot = state[4][0]
        theta = state[2][0]
        thetadot = state[5][0]
        
        # do eoms
        zddot = (-1*(fr + fl)*np.sin(theta) - self.u*zdot)/(self.mc + 2*self.mr)
        hddot = ((fr + fl)*np.cos(theta))/(self.mc + 2*self.mr) - self.g
        thetaddot = (self.d*(fr - fl))/(self.Jc + 2*self.mr*self.d**2)
        
        # build xdot
        xdot = np.array([[zdot], [hdot], [thetadot], [zddot], [hddot], [thetaddot]])
        return xdot
    
    
    def h(self):
        """
        TODO
        """
        z = self.state[0][0]
        h = self.state[1][0]
        theta = self.state[2][0]
        y = np.array([[z], [h], [theta]])
        return y
    
    
    def rk4_step(self, fr, fl):
        """
        TODO
        """
        F1 = self.f(self.state, fr, fl)
        F2 = self.f(self.state + self.Ts/2*F1, fr, fl)
        F3 = self.f(self.state + self.Ts/2*F2, fr, fl)
        F4 = self.f(self.state + self.Ts*F3, fr, fl)
        self.state += self.Ts/6*(F1 + 2*F2 + 2*F3 + F4)