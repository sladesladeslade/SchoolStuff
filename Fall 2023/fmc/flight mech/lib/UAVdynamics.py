# UAV Dynamics Class
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import UAVparams as P
from numpy import sin as s
from numpy import cos as c
from numpy import tan as t


class UAVdynamics():
    """
    TODO
    """
    
    def __init__(self):
        self.state = P.states0
        self.Ts = P.Ts
        self.mass = P.m
        self.jxz = P.jxz
        self.jx = P.jx
        self.jy = P.jy
        self.jz = P.jz
        self.gam = self.jx*self.jz - self.jxz**2
        self.gam1 = (self.jxz*(self.jx - self.jy + self.jz))/self.gam
        self.gam2 = (self.jz*(self.jz - self.jy) + self.jxz**2)/self.gam
        self.gam3 = self.jz/self.gam
        self.gam4 = self.jxz/self.gam
        self.gam5 = (self.jz - self.jx)/self.jy
        self.gam6 = self.jxz/self.jy
        self.gam7 = (self.jx*(self.jx - self.jy) + self.jxz**2)
        self.gam8 = self.jx/self.gam
        
    
    def update(self, fx, fy, fz, l, m, n):
        """
        TODO
        """
        self.rk4_step(fx, fy, fz, l, m, n)
        y = self.h()
        return y
    
    
    def f(self, state, fx, fy, fz, l, m, n):
        """
        TODO
        """
        # get states
        pn = state[0][0]
        pe = state[1][0]
        pd = state[2][0]
        u = state[3][0]
        v = state[4][0]
        w = state[5][0]
        phi = state[6][0]
        theta = state[7][0]
        psi = state[8][0]
        p = state[9][0]
        q = state[10][0]
        r = state[11][0]
        
        # do EOMS
        eom1 = np.array([[c(theta)*c(psi), s(phi)*s(theta)*c(phi) - c(phi)*s(psi), c(phi)*s(theta)*c(psi) + s(phi)*s(psi)],
                         [c(theta)*s(psi), s(phi)*s(theta)*s(psi) + c(phi)*c(psi), c(phi)*s(theta)*s(psi) - s(phi)*c(psi)],
                         [-s(theta), s(phi)*s(theta), c(phi)*c(theta)]])
        eom1x = np.array([[u], [v], [w]])
        eom1s = np.dot(eom1, eom1x)
        eom2s = np.array([[r*v - q*w], [p*w - r*u], [q*u - p*v]]) + (1/self.mass)*np.array([[fx], [fy], [fz]])
        eom3 = np.array([[1, s(phi)*t(theta), c(phi)*t(theta)],
                         [0, c(phi), -s(phi)],
                         [0, s(phi)/c(theta), c(phi)/c(theta)]])
        eom3x = np.array([[p], [q], [r]])
        eom3s = np.dot(eom3, eom3x)
        eom4s = np.array([[self.gam1*p*q - self.gam2*q*r], [self.gam5*p*r - self.gam6*(p**2 - r**2)],
            [self.gam7*p*q - self.gam1*q*r]]) + np.array([[self.gam3*l + self.gam4*n], [m/self.jy],
            [self.gam4*l + self.gam8*n]])
            
        # build and return output
        xdot = np.array([[eom1s[0][0]], [eom1s[1][0]], [eom1s[2][0]], [eom2s[0][0]], [eom2s[1][0]], [eom2s[2][0]],
                         [eom3s[0][0]], [eom3s[1][0]], [eom3s[2][0]], [eom4s[0][0]], [eom4s[1][0]], [eom4s[2][0]]])
        return xdot
    
    
    def h(self):
        """
        TODO
        """
        pn = self.state[0][0]
        pe = self.state[1][0]
        pd = self.state[2][0]
        u = self.state[3][0]
        v = self.state[4][0]
        w = self.state[5][0]
        phi = self.state[6][0]
        theta = self.state[7][0]
        psi = self.state[8][0]
        p = self.state[9][0]
        q = self.state[10][0]
        r = self.state[11][0]
        y = np.array([[pn], [pe], [pd], [u], [v], [w], [phi], [theta], [psi], [p], [q], [r]])
        return y
    
    
    def rk4_step(self, fx, fy, fz, l, m, n):
        """
        TODO
        """
        F1 = self.f(self.state, fx, fy, fz, l, m, n)
        F2 = self.f(self.state + self.Ts/2*F1, fx, fy, fz, l, m, n)
        F3 = self.f(self.state + self.Ts/2*F2, fx, fy, fz, l, m, n)
        F4 = self.f(self.state + self.Ts*F3, fx, fy, fz, l, m, n)
        self.state += self.Ts/6*(F1 + 2*F2 + 2*F3 + F4)