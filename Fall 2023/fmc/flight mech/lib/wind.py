# Wind Generation Class
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
from rotations import Rvb
rng = np.random.default_rng()


class wind():
    """
    TODO
    """
    
    def __init__(self, Vs, Vgust):
        self.Vs = Vs
        self.Vgust = Vgust
        
    
    def gust(self):
        uwg = rng.random()*self.Vgust[0][0]
        vwg = rng.random()*self.Vgust[1][0]
        wwg = rng.random()*self.Vgust[2][0]
        return np.array([[uwg],[vwg],[wwg]])
        
    
    def windout(self, states):
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        gusts = self.gust()
        windtot = Rvb(phi, theta, psi)*self.Vs + gusts
        Var = np.array([[u - windtot[0][0]],
                        [v - windtot[1][0]],
                        [w - windtot[2][0]]])
        ur, vr, wr = Var.flatten()
        Va = np.sqrt(ur**2 + vr**2 + wr**2)
        alpha = np.arctan(wr/ur)
        beta = np.arcsin(vr/Va)
        return Va, alpha, beta