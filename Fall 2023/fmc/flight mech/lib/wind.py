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
    
    def __init__(self, Vs):
        self.Vs = Vs
        
    
    def drydenGust(self, Va):
        s = rng.random()
        Lu = Lv = 200
        Lw = 50
        ou = ov = 1.06
        ow = 0.7
        uwg = ou*np.sqrt(2*Va/Lu)*(1/(s + Va/Lu))
        vwg = ov*np.sqrt(3*Va/Lv)*((s + Va/(np.sqrt(3)*Lv))/((s + Va/Lv)**2))
        wwg = ow*np.sqrt(3*Va/Lw)*((s + Va/(np.sqrt(3)*Lw))/((s + Va/Lw)**2))
        return np.array([[uwg],[vwg],[wwg]])
        
    
    def windout(self, states, Va):
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        gusts = self.drydenGust(Va)
        windtot = Rvb(phi, theta, psi)*self.Vs + gusts
        Var = np.array([[u - windtot[0][0]],
                        [v - windtot[1][0]],
                        [w - windtot[2][0]]])
        ur, vr, wr = Var.flatten()
        Va = np.sqrt(ur**2 + vr**2 + wr**2)
        alpha = np.arctan(wr/ur)
        beta = np.arcsin(vr/Va)
        return Va, alpha, beta