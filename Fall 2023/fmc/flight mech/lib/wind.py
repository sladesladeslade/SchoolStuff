# Wind Generation Class
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
from rotations import Rvb
import control.matlab as mat
import simparams as SIM


class wind():
    """
    TODO
    """
    
    def __init__(self, Vs):
        self.Vs = Vs
        
    
    def drydenGust(self, Va, t):
        # get white noise vals
        su = np.random.normal(0, 1, 1)[0]
        sv = np.random.normal(0, 1, 1)[0]
        sw = np.random.normal(0, 1, 1)[0]
        
        # set dryden gust model params
        Lu = Lv = 200
        Lw = 50
        ou = ov = 1.06
        ow = 0.7

        # make coeffs for tfs
        C1 = ou*np.sqrt(2*Va/Lu)
        C2 = ov*np.sqrt(3*Va/Lv)
        C3 = ow*np.sqrt(3*Va/Lw)
        
        # set up tfs
        Hu = mat.tf([0, C1],[1, Va/Lu])
        Hv = mat.tf([C2, C2*Va/(np.sqrt(3)*Lv)],[1, 2*Va/Lv, (Va/Lv)**2])
        Hw = mat.tf([C3, C3*Va/(np.sqrt(3)*Lw)],[1, 2*Va/Lw, (Va/Lw)**2])
        
        # calc gusts from tfs
        T = [0, t]
        uwg, _, _ = mat.lsim(Hu, su, T)
        vwg, _, _ = mat.lsim(Hv, sv, T)
        wwg, _, _ = mat.lsim(Hw, sw, T)
        return np.array([[uwg[1]],[vwg[1]],[wwg[1]]])
        
    
    def windout(self, states, Va, simtime):
        pn, pe, pd, u, v, w, phi, theta, psi, p, q, r = states.flatten()
        gusts = self.drydenGust(Va, simtime)
        windtot = np.matmul(Rvb(phi, theta, psi), self.Vs) + gusts
        Var = np.array([[u - windtot[0][0]],
                        [v - windtot[1][0]],
                        [w - windtot[2][0]]])
        ur, vr, wr = Var.flatten()
        Va = np.sqrt(ur**2 + vr**2 + wr**2)
        alpha = np.arctan(wr/ur)
        beta = np.arcsin(vr/Va)
        return Va, alpha, beta