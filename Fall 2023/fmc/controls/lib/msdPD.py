import mssimparams as P
import numpy as np


class controller:
    def __init__(self, Fmax):
        self.m = P.m1
        self.k = P.k
        self.b = P.b
        self.kP = 4.5
        self.kD = 12
        self.Fmax = Fmax

    def update(self, zc, state):
        z = state[0]
        zdot = state[1]
        feq = P.k*z
        fc = self.kP*(zc - z)  - self.kD*zdot
        f = feq + fc
        if abs(f[0]) > self.Fmax:
            f[0] = self.Fmax*np.sign(f[0])
        return f[0]