import bbsimparams2 as P


class controller:
    def __init__(self):
        self.m1 = P.m1
        self.m2 = P.m2
        self.l = P.l
        self.g = P.g
        self.kPt = 1.825
        self.kDt = 1.173
        self.kPz = -0.004939
        self.kDz = -0.031743

    def update(self, zc, state):
        z, theta, zdot, thetadot = state.flatten()
        feq = self.m2*self.g/2 + self.m1*self.g*z/self.l
        zs = self.kPz*(zc - z) - self.kDz*zdot
        fc = self.kPt*(zs - theta)  - self.kDt*thetadot
        f = feq + fc
        if f > 15:
            f = 15
        return f