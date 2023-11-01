import VTOLsimparams2 as P


class controller:
    def __init__(self):
        self.mr = P.mr
        self.mc = P.mc
        self.d = P.d
        self.g = P.g
        self.u = P.u
        self.jc = P.Jc
        self.d = P.d
        self.kPt = 0.3721
        self.kDt = 0.1913
        self.kPz = -0.0077095
        self.kDz = -0.032858
        self.kPh = 0.11345
        self.kDh = 0.5835

    def update(self, zc, hc, state):
        z, h, theta, zdot, hdot, thetadot = state.flatten()
        
        # do theta part
        taueq = 0.
        zs = self.kPz*(zc - z) - self.kDz*zdot
        tauc = self.kPt*(zs - theta)  - self.kDt*thetadot
        tau = taueq + tauc
        
        # do height part
        feq = self.g*(self.mc + 2*self.mr)
        fc = self.kPh*(hc - h)  - self.kDh*hdot
        f = feq + fc
        
        # get indiv fs from answers
        fr = (tau + self.d*f)/(2*self.d)
        fl = f - fr

        return fr, fl