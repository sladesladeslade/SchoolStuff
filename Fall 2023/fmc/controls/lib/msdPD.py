import mssimparams as P


class controller:
    def __init__(self):
        self.m = P.m1
        self.k = P.k
        self.b = P.b
        self.kP = 4.5
        self.kD = 12

    def update(self, zc, state):
        z = state[0]
        zdot = state[1]
        feq = P.k*z
        fc = self.kP*(zc - z)  - self.kD*zdot
        f = feq + fc
        return f[0]