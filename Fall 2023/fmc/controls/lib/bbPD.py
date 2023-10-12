import bbsimparams as P


class controller:
    def __init__(self):
        self.kP = 4.5
        self.kD = 12

    def update(self, zc, state):
        z = state[0]
        zdot = state[2]
        feq = 13.23
        fc = self.kP*(zc - z)  - self.kD*zdot
        f = feq + fc
        return f[0]