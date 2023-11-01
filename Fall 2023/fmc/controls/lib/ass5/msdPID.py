import mssimparams as P
import numpy as np
import mssimparams as P


class controller:
    def __init__(self, Fmax, sigma, flag):
        self.m = P.m1
        self.k = P.k
        self.b = P.b
        self.kp = 7.2
        self.kd = 3.05
        self.ki = 0.4
        self.limit = Fmax
        self.sigma = sigma
        self.flag = flag
        self.beta = (2.0*sigma-P.Ts)/(2.0*sigma+P.Ts)
        self.Ts = P.Ts # sample rate
        self.y_dot = 0.0 # estimated derivative of y
        self.y_d1 = 0.0 # Signal y delayed by one sample
        self.error_dot = 0.0 # estimated derivative of error
        self.error_d1 = 0.0 # Error delayed by one sample
        self.integrator = 0.0 # integrator
    
    
    def update(self, zr, h):
        z = h
        f_eq = P.k * z
        f_s = self.PID(zr, z)
        f = f_eq + f_s
        f = self.saturate(f)
        return f
    
    
    def PID(self, y_r, y):
        # Compute the current error
        error = y_r - y
        # integrate error using trapazoidal rule
        self.integrator = self.integrator \
        + (self.Ts/2) * (error + self.error_d1)
        # PID Control
        if self.flag is True:
            # differentiate error
            self.error_dot = self.beta * self.error_dot \
                + (1-self.beta)/self.Ts * (error - self.error_d1)
            # PID control
            u_unsat = self.kp*error \
                + self.ki*self.integrator \
                + self.kd*self.error_dot
        else:
            # differentiate y
            self.y_dot = self.beta * self.y_dot \
            + (1-self.beta)/self.Ts * (y - self.y_d1)
            # PID control
            u_unsat = self.kp*error \
                + self.ki*self.integrator \
                - self.kd*self.y_dot
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.ki != 0.0:
            self.integrator = self.integrator \
                + 1.0 / self.ki * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1 = error
        self.y_d1 = y
        return u_sat
    
    def saturate(self,u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u