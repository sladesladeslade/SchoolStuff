import mssimparams as P
import numpy as np
import VTOLsimparams as P


class controller:
    def __init__(self, Fmax, sigma, flag):
        self.mr = P.mr
        self.mc = P.mc
        self.Jc = P.Jc
        self.u = P.u
        self.d = P.d
        self.g = P.g
        self.kpt = 0.3721
        self.kdt = 0.1913
        self.kpz = -0.0077095
        self.kdz = -0.032858
        self.kph = 0.11345
        self.kdh = 0.5835
        self.kih = 0.
        self.kit = 0.
        self.kiz = 0.
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
    
    
    def update(self, hc, zc, h, z, theta):
        # do theta part
        taueq = 0.
        zs = self.PIDz(zc, z)
        tauc = self.PIDt(zs, theta)
        tau = taueq + tauc
        
        # do height part
        feq = self.g*(self.mc + 2*self.mr)
        fc = self.PIDh(hc, h)
        f = feq + fc
        
        # get indiv fs from answers
        fr = (tau + self.d*f)/(2*self.d)
        fl = f - fr
        
        # saturate them bois
        fr = self.saturate(fr)
        fl = self.saturate(fl)
        
        return fr, fl
    
    
    def PIDz(self, y_r, y):
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
            u_unsat = self.kpz*error \
                + self.kiz*self.integrator \
                + self.kdz*self.error_dot
        else:
            # differentiate y
            self.y_dot = self.beta * self.y_dot \
            + (1-self.beta)/self.Ts * (y - self.y_d1)
            # PID control
            u_unsat = self.kpz*error \
                + self.kiz*self.integrator \
                - self.kdz*self.y_dot
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.kiz != 0.0:
            self.integrator = self.integrator \
                + 1.0 / self.kiz * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1 = error
        self.y_d1 = y
        return u_sat
    
    
    def PIDt(self, y_r, y):
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
            u_unsat = self.kpt*error \
                + self.kit*self.integrator \
                + self.kdt*self.error_dot
        else:
            # differentiate y
            self.y_dot = self.beta * self.y_dot \
            + (1-self.beta)/self.Ts * (y - self.y_d1)
            # PID control
            u_unsat = self.kpt*error \
                + self.kit*self.integrator \
                - self.kdt*self.y_dot
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.kit != 0.0:
            self.integrator = self.integrator \
                + 1.0 / self.kit * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1 = error
        self.y_d1 = y
        return u_sat
    
    
    def PIDh(self, y_r, y):
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
            u_unsat = self.kph*error \
                + self.kih*self.integrator \
                + self.kdh*self.error_dot
        else:
            # differentiate y
            self.y_dot = self.beta * self.y_dot \
            + (1-self.beta)/self.Ts * (y - self.y_d1)
            # PID control
            u_unsat = self.kph*error \
                + self.kih*self.integrator \
                - self.kdh*self.y_dot
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.kih != 0.0:
            self.integrator = self.integrator \
                + 1.0 / self.kih * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1 = error
        self.y_d1 = y
        return u_sat
    
    
    def saturate(self,u):
        if u > self.limit:
            u = self.limit
        elif u < 0:
            u = 0
        return u