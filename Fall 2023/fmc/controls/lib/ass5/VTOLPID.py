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
        self.y_dotz = 0.0 # estimated derivative of y
        self.y_d1z = 0.0 # Signal y delayed by one sample
        self.error_dotz = 0.0 # estimated derivative of error
        self.error_d1z = 0.0 # Error delayed by one sample
        self.integratorz = 0.0 # integrator
        self.y_dott = 0.0 # estimated derivative of y
        self.y_d1t = 0.0 # Signal y delayed by one sample
        self.error_dott = 0.0 # estimated derivative of error
        self.error_d1t = 0.0 # Error delayed by one sample
        self.integratort = 0.0 # integrator
        self.y_doth = 0.0 # estimated derivative of y
        self.y_d1h = 0.0 # Signal y delayed by one sample
        self.error_doth = 0.0 # estimated derivative of error
        self.error_d1h = 0.0 # Error delayed by one sample
        self.integratorh = 0.0 # integrator
    
    
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
        self.integratorz = self.integratorz \
        + (self.Ts/2) * (error + self.error_d1z)
        # PID Control
        if self.flag is True:
            # differentiate error
            self.error_dotz = self.beta * self.error_dotz \
                + (1-self.beta)/self.Ts * (error - self.error_d1z)
            # PID control
            u_unsat = self.kpz*error \
                + self.kiz*self.integratorz \
                + self.kdz*self.error_dotz
        else:
            # differentiate y
            self.y_dotz = self.beta * self.y_dotz \
            + (1-self.beta)/self.Ts * (y - self.y_d1z)
            # PID control
            u_unsat = self.kpz*error \
                + self.kiz*self.integratorz \
                - self.kdz*self.y_dotz
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.kiz != 0.0:
            self.integratorz = self.integratorz \
                + 1.0 / self.kiz * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1z = error
        self.y_d1z = y
        return u_sat
    
    
    def PIDt(self, y_r, y):
        # Compute the current error
        error = y_r - y
        # integrate error using trapazoidal rule
        self.integratort = self.integratort \
        + (self.Ts/2) * (error + self.error_d1t)
        # PID Control
        if self.flag is True:
            # differentiate error
            self.error_dott = self.beta * self.error_dott \
                + (1-self.beta)/self.Ts * (error - self.error_d1t)
            # PID control
            u_unsat = self.kpt*error \
                + self.kit*self.integratort \
                + self.kdt*self.error_dott
        else:
            # differentiate y
            self.y_dott = self.beta * self.y_dott \
            + (1-self.beta)/self.Ts * (y - self.y_d1t)
            # PID control
            u_unsat = self.kpt*error \
                + self.kit*self.integratort \
                - self.kdt*self.y_dott
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.kit != 0.0:
            self.integratort = self.integratort \
                + 1.0 / self.kit * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1t = error
        self.y_d1t = y
        return u_sat
    
    
    def PIDh(self, y_r, y):
        # Compute the current error
        error = y_r - y
        # integrate error using trapazoidal rule
        self.integratorh = self.integratorh \
        + (self.Ts/2) * (error + self.error_d1h)
        # PID Control
        if self.flag is True:
            # differentiate error
            self.error_doth = self.beta * self.error_doth \
                + (1-self.beta)/self.Ts * (error - self.error_d1h)
            # PID control
            u_unsat = self.kph*error \
                + self.kih*self.integratorh \
                + self.kdh*self.error_doth
        else:
            # differentiate y
            self.y_doth = self.beta * self.y_doth \
            + (1-self.beta)/self.Ts * (y - self.y_d1h)
            # PID control
            u_unsat = self.kph*error \
                + self.kih*self.integratorh \
                - self.kdh*self.y_doth
        # return saturated control signal
        u_sat = self.saturate(u_unsat)
        # integrator anti - windup
        if self.kih != 0.0:
            self.integratorh = self.integratorh \
                + 1.0 / self.kih * (u_sat - u_unsat)
        # update delayed variables
        self.error_d1h = error
        self.y_d1h = y
        return u_sat
    
    
    def saturate(self,u):
        if u > self.limit:
            u = self.limit
        elif u < 0:
            u = 0
        return u