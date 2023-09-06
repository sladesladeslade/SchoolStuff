# Mass-Spring Dynamics Class
# Slade Brooks
# brooksl@mail.uc.edu
# bruh

import numpy as np
import simparams as P


class massSpringDynamics():
    """
    Mass spring damper system dynamics class.
    
    Methods
    -------
    # TODO
    """
    
    def __init__(self):
        self.state = np.array([P.z0, P.zdot0])      # initial states
        self.Ts = P.Ts                              # sim time step
        self.m1 = P.m1                              # mass of mass
        self.k = P.k                                # spring coef
        self.b = P.b                                # damp coef
        self.g = P.g                                # gravity
        self.force_limit = P.F_max                  # max force

    
    def update(self, u):
        """
        Updates output of system based on force.
        
        Parameters
        ----------
        u : float
            Force.
            
        Returns
        -------
        y : np.ndarray
            Output array.
        """
        u = self.saturate(u, self.force_limit)
        self.rk4_step(u)
        y = self.h()
        return y
    
    
    def f(self, state, u):
        """
        Method to create state space xdot from EOM.
        Ax + Bu
        
        Parameters
        ----------
        state : np.ndarray
            Array of current state.
        u : float
            Force.
            
        Returns
        -------
        xdot : np.ndarray
            State space xdot.
        """
        # get states
        z = state[0]
        zdot = state[1]
        F = u
        
        # get zddot from EOM
        zddot = (-self.b/self.m1)*zdot - (k/self.m1)*z + F/self.m1
        
        # build xdot
        xdot = np.array([[zdot], [zddot]])
        return xdot
    
    
    def h(self):
        """
        Method to create state space output from state.
        y = h(x)
        
        Returns
        -------
        y : np.ndarray
            State space output array.
        """
        z = self.state[0]
        y = np.array([z])
        return y
    
    
    def rk4_step(self, u):
        """
        Integrate ODE w/ Runge-Kutta Rk4.
        
        Parameters
        ----------
        u : np.ndarray
            Force.
        """
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts/2*F1, u)
        F3 = self.f(self.state + self.Ts/2*F2, u)
        F4 = self.f(self.state + self.Ts*F3, u)
        self.state += self.Ts/6*(F1 + 2*F2 + 2*F3 + F4)

    
    @staticmethod
    def saturate(u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u