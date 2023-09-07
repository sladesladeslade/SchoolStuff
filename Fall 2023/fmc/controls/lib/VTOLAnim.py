# VTOL Animation Class
# Slade Brooks
# brooksl@mail.uc.edu
# pls work part 2

import matplotlib.pyplot as plt
import matplotlib.patches as mpat
import VTOLsimparams as P
import numpy as np


class VTOLAnim():
    """
    TODO
    """
    
    def __init__(self, limits, flag=False):
        
        # set up plot
        self.fig = plt.figure(1)
        if flag == True:
            self.ax = self.fig.add_subplot(1, 2, 1)
        else:
            self.ax = self.fig.add_subplot(1, 1, 1)
        
        # draw ground
        plt.plot([-limits, limits], [0, 0], "k-", linewidth=3)
        
        # for lsit of objects
        self.handle = []
        
        # init flag
        self.flag_init = True
        
        # uh
        self.limits = limits
        
    
    def update(self, state):
        """
        TODO
        """
        # update state
        z = state[0][0]
        h = state[1][0]
        theta = np.rad2deg(state[2][0])
        
        # draw elements
        self.draw_vtol(z, h, theta)
        
        # plot formatting
        self.ax.set_aspect("equal")
        self.ax.set_ylim(top=self.limits)
        self.ax.set_xlim(left=-self.limits, right=self.limits)
        
        # check for first time
        if self.flag_init == True:
            self.flag_init = False
            
    
    def draw_vtol(self, z, h, theta):
        """
        TODO
        """
        # corner of rectangle
        x = z - P.w/2
        y = h - P.h/2
        corner = (x, y)
        
        # location of rotors
        xr = z + (P.w/2 + P.d*np.cos(-theta))
        yr = h - P.d*np.sin(-theta)
        rr = (xr, yr)
        xl = z - (P.w/2 + P.d*np.cos(-theta))
        yl = h + P.d*np.sin(-theta)
        lr = (xl, yl)
        
        # draw body
        if self.flag_init == True:
            self.handle.append(mpat.Rectangle(corner, P.w, P.h, fc="black", ec="black"))
            self.handle.append(mpat.CirclePolygon(rr, 0.1, resolution=15, fc="black", ec="black"))
            self.handle.append(mpat.CirclePolygon(lr, 0.1, resolution=15, fc="black", ec="black"))
            self.ax.add_patch(self.handle[0])
            self.ax.add_patch(self.handle[1])
            self.ax.add_patch(self.handle[2])
        else:
            self.handle[0].set_xy(corner)
            self.handle[1].xy = rr
            self.handle[2].xy = lr