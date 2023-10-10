# Ball-Beam Animation Class
# Slade Brooks
# brooksl@mail.uc.edu
# 2 ez

import matplotlib.pyplot as plt
import matplotlib.patches as mpat
import bbsimparams as P
import numpy as np


class ballBeamAnim():
    """"""
    
    def __init__(self, limits, flag=False):
        
        # set up plot
        self.fig = plt.figure(1)
        if flag == True:
            self.ax = self.fig.add_subplot(1, 2, 1)
        else:
            self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax.set_xlim(left=-0.2, right=limits)

        # for list of objects
        self.handle = []
        
        # init flag
        self.flag_init = True
        
        # uh
        self.limits = limits
        
    
    def update(self, state):
        """"""
        # update state
        z = state[0][0]
        theta = state[1][0]
        
        # draw elements
        self.drawBallBeam(z, theta)
        
        # plot fromatting
        self.ax.set_aspect("equal")
        self.ax.set_xlim(left=-0.2, right=self.limits)
        self.ax.set_ylim(bottom=-0.5, top=0.5)
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")
        
        # check for first time runthrough
        if self.flag_init == True:
            self.flag_init = False
            
    
    def drawBallBeam(self, z, theta):
        """"""
        # find position of ball
        c = (z*np.cos(theta) + P.r*np.cos(theta + np.pi/2), z*np.sin(theta) + P.r*np.sin(theta + np.pi/2))
        
        # draw ball
        if self.flag_init == True:
            self.handle.append(mpat.CirclePolygon(c, P.r, resolution=15, fc="blue", ec="blue"))
            self.handle.append(mpat.Rectangle((0, 0), P.l, -0.01, theta, fc="white", ec="black"))
            self.ax.add_patch(self.handle[0])
            self.ax.add_patch(self.handle[1])
        else:
            self.handle[0].xy = c
            self.handle[1].set_angle(np.rad2deg(theta))