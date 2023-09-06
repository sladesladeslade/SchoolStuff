# Mass-Spring Animation Class
# Slade Brooks
# brooksl@mail.uc.edu
# pls work

import matplotlib.pyplot as plt
import matplotlib.patches as mpat
import numpy as np
import simparams as P


class massSpringAnim():
    """
    # TODO
    """
    
    def __init__(self, limits):
        
        # set up plot
        self.fig, self.ax = plt.subplots()
        plt.axis([-limits, limits, -0.1, 2])
        # draw ground
        plt.plot([-limits, limits], [0, 0], "k-")
            
        # for list of objects
        self.handle = []
        
        # init flag
        self.flag_init = True
        
    
    def update(self, state):
        """
        Update animation based on state.
        TODO
        """
        # update state
        z = state[0]
        
        # draw elements
        self.draw_mass(z)
        # self.draw_springdamp(z)
        self.ax.axis("equal")
        
        # check for first time runthrough
        if self.flag_init == True:
            self.flag_init = False

    
    def draw_mass(self, z):
        """
        TODO
        """
        # corner of rectangle
        x = z - P.w/2
        corner = (x, 0.01)
        
        # draw square
        if self.flag_init == True:
            self.handle.append(mpat.Rectangle(corner, P.w, P.w, fc="blue", ec="black"))
            self.ax.add_patch(self.handle[0])
        else:
            self.handle[0].set_xy(corner)
            
    
    def draw_springdamp(self, z):
        """
        TODO
        """
        self.ax.plot([0, z], [0.5, 0.5], "r-")
        self.ax.plot([0, z], [0.25, 0.25], "b-")