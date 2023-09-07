import sys
sys.path.append("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\controls\\examples\\")
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np 
import parameters.pendulumParam as P
# if you are having difficulty with the graphics, 
# try using one of the following backends  
# See https://matplotlib.org/stable/users/explain/backends.html
# import matplotlib
# matplotlib.use('qtagg')  # requires pyqt or pyside
# matplotlib.use('ipympl')  # requires ipympl
# matplotlib.use('gtk3agg')  # requires pyGObject and pycairo
# matplotlib.use('gtk4agg')  # requires pyGObject and pycairo
# matplotlib.use('gtk3cairo')  # requires pyGObject and pycairo
# matplotlib.use('gtk4cairo')  # requires pyGObject and pycairo
# matplotlib.use('tkagg')  # requires TkInter
# matplotlib.use('wxagg')  # requires wxPython


class pendulumAnimation:
    def __init__(self):
        self.flag_init = True  # Used to indicate initialization
        # Initialize a figure and axes object
        self.fig, self.ax = plt.subplots()
        # Initializes a list of objects (patches and lines)
        self.handle = []
        # Specify the x,y axis limits
        plt.axis([-3*P.ell, 3*P.ell, -0.1, 3*P.ell])
        # Draw line for the ground
        plt.plot([-2*P.ell, 2*P.ell], [0, 0], 'b--')
        # label axes
        plt.xlabel('z')

    def update(self, state):
        z = state[0][0]  # Horizontal position of cart, m
        theta = state[1][0]  # Angle of pendulum, rads
        # draw plot elements: cart, bob, rod
        self.draw_cart(z)
        self.draw_bob(z, theta)
        self.draw_rod(z, theta)
        self.ax.axis('equal')
        # Set initialization flag to False after first call
        if self.flag_init == True:
            self.flag_init = False

    def draw_cart(self, z):
        # specify bottom left corner of rectangle
        x = z-P.w/2.0
        y = P.gap
        corner = (x, y)
        # create rectangle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the Rectangle patch and append its handle
            # to the handle list
            self.handle.append(
                mpatches.Rectangle(corner, P.w, P.h, fc='blue', ec='black'))
            # Add the patch to the axes
            self.ax.add_patch(self.handle[0])
        else:
            self.handle[0].set_xy(corner)  # Update patch

    def draw_bob(self, z, theta):
        # specify center of circle
        x = z+(P.ell+P.radius)*np.sin(theta)
        y = P.gap+P.h+(P.ell+P.radius)*np.cos(theta)
        center = (x, y)
        # create circle on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the CirclePolygon patch and append its handle
            # to the handle list
            self.handle.append(
                mpatches.CirclePolygon(center, radius=P.radius,
                                       resolution=15, fc='limegreen', ec='black'))
            # Add the patch to the axes
            self.ax.add_patch(self.handle[1])
        else:
            self.handle[1].xy = center

    def draw_rod(self, z, theta):
        # specify x-y points of the rod
        X = [z, z+P.ell*np.sin(theta)]
        Y = [P.gap+P.h, P.gap+P.h+P.ell*np.cos(theta)]
        # create rod on first call, update on subsequent calls
        if self.flag_init is True:
            # Create the line object and append its handle
            # to the handle list.
            line, = self.ax.plot(X, Y, lw=1, c='black')
            self.handle.append(line)
        else:
            self.handle[2].set_xdata(X)
            self.handle[2].set_ydata(Y)