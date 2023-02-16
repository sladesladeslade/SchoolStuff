# Slade Brooks
# brooksl@mail.uc.edu
# 02.07.2023

# stream func class


import math
import numpy as np


def uniform(y2, V):
    """
    Calculates the stream function (psi) of a uniform flow
    
    :param y2: y pos of point
    :param V: free stream velocity
    :returns: strength of V @ pt (x2, y2)
    """

    V = V*y2

    return V


def sourcesink(x1, y1, x2, y2, A):
    """
    Calculates the stream function (psi) of a source/sink
    
    :param x1: x pos of source
    :param y1: y pos of source
    :param x2: x pos of point
    :param y2: y pos of point
    :param A: lambda (strength of source/sink)
    :returns: strength of V @ pt (x2, y2)
    """

    # calculate angle
    theta = np.arctan2((y2 - y1),(x2 - x1))

    # calculate strength
    V = (A/(2*(math.pi)))*theta

    return V


def vortex(x1, y1, x2, y2, gamma):
    """
    Calculates the stream function (psi) of a vortex
    
    :param x1: x pos of source
    :param y1: y pos of source
    :param x2: x pos of point
    :param y2: y pos of point
    :param gamma: constant of vortex
    :returns: strength of V @ pt (x2, y2)
    """

    # calculate radius
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    # calculate strength
    V = (gamma/(2*math.pi))*np.log(r)

    return V


def doublet(x1, y1, x2, y2, kappa):
    """
    Calculates the stream function (psi) of a doublet
    
    :param x1: x pos of source
    :param y1: y pos of source
    :param x2: x pos of point
    :param y2: y pos of point
    :param kappa: constant of doublet
    :returns: strength of V @ pt (x2, y2)
    """

    # calculate radius
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    # calculate angle
    theta = np.arctan2((y2 - y1),(x2 - x1))

    # calculate V
    V = (-kappa/(2*math.pi))*(np.sin(theta)/r)

    return V


if __name__ == "__main__":
    print("test")