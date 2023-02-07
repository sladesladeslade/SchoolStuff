# Slade Brooks
# brooksl@mail.uc.edu
# 02.07.2023

# stream func class


import math
import numpy as np


class streamFunc():
    """
    Collection of functions to calculate
    """

    def sourcesinkPsi(self, x1, y1, x2, y2, A):
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


if __name__ == "__main__":
    test = streamFunc()
    #print(test.sourcesinkPSI(1, 1, 2, 2, 0))       # tested angle