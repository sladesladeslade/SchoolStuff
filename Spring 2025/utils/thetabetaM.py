# Theta-Beta-M Relation Function
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
from scipy.optimize import fsolve


def shockRelation(theta:float=None, M:float=None, B:float=None, y=1.4):
    """
    Function to calculate theta/beta/mach depending on knowns. Requires two of the 3 inputs (+ gamma).
    
    Parameters
    ----------
    theta : float, optional
        Wedge angle (deg).
    M : float, optional
        Mach number.
    B : float, optional
        Mach wave angle (deg).
    y : float=1.4
        Specific heat ratio gamma.
        
    Returns
    -------
    Missing of the 3 parameters based on inputs (deg).
    """
    if M != None and B != None:
        B = np.deg2rad(B)
        theta = np.atan(2/np.tan(B)*((M**2)*(np.sin(B))**2 - 1)/((M**2)*(y + np.cos(2*B)) + 2))
        return np.rad2deg(theta)
    elif M != None and theta != None:
        def f(B):
            B = np.deg2rad(B)
            return np.rad2deg(np.atan(2/np.tan(B)*((M**2)*(np.sin(B))**2 - 1)/((M**2)*(y + np.cos(2*B)) + 2))) - theta
        B = fsolve(f, theta + 0.1)[0]
        return B
    elif B != None and theta != None:
        B = np.deg2rad(B)
        def f(M):
            return np.rad2deg(np.atan(2/np.tan(B)*((M**2)*(np.sin(B))**2 - 1)/((M**2)*(y + np.cos(2*B)) + 2))) - theta
        M = fsolve(f, 1.1)[0]
        return M