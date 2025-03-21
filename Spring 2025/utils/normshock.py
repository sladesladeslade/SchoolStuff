# Normal Shock Functions
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np


def normshock(Mx:float, y:float):
    """
    Function to calculate change in flow parameters across normal shock.
    
    Parameters
    ----------
    Mx : float
        Upstream Mach number.
    y : float
        Specific heat ratio.
        
    Returns
    -------
    My : float
        Downstream Mach number.
    Pty_Ptx : float
        Total pressure ratio (= sonic area ratio).
    Py_Px : float
        Static pressure ratio.
    rhoy_rhox : float
        Density ratio (= velocity ratio).
    Ty_Tx : float
        Static temperature ratio.
    """
    # calculate My
    My = np.sqrt((Mx**2 + 2/(y - 1))/(2*y/(y - 1)*Mx**2 - 1))
    
    # calculate total pressure ratio
    Pty_Ptx = ((((y + 1)/2*Mx**2)/(1 + (y - 1)/2*Mx**2))**(y/(y - 1)))*((2*y/(y + 1)*Mx**2 - (y - 1)/(y + 1))**(-1/(y - 1)))
    
    # calculate static pressure ratio
    Py_Px = 2*y/(y + 1)*Mx**2 - (y - 1)/(y + 1)
    
    # calculate static temperature ratio
    Ty_Tx = (2*y*Mx**2 - (y - 1))*((y - 1)*Mx**2 + 2)/(((y + 1)**2)*(Mx**2))
    
    # calculate density ratio
    rhoy_rhox = Py_Px/Ty_Tx
    
    return My, Pty_Ptx, Py_Px, rhoy_rhox, Ty_Tx