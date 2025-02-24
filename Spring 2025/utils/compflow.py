# Compressible Flow Functions
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np


def compflow(M:float, y:float):
    """
    Function to calculate compressible flow parameters.
    
    Parameters
    ----------
    M : float
        Flow Mach number.
    y : float
        Specific heat ratio.
        
    Returns
    -------
    T_Tt : float
        Static to total temperature ratio.
    P_Pt : float
        Static to total pressure ratio.
    rho_rhot : float
        Static to total density ratio.
    A_Astar : float
        Sonic area ratio.
    MFPsRgc : float
        Mass flow parameter * sqrt(R/gc).
    mu : float
        Mach wave angle (deg).
    v : float
        P-M angle (deg).
        
    Notes
    -----
    sqrt(R/gc) = 1.28758 lbf-s/lbm-sqrt(R) or 16.9115 N-s/kg-sqrt(K)
    """
    # calculate temperature ratio
    T_Tt = (1 + (y - 1)/2*M**2)**(-1)
    
    # calculate pressure ratio
    P_Pt = (1 + (y - 1)/2*M**2)**(-y/(y - 1))
    
    # calculate density ratio
    rho_rhot = (1 + (y - 1)/2*M**2)**(-1/(y - 1))

    # calculate sonic area ratio
    A_Astar = 1/M*(2/(y + 1)*(1 + (y - 1)/2*M**2))**((y + 1)/(2*(y - 1)))
    
    # calculate gross MFP*sqrt(R/gc)
    MFPsRgc = M*np.sqrt(y)*(1 + (y - 1)/2*M**2)**(-(y + 1)/(2*(y - 1)))
    
    # check for supersonic
    if M >= 1:
        # calculate mach wave angle
        mu = np.rad2deg(np.asin(1/M))
        
        # calculate P-M function
        v = np.rad2deg(np.sqrt((y + 1)/(y - 1))*np.atan(np.sqrt((y - 1)/(y + 1)*(M**2 - 1))) -\
            np.atan(np.sqrt(M**2 - 1)))
    else:
        mu = 0.
        v = 0.
        
    return np.round(T_Tt, 6), np.round(P_Pt, 6), np.round(rho_rhot, 6), np.round(A_Astar, 6), np.round(MFPsRgc, 6), \
        np.round(mu, 6), np.round(v, 6)