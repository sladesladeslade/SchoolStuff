# ATAVP Project 2
# Slade Brooks
# brooksl@mail.uc.edu

import sys, os
sys.path.append(os.getcwd() + r"/Spring 2025/")
from utils.normshock import normshock
from utils.compflow import compflow
from utils.thetabetaM import shockRelation
from utils.stdatmos import stdAtmos
from scipy.optimize import fsolve
import utils.units as uu
import numpy as np
import matplotlib.pyplot as plt
y = 1.4
std = stdAtmos()


def Ax_Ay(Mx, My, Ptx_Pty):
    Arat = 1/Ptx_Pty*My/Mx*((1 + (y-1)/2*Mx**2)/(1 + (y-1)/2*My**2))**((y + 1)/(2*(y - 1)))
    return Arat


def inlet(A2, M0sub, M0sup, M2f, AR, LD, flag):
    """
    Parameters
    ----------
    A2 : float
        Fan face area.
    M0sub : float
        Subsonic design Mach number.
    M0sup : float
        Supersonic design Mach number.
    M2 : float
        Fan face Mach number.
    AR : float
        Ramp aspect ratio (h/w).
    LD : float
        Inlet length ratio (L/D2).
    """
    # calculate the dimensions of A2
    r2 = np.sqrt(A2/np.pi)
    
    # determine subsonic area ratio A1/A2
    A1_A2 = Ax_Ay(M0sub, M2f, 1.02)
    A1 = A1_A2*A2
    
    # get dimensions of A1
    h1 = np.sqrt(AR*A1)
    hr = r2 - h1/2.
    
    # calculate total length and diffuser length
    L = LD*2*r2
    L12 = (r2 - h1/2.)/np.tan(np.deg2rad(7.))
    L01 = L - L12
    
    # solve for required beta1 and theta1 for geometry
    beta1 = np.atan(((h1/2) + r2)/L01)
    theta1 = np.deg2rad(shockRelation(M=M0sup, B=np.rad2deg(beta1)))
    
    # solve through the first shock
    M1n = M0sup*np.sin(beta1)
    M2n, Pty1_Ptx1, _, _, _ = normshock(M1n, y)
    M2 = M2n/np.sin(beta1 - theta1)
    
    # fsolve gross stuff to find L015, theta2, ad beta2
    def tb2(theta2):
        beta2 = np.deg2rad(shockRelation(theta=theta2, M=M2))
        theta2 = np.deg2rad(theta2)
        
        # geometric restriction
        term1 = (hr - L01*np.tan(theta1 + theta2))/(np.tan(theta1) - np.tan(theta1 + theta2))
        
        # shock angle restriction
        term2 = (h1 + hr - L01*np.tan(beta2 + theta1))/(np.tan(theta1) - np.tan(beta2 + theta1))
        
        return term1 - term2
    
    # solve for theta 2 with function and fsolve
    theta2 = np.deg2rad(fsolve(tb2, 0.1)[0])
    beta2 = np.deg2rad(shockRelation(theta=np.rad2deg(theta2), M=M2))
    L0151 = (hr - L01*np.tan(theta1 + theta2))/(np.tan(theta1) - np.tan(theta1 + theta2))
    L0152 = (h1 + hr - L01*np.tan(beta2 + theta1))/(np.tan(theta1) - np.tan(beta2 + theta1))
    
    # pressure recovery and norm shock
    M3n = M2*np.sin(beta2)
    M4n, Pty2_Ptx2, _, _, _ = normshock(M3n, y)
    M3 = M4n/np.sin(beta2 - theta2)
    M4, Pty3_Ptx3, _, _, _ = normshock(M3, y)
    Prec = Pty1_Ptx1*Pty2_Ptx2*Pty3_Ptx3
    
    # prints
    print("")
    print(f"r2={r2:.3f} in; A1_A2={A1_A2:.2f}; A1={A1:.3f} in^2; h1={h1:.2f} in; hr={hr:.2f} in")
    print(f"L={L:.2f} in; L12={L12:.2f} in; L01={L01:.2f} in")
    print("")
    print(f"beta1={np.rad2deg(beta1):.3f} deg; theta1={np.rad2deg(theta1):.3f} deg")
    print(f"M2={M2:.3f}")
    print(f"L015(1)={L0151:.4f} in; L015(2)={L0152:.4f} in")
    print(f"beta2={np.rad2deg(beta2):.3f} deg; theta2={np.rad2deg(theta2):.3f} deg")
    print("")
    
    # plotting
    if flag == True:
        plt.figure()
        plt.plot([0, L01], [0, 0], "k-")
        plt.plot([L01, L01 + L12], [0, 0], "k-")
        plt.plot([L01, L01], [hr, hr + h1], "r--", label="A1")
        plt.plot([L01, L], [hr, 0], "k-")
        plt.plot([L01, L], [hr + h1, 2*r2], "k-")
        plt.plot([L, L], [0, 2*r2], "r-", label="A2")
        plt.plot([0, L01], [0, L01*np.tan(beta1)])
        plt.plot([0, L0152], [0, L0152*np.tan(theta1)], "k-")
        plt.plot([L0152, L01], [L0152*np.tan(theta1), L0152*np.tan(theta1) + (L01 - L0152)*np.tan(theta2 + theta1)], "k-")
        plt.plot([L0152, L01], [L0152*np.tan(theta1), L0152*np.tan(theta1) + (L01 - L0152)*np.tan(theta1 + beta2)])
        plt.xlim(-4, 180); plt.ylim(0, 50)
        plt.xlabel("Length (in.)"); plt.ylabel("Height (in.)")
        plt.gca().set_aspect("equal")
        plt.legend(loc="upper left")
        plt.grid(); plt.tight_layout()
    
    return np.abs(L0151 - L0152), Prec, L0151, L0152, np.rad2deg(theta2), L01
    

if __name__ == "__main__":
    # engine fan face area calculation
    MFP = compflow(0.9, y)[4]/1.28758
    mdot = 200.
    Tt0 = uu.degF2r(std.T(40000.))/compflow(0.9, y)[0]
    Pt0 = std.P(40000.)/compflow(0.9, y)[1]*uu.psf2psi
    A = mdot*np.sqrt(Tt0)/MFP/Pt0
    A2 = A
    print(f"A2: {A2:.3f}")
    
    # design point
    M0sub = 0.8
    M0sup = 2.5
    M2f = 0.5
    AR = 0.237
    LD = 3.65
    tol, Prec, L0151, L0152, theta2, L01 = inlet(A2, M0sub, M0sup, M2f, AR, LD, True)
    print(tol, Prec)
    plt.show()
    
    flag = False
    # flag = True
    if flag == True:
        # range of vals to study
        LDs = np.linspace(1., 5., 50)
        ARs = np.linspace(0.2, 2., 51)
        
        # loop through and determine prec
        precs = np.empty((len(LDs), len(ARs)))
        for i, LD in enumerate(LDs):
            for k, AR in enumerate(ARs):
                # solve inlet
                tol, Prec, L0151, L0152, theta2, L01 = inlet(A2, M0sub, M0sup, M2f, AR, LD, False)
                
                # check that fsolve isn't a lying turd
                if tol > 0.000001:
                    Prec = 0.0
                elif L0151 < 0.0:
                    Prec = 0.0
                elif L0152 < 0.0:
                    Prec = 0.0
                elif L0151 > L01:
                    Prec = 0.0
                elif L0152 > L01:
                    Prec = 0.0
                elif theta2 > 35.0:
                    Prec = 0.0
                elif theta2 < 0.:
                    Prec = 0.0

                # store vals
                precs[i, k] = Prec
                
        # contour
        X, Y = np.meshgrid(LDs, ARs)
        plt.figure()
        plt.contourf(X, Y, precs.T, 100, cmap="magma")
        cbar = plt.colorbar()
        cbar.set_label(r"$P_{t1}/P_{t0}$")
        plt.xlabel(r"$L/D_2$"); plt.ylabel("AR")
        plt.title("Inlet Design Limits; $M_{0_{sub}}=0.8$, $M_{0_{sup}}=2.5$")
        plt.xlim(1, 5); plt.ylim(0.2, 2)
        plt.grid(); plt.tight_layout()
        plt.show()