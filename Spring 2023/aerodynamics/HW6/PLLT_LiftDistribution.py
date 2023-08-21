# aero PLLT code
# adapted from Brice <3

import sys
import numpy as np
import scipy.integrate as spint 
import matplotlib.pyplot as plt

sys.path.append("C:\\Users\\spbro\\SchoolStuff\\aeroComputing\\lib\\")
import pyvot as pyv


def TAFT(NACA, c=1):
    # airfoil info from name
    m = int(NACA[0])*0.01
    p = int(NACA[1])*0.1
    h = int(NACA[2:])*0.01

    # Define Piece wise functions
    def dz_dx_Cl(theta):
        if theta <= np.arccos(1-2*p/c):
            return ((2*m/p)-(m*c/p**2)*(1-np.cos(theta)))*(np.cos(theta)-1)
        else:
            return ((2*p*m - m*c +m*c*np.cos(theta))/(1-p)**2)*(np.cos(theta)-1)

    def dz_dx_A1(theta):
        if theta <= np.arccos(1-2*p/c):
            return ((2*m/p)-(m*c/p**2)*(1-np.cos(theta)))*np.cos(theta)
        else:
            return ((2*p*m - m*c +m*c*np.cos(theta))/(1-p)**2)*np.cos(theta)

    def dz_dx_A2(theta):
        if theta <= np.arccos(1-2*p/c):
            return ((2*m/p)-(m*c/p**2)*(1-np.cos(theta)))*np.cos(2*theta)
        else:
            return ((2*p*m - m*c +m*c*np.cos(theta))/(1-p)**2)*np.cos(2*theta)

    # Compute integrals/coeffs
    a_L0 = - spint.quad(dz_dx_Cl,0,np.pi)[0]/np.pi
    A_1 = (2/np.pi)*spint.quad(dz_dx_A1,0,np.pi)[0]
    A_2 = (2/np.pi)*spint.quad(dz_dx_A2,0,np.pi)[0]

    # solve TAFT and make fxns
    def CL(alpha):
        return 2*np.pi*(alpha - a_L0)
    def CM_LE(alpha): 
        return - (CL(alpha)/4 + np.pi*(A_1-A_2)/4)
    
    return a_L0, CL, CM_LE


def LiftDistribution(a, a0, b, c, Vinf, S=None, N=50, N_p=100): 
    # make chord a function
    chord_length = c
    def def_c(theta):
        return chord_length
    c = def_c
    
    # Convert for calcs
    a_rad = np.radians(a)
    a0_rad = np.radians(a0)
    
    # Set theta
    theta = np.linspace(0,np.pi,N)
    theta_p = np.linspace(0,np.pi,N_p)

    # set up span stuff
    y = np.array([-b/2*np.cos(t) for t in theta])
    bs = np.linspace(-b/2, b/2, N_p)
    
    # make vars
    Gamma = np.zeros(N_p)
    A = np.zeros((N,N),dtype=float)
    sol = np.empty(N)
    sol = np.full(N, a_rad - a0_rad)
    
    # get coeffs matrix
    for i, th in enumerate(theta):
        for n in range(1,N+1):
            A[i, n-1] = 2*b /(np.pi*c(y[i])) * np.sin(n*th)
            if np.sin(th) != 0:
                A[i, n-1] += n*np.sin(n*th)/np.sin(th)
            else: 
                A[i, n-1] += n*n*np.cos(n*th)/np.cos(th)

    # PLLT coeffs vector
    An = pyv.gaussPivot(A.copy(), sol.copy())
    
    # get gamma
    for i, th in enumerate(theta_p):
        sig=0
        for n in range(1, N+1):
            sig += An[n-1]*np.sin(n*th)
        Gamma[i] = 2*b*Vinf*sig

    # get CL
    Cl = An[0]*np.pi*(b**2)/S
    Cl2 = 2*np.trapz(Gamma, x=bs)/(Vinf*S)
    Cl = (Cl2+Cl)/2

    # get CDi
    sumAs = 0
    for k in range(1, N):
        sumAs += k*(An[k]/An[0])**2
    Cdi = (np.pi*(b**2)/S)*(An[0]**2)*(1 + sumAs)

    return [theta_p, Gamma, Cl, bs, Cdi]


if __name__ == '__main__':
    
    Airfoil = '4412'
    a = 5   # test aoa
    a0 = np.degrees(TAFT(Airfoil)[0])  # 0 lift aoa from TAFT
    b = 5
    c = 1
    Vinf = 10
    S = b*c
    
    t, g, Cl, bs, Cdi = LiftDistribution(a, a0, b, c, Vinf, S)

    # Plotting
    plt.figure()
    plt.plot(bs, g)
    plt.title('NACA{0}\n b = {1:.1f}, c = {2:.1f}, $C_L$ = {3:.4f}, $CD_i$ = {4:.4f}'.format(Airfoil, b, c, Cl, Cdi))
    plt.xlabel('y')
    plt.ylabel('Gamma')
    plt.grid()
    plt.show()