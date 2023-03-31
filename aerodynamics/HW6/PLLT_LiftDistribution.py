# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 01:23:34 2023

@author: cycon
stolen by slade :)
"""
import sys
import numpy as np
import scipy.integrate as spint 
import matplotlib.pyplot as plt

sys.path.append("C:\\Users\\spbro\\SchoolStuff\\aeroComputing\\lib\\")
import pyvot as MM


def TAFT(NACA, c=1):
    # Cant check NACA for correcft input yet
    # Set up parameters 
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

    # Compute Integrals
    a_L0 = - spint.quad(dz_dx_Cl,0,np.pi)[0]/np.pi
    A_1 = (2/np.pi)*spint.quad(dz_dx_A1,0,np.pi)[0]
    A_2 = (2/np.pi)*spint.quad(dz_dx_A2,0,np.pi)[0]

    # Compute Coefficients
    def CL(alpha):
        return 2*np.pi*(alpha - a_L0)

    def CM_LE(alpha): 
        return - (CL(alpha)/4 + np.pi*(A_1-A_2)/4)
    
    return a_L0, CL, CM_LE


# Function assuming a constant airfoil acros cross section, sym about midpoint
def LiftDistribution(a, a0, b, c, Vinf, S=None, N=50, N_p=100): 
    # If c is passed as a scalar, covnert to a function
    if not callable(c):
        chord_length = c;
        def def_c(theta):
            return chord_length
        c = def_c
    
    # Convert to radians for calculations
    a_rad = np.radians(a)
    a0_rad = np.radians(a0)
    
    # Set up our thetas
    theta = np.linspace(0,np.pi,N)
    theta_p = np.linspace(0,np.pi,N_p)
    
    y = np.linspace(-b/2, b/2, N)
    
    # Initialize arrays/Matrices
    Gamma = np.zeros(N_p)
    A = np.zeros((N,N),dtype=float)
    sol = np.empty(N)
    sol = np.full(N, a_rad - a0_rad)
    
    # Calculate the matrix to find the coefficients An
    for i, th in enumerate(theta):
        for n in range(1,N+1):
            A[i, n-1] = 2*b /(np.pi*c(y[i])) * np.sin(n*th)
            if np.sin(th) != 0:
                A[i, n-1] += n*np.sin(n*th)/np.sin(th)
            else: 
                A[i, n-1] += n*n*np.cos(n*th)/np.cos(th)

    # Get the coefficient vector
    An = MM.gaussPivot(A.copy(), sol.copy())
   
    # # Get error in case it is needed
    # res = MM.accuracy(A, An, sol)
    
    # Calculate Gamma
    for i, th in enumerate(theta_p):
        sig=0
        for n in range(1,N+1):
            sig += An[n-1]*np.sin(n*th)
        Gamma[i] = 2*b*Vinf*sig
    
    if S != None:
        # Calculate Lift Coefficient
        CL = 2*np.trapz(Gamma, dx=np.pi/N_p)/(Vinf*S)
        
    return [theta_p, Gamma] if S == None else [theta_p, Gamma, CL]

def chord_Lin_func_maker(max_c, min_c, b):
    # will work according to 
    
    def chord_funct_y(y):
        c = abs(y)*(-min_c/(2*b)) + max_c
        return c
    
    def chord_funct_theta(th):
        c = (max_c-min_c)*np.sin(th) + min_c
        return c
    
    return chord_funct_theta, chord_funct_y


def plotPlaniform(b,c_func):
    # Assuming c is callable for now
    # Works on theta:
    th = np.linspace(0,np.pi,100,endpoint=True)
    # y = -b*np.cos(th)/2
    y = np.linspace(-b/2, b/2, 100, endpoint=True)
    c = c_func(y) # Chord lengths
    
    upperC = c/2
    lowerC = -c/2
    leftTip = np.array([[y[0], upperC[0]],[y[0], lowerC[0]]])
    rightTip = np.array([[y[-1], upperC[-1]],[y[-1], lowerC[-1]]])
    
    plt.figure()
    plt.plot(y,upperC, y, lowerC, leftTip[:,0], leftTip[:,1], rightTip[:,0], rightTip[:,1])
    c_max = np.max(c)
    plt.ylim([-c_max,c_max])
    # plt.xlim([-3*b/4,3*b])
    plt.grid()


if __name__ == '__main__':
    
    Airfoil = '8412'
    # Set up some standard constants
    a = 2   # deg - angle of attack of interest
    a0 = TAFT(Airfoil)[0]  # deg - zero lift angle of attack
    b = 8  # m - span length
    c = 1   # m - chord length 
    Vinf = 1.56 # m/s - freestream velocity
    S = b*c # m^2 - Ref Area (different if chord varies)
    
    t, g, Cl = LiftDistribution(a, a0, b, c, Vinf,S)

    # ---- Plotting Gamma -----
    if True:
        # Plotting
        plt.figure()
        plt.plot(t, g)
        # plt.yticks(np.arange(0,15,2.5))
        plt.title('NACA{}\n b = {}, c = {}, $C_L$ = {}'.format(Airfoil, b, c, Cl))
        plt.xlabel('Theta (rad)')
        plt.ylabel('Gamma')
        plt.grid()
        plt.show()
    
    if False:
        for n in range(5,50):
            t, g, Cl = LiftDistribution(a, a0, b, 1, Vinf,S,n, 100)
        
            # ---- Plotting Gamma -----
            if True:
                # Plotting
                plt.plot(t, g)
                # plt.yticks(np.arange(0,15,2.5))
                plt.xlabel('Theta (rad)')
                plt.ylabel('Gamma')
                plt.title('$C_L$ = {}\nN={}'.format(Cl, n))
                plt.grid()
                plt.pause(0.5)
                plt.clf()
        
        
