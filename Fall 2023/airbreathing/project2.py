# AEEM4063 Project 2 Turbine Design Code
# Slade Brooks
# brooksl@mail.uc.edu
# this thing doin numbers (i hope)

import numpy as np
from numpy import sin, cos, tan, arcsin as asin, arccos as acos, arctan as atan


# given vals
cp = 1148                               # cp of hot air, J/kgK
y = 1.333                               # specific heat ratio of hot air
iT01 = 2275.                            # turbine inlet temp, K
iP01 = 0.                               # turbine inlet pressure, bar
etaInfT = 0.92                          # turbine polytropic eff
etaM = 0.99                             # mechanical efficiency
wComp = 0.                              # compressor work, J/kg
mdotA = 280.                            # air mass flow, kg/s
mdotF = 0.                              # fuel mass flow, kg/s
pTurb = wComp*mdotA*(1 - etaM)          # power required, W
wTurb = pTurb/(mdotA + mdotF)           # turbine specific work, J/kg
n = 0.                                  # rotational speed, rev/s
omega = 2*np.pi*n                       # rotational speed, rad/s
Ca1 = 275.                              # axial speed into turbine, m/s
alpha1 = np.deg2rad(0.)                 # inlet rel gas angle, rad


# calc gas angles from design params
def gasAngles(phi, psi, lam):
    b2 = atan(1/(2*phi)*(psi/2 - 2*lam))
    a2 = atan(tan(b2) + 1/phi)
    b3 = atan(1/(2*phi)*(psi/2 + 2*lam))
    a3 = atan(tan(b3) - 1/phi)
    return b2, a2, b3, a3


# mean-line calcs for any stage
def meanCalcs(alpha1, Ca1, phi, psi, lam, T01, P01, omega):
    # solve for inlet speeds
    C1 = Ca1/cos(alpha1)
    if alpha1 != 0:
        Cw1 = Ca1/sin(alpha1)
    else:
        Cw1 = 0.

    # solve for U
    U = Ca1/phi

    # solve for delta T
    dT013 = (U**2)*psi/(2*cp)

    # solve stage specific work
    w = cp*dT013

    # get gas angles
    b2, a2, b3, a3 = gasAngles(phi, psi, lam)
    
    # solve rest of velocity trinangle
    Ca2 = Ca1
    Cw2 = Ca2*tan(a2)
    C2 = np.sqrt(Ca2**2 + Cw2**2)
    Ca3 = Ca1
    Cw3 = Ca3*tan(a3)
    C3 = np.sqrt(Ca3**2 + Cw3**2)
    V2 = Ca2/cos(b2)
    V3 = Ca3/cos(b3)
    
    # use polytropic eff to solve m1m, then for Trat, then for Prat
    m1m = etaInfT*(y - 1)/y
    T03 = T01 - dT013
    T0rat = T03/T01
    P0rat = T0rat**(1/m1m)
    P03 = P0rat*P01
    
    # calc rm from U
    rm = U/omega
    
    return C1, Ca1, Cw1, C2, Ca2, Cw2, C3, Ca3, Cw3, V2, V3, b2, a2, b3, a3, T03, P03, w, rm


# mean-line calcs for stage with axial outflow
def lastStageCalc(alpha1, Ca1, w, T01, P01, omega, r):
    # solve for inlet speeds
    C1 = Ca1/cos(alpha1)
    if alpha1 != 0:
        Cw1 = Ca1/sin(alpha1)
    else:
        Cw1 = 0.

    # solve for U
    U = omega*r
    
    # solve for deltaT to get required word
    dT013 = w/cp
    
    # solve gas angles and velocity triangles while enforcing no whirl at exit
    a3 = 0.
    Cw3 = 0.
    Ca3 = Ca1
    C3 = Ca3
    Ca2 = Ca1
    a2 = atan(w/(U*Ca2))
    Cw2 = Ca2*tan(a2)
    C2 = np.sqrt(Ca2**2 + Cw2**2)
    b2 = atan((Cw2 - U)/Ca2)
    V2 = np.sqrt(Ca2**2 + (Cw2 - U)**2)
    b3 = atan(U/Ca3)
    V3 = np.sqrt(U**2 + Ca3**2)
    
    # use polytropic eff to solve m1m, then for Trat, then for Prat
    m1m = etaInfT*(y - 1)/y
    T03 = T01 - dT013
    T0rat = T03/T01
    P0rat = T0rat**(1/m1m)
    P03 = P0rat*P01
    
    # check phi, psi, and degree of reaction
    phi = Ca1/U
    psi = (2*cp*dT013)/(U**2)
    lam = Ca1/(2*U)*(tan(b3) - tan(b2))
    
    return C1, Ca1, Cw1, C2, Ca2, Cw2, C3, Ca3, Cw3, V2, V3, b2, a2, b3, a3, T03, P03, phi, psi, lam


# -------------------- Stage 1 Mean --------------------
# design params
lambdam = 0.5                           # degree of reaction @ rm
phi = 0.78                              # flow coeff @ rm
psi = 3.3                               # blade loading coeff @ rm

# run meanline calcs
s1C1, s1Ca1, s1Cw1, s1C2, s1Ca2, s1Cw2, s1C3, s1Ca3, s1Cw3, s1V2, s1V3, s1b2, s1a2, s1b3, s1a3, s1T03, s1P03, s1w, s1rm\
    = meanCalcs(alpha1, Ca1, phi, psi, lambdam, iT01, iP01, omega)


# -------------------- Stage 2 Mean --------------------
# design params
lambdam = 0.5                           # degree of reaction @ rm
phi = 0.78                              # flow coeff @ rm
psi = 3.3                               # blade loading coeff @ rm

# run meanline calcs
s2C1, s2Ca1, s2Cw1, s2C2, s2Ca2, s2Cw2, s2C3, s2Ca3, s2Cw3, s2V2, s2V3, s2b2, s2a2, s2b3, s2a3, s2T03, s2P03, s2w, s2rm\
    = meanCalcs(s1a3, s1Ca3, phi, psi, lambdam, s1T03, s1P03, omega)


# -------------------- Stage 3 Mean --------------------
wleft = wTurb - (s1w + s2w)
rinc = s2rm/s1rm
s3rm = rinc*s2rm
s3C1, s3Ca1, s3Cw1, s3C2, s3Ca2, s3Cw2, s3C3, s3Ca3, s3Cw3, s3V2, s3V3, s3b2, s3a2, s3b3, s3a3, s3T03, s3P03, s3phi, \
    s3psi, s3lam = lastStageCalc(s2a3, s2Ca3, wleft, s2T03, s2P03, omega, s3rm)