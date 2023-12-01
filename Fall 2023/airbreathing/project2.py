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
iP01 = 21.278                           # turbine inlet pressure, bar
etaInfT = 0.92                          # turbine polytropic eff
etaM = 0.99                             # mechanical efficiency
wComp = 437774.0                        # compressor work, J/kg
wFan = 18200.6                          # fan work, J/kg
mdotA = 280.                            # fan mass flow, kg/s
mdotH = 40.                             # core mass flow, kg/s
mdotF = 1.862                           # fuel mass flow, kg/s
mdotT = mdotF + mdotH                   # total mass flow, kg/s
n = 172.32                              # rotational speed, rev/s
omega = 2*np.pi*n                       # rotational speed, rad/s
Ca1 = 275.                              # axial speed into turbine, m/s
alpha1 = np.deg2rad(0.)                 # inlet rel gas angle, rad
R = 287.                                # r, J/kgK

# calc turb specific work required, J/kg
wTurb = (mdotA*wFan + mdotH*wComp)/(etaM*mdotT)


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
    Cw1 = Ca1*tan(alpha1)

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
    Cw1 = Ca1*tan(alpha1)

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
    lam = (Ca1/(2*U))*(tan(b3) - tan(b2))
    
    return C1, Ca1, Cw1, C2, Ca2, Cw2, C3, Ca3, Cw3, V2, V3, b2, a2, b3, a3, T03, P03, phi, psi, lam


# find root and tip radii from annulus area and height
def rootTipRadii(mdot, C, P0, T0, rm):
    # get static T
    Ts = T0 - (C**2)/(2*cp)
    
    # get static P
    Ps = P0/((T0/Ts)**(y/(y - 1)))
    
    # find area then height then radii
    area = (mdot*R*Ts)/((Ps*100000)*C) # m^2
    h = area/(2*np.pi*rm) # m
    rr = rm - h/2
    rt = rm + h/2
    
    return rr, rt


# scale stuff for root or tip
def rootTipCalc(ri, rm, omega, Cwm, Ca, b3):
    # find U
    U = omega*ri
    
    # mean/i ratio
    rrat = rm/ri
    
    # find whirl from ratio
    Cw2 = Cwm*rrat
    
    # get a2, b2, c2, v2
    Ca2 = Ca
    a2 = atan(Cw2/Ca2)
    b2 = atan((Cw2 - U)/Ca2)
    C2 = np.sqrt(Ca2**2 + Cw2**2)
    V2 = np.sqrt((Cw2 - U)**2 + Ca2**2)
    
    # find phi, psi, lambda  
    lam = (Ca/(2*U))*(tan(b3) - tan(b2))
    phi = Ca/U
    psi = ((2*Ca)/U)*(tan(b2) + tan(b3))
    
    return C2, Ca2, Cw2, V2, b2, a2, phi, psi, lam


# -------------------- Stage 1 Mean --------------------
# design params
lambdam = 0.5                           # degree of reaction @ rm
phi = 0.78                              # flow coeff @ rm
psi = 3.3                               # blade loading coeff @ rm

# run meanline calcs
s1C1, s1Ca1, s1Cw1, s1C2, s1Ca2, s1Cw2, s1C3, s1Ca3, s1Cw3, s1V2, s1V3, s1b2, s1a2, s1b3, s1a3, s1T03, s1P03, s1w, s1rm \
    = meanCalcs(alpha1, Ca1, phi, psi, lambdam, iT01, iP01, omega)


# -------------------- Stage 2 Mean --------------------
# design params
lambdam = 0.5                           # degree of reaction @ rm
phi = 0.78                              # flow coeff @ rm
psi = 3.3                               # blade loading coeff @ rm

# run meanline calcs
s2C1, s2Ca1, s2Cw1, s2C2, s2Ca2, s2Cw2, s2C3, s2Ca3, s2Cw3, s2V2, s2V3, s2b2, s2a2, s2b3, s2a3, s2T03, s2P03, s2w, s2rm \
    = meanCalcs(s1a3, s1Ca3, phi, psi, lambdam, s1T03, s1P03, omega)


# -------------------- Stage 3 Mean --------------------
wleft = wTurb - (s1w + s2w)
print(wTurb)
print(s1w, s2w)
print(wleft)
rinc = s2rm/s1rm
s3rm = rinc*s2rm
s3C1, s3Ca1, s3Cw1, s3C2, s3Ca2, s3Cw2, s3C3, s3Ca3, s3Cw3, s3V2, s3V3, s3b2, s3a2, s3b3, s3a3, s3T03, s3P03, s3phi, \
    s3psi, s3lam = lastStageCalc(s2a3, s2Ca3, wleft, s2T03, s2P03, omega, s3rm)
    

# -------------------- Root and Tips --------------------
# get root and tip radii for each
s1rr, s1rt = rootTipRadii(mdotT, s1C1, iP01, iT01, s1rm)
s2rr, s2rt = rootTipRadii(mdotT, s2C1, s1P03, s1T03, s2rm)
s3rr, s3rt = rootTipRadii(mdotT, s3C1, s2P03, s2T03, s3rm)

# stage 1 root and tip
s1C2r, s1Ca2r, s1Cw2r, s1V2r, s1b2r, s1a2r, s1phir, s1psir, s1lamr = rootTipCalc(s1rr, s1rm, omega, s1Cw2, s1Ca1, s1b3)
s1C2t, s1Ca2t, s1Cw2t, s1V2t, s1b2t, s1a2t, s1phit, s1psit, s1lamt = rootTipCalc(s1rt, s1rm, omega, s1Cw2, s1Ca1, s1b3)

# stage 2 root and tip
s2C2r, s2Ca2r, s2Cw2r, s2V2r, s2b2r, s2a2r, s2phir, s2psir, s2lamr = rootTipCalc(s2rr, s2rm, omega, s2Cw2, s2Ca1, s2b3)
s2C2t, s2Ca2t, s2Cw2t, s2V2t, s2b2t, s2a2t, s2phit, s2psit, s2lamt = rootTipCalc(s2rt, s2rm, omega, s2Cw2, s2Ca1, s2b3)

# stage 3 root and tip
s3C2r, s3Ca2r, s3Cw2r, s3V2r, s3b2r, s3a2r, s3phir, s3psir, s3lamr = rootTipCalc(s3rr, s3rm, omega, s3Cw2, s3Ca1, s3b3)
s3C2t, s3Ca2t, s3Cw2t, s3V2t, s3b2t, s3a2t, s3phit, s3psit, s3lamt = rootTipCalc(s3rt, s3rm, omega, s3Cw2, s3Ca1, s3b3)

# design lims check
phis = [s1phir, s1phit, s2phir, s2phit, s3phir, s3phit, s3phi]
phiss = ["s1phir", "s1phit", "s2phir", "s2phit", "s3phir", "s3phit", "s3phi"]
psis = [s1psir, s1psit, s2psir, s2psit, s3psir, s3psit, s3psi]
psiss = ["s1psir", "s1psit", "s2psir", "s2psit", "s3psir", "s3psit", "s3psi"]
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Design Params Check")
for k in range(len(phis)):
    if phis[k] < 0.78:
        print(f"{phiss[k]} = {phis[k]:.3f} (<0.78)")
    if psis[k] > 3.3:
        print(f"{psiss[k]} = {psis[k]:.2f} (>3.3)")
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


# -------------------- Print Results --------------------
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('stage 1 mean','\n',\
      '\n-- gas angles --', \
      '\nb2_s1:', np.degrees(s1b2), \
      '\nb3_s1:', np.degrees(s1b3), \
      '\na2_s1:', np.degrees(s1a2), \
      '\na3_s1:', np.degrees(s1a3), \
      '\n-- nozzle inlet --', \
      '\nCa1_s1:', s1Ca1, \
      '\nCw1_s1:', s1Cw1, \
      '\nC1_s1:', s1C1, \
      '\n-- rotor inlet --', \
      '\nCa2_s1:', s1Ca2, \
      '\nCw2_s1:', s1Cw2, 
      '\nC2_s1:', s1C2, \
      '\nV2_s1:', s1V2, \
      '\n-- stage 1 exit --', \
      '\nCa3_s1:', s1Ca3, \
      '\nCw3_s1:', s1Cw3, \
      '\nC3_s1:', s1C3, \
      '\nV3_s1:', s1V3, \
      '\n-- stage 1 exit flow conditions --', \
      '\nT03_s1:', s1T03, \
      '\nP03_s1:', s1P03)
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('--stage 1 root--', \
      '\nCw2r_s1:', s1Cw2r, \
      '\nC2r_s1:', s1C2r, \
      '\nV2r_s1:', s1V2r, \
      '\nphir_s1:', s1phir, \
      '\npsir_s1:', s1psir, \
      '\nlamr_s1:', s1lamr, \
      '\n-- stage 1 tip--', \
      '\nCw2t_s1:', s1Cw2t, \
      '\nC2t_s1:', s1C2t, \
      '\nV2t_s1:', s1V2t, \
      '\nphit_s1:', s1phit, \
      '\npsit_s1:', s1psit, \
      '\nlamt_s1:', s1lamt, \
      '\n-- radii --', \
      '\nr_r_s1:', s1rr, \
      '\nr_t_s1:', s1rt)
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('stage 2 mean','\n',
      '\n-- gas angles --', \
      '\nb2_s2:', np.degrees(s2b2), \
      '\nb3_s2:', np.degrees(s2b3), \
      '\na2_s2:', np.degrees(s2a2), \
      '\na3_s2:', np.degrees(s2a3), \
      '\n-- nozzle inlet --', \
      '\nCa2_s2:', s2Ca1, \
      '\nCw1_s2:', s2Cw1, \
      '\nC1_s2:', s2C1, \
      '\n-- rotor inlet --', \
      '\nCa2_s2:', s2Ca2, \
      '\nCw2_s2:', s2Cw2, \
      '\nC2_s2:', s2C2, \
      '\nV2_s2:', s2V2, \
      '\n-- stage 2 exit --', \
      '\nCa3_s2:', s2Ca3, \
      '\nCw3_s2:', s2Cw3, \
      '\nC3_s2:', s2C3, \
      '\nV3_s2:', s2V3, \
      '\n-- stage 2 exit flow conditions --', \
      '\nT03_s2:', s2T03, \
      '\nP03_s2:', s2P03)
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('-- stage 2 root --',\
      '\nCw2r_s2:', s2Cw2r, \
      '\nC2r_s2:', s2C2r, \
      '\nV2r_s2:', s2V2r, \
      '\nphir_s2:', s2phir, \
      '\npsir_s2:', s2psir, \
      '\nlamr_s2:', s2lamr, \
      '\n-- stage 2 tip --',\
      '\nCw2t_s2:', s2Cw2t, \
      '\nC2t_s2:', s2C2t, \
      '\nV2t_s2:', s2V2t, \
      '\nphit_s2:', s2phit, \
      '\npsit_s2:', s2psit, \
      '\nlamt_s2:', s2lamt, \
      '\n-- radii --',\
      '\nr_r_s2:', s2rr, \
      '\nr_t_s2:', s2rt)
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('final stage mean','\n',
      '\n-- gas angles --', \
      '\nb2_s3:', np.degrees(s3b2), \
      '\nb3_s3:', np.degrees(s3b3), \
      '\na2_s3:', np.degrees(s3a2), \
      '\na3_s3:', np.degrees(s3a3), \
      '\n-- nozzle inlet --', \
      '\nCw1_s3:', s3Cw1, \
      '\nC1_s3:', s3C1, \
      '\n-- rotor inlet --', \
      '\nCa2_s3:', s3Ca2, \
      '\nCw2_s3:', s3Cw2, \
      '\nC2_s3:', s3C2, \
      '\nV2_s3:', s3V2, \
      '\n-- stage 3 exit --', \
      '\nCa3_s3:', s3Ca3, \
      '\nCw3_s3:', s3Cw3, \
      '\nC3_s3:', s3C3, \
      '\nV3_s3:', s3V3, \
      '\n-- stage 3 exit flow conditions --', \
      '\nT03_s3:', s3T03, \
      '\nP03_s3:', s3P03, \
      '\n-- parameters --', \
      '\nphi_s3:', s3phi, \
      '\npsi_s3:', s3psi, \
      '\nlam_s3:', s3lam)
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('-- stage 3 root --',\
      '\nCw2r_s3:', s3Cw2r, \
      '\nC2r_s3:', s3C2r, \
      '\nV2r_s3:', s3V2r, \
      '\nphir_s3:', s3phir, \
      '\npsir_s3:', s3psir, \
      '\nlamr_s3:', s3lamr, \
      '\n-- stage 3 tip --',\
      '\nCw2t_s3:', s3Cw2t, \
      '\nC2t_s3:', s3C2t, \
      '\nV2t_s3:', s3V2t, \
      '\nphit_s3:', s3phit, \
      '\npsit_s3:', s3psit, \
      '\nlamt_s3:', s3lamt, \
      '\n-- radii --',\
      '\nr_r_s3:', s3rr, \
      '\nr_t_s3:', s3rt)
print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~')