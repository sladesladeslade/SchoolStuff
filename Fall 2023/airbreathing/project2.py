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
iP01 = 31.045                           # turbine inlet pressure, bar
etaInfT = 0.92                          # turbine polytropic eff
etaM = 0.99                             # mechanical efficiency
wComp = 419350.5                        # compressor work, J/kg
wFan = 18200.6                          # fan work, J/kg
mdotA = 280.                            # fan mass flow, kg/s
mdotH = 40.                             # core mass flow, kg/s
mdotF = 1.862                           # fuel mass flow, kg/s
mdotT = mdotF + mdotH                   # total mass flow, kg/s
n = 172.32                              # rotational speed, rev/s
omega = 2*np.pi*n                       # rotational speed, rad/s
Ca1 = 300.                              # axial speed into turbine, m/s
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
phi = 0.804                             # flow coeff @ rm
psi = 3.0609                            # blade loading coeff @ rm

# run meanline calcs
s1C1, s1Ca1, s1Cw1, s1C2, s1Ca2, s1Cw2, s1C3, s1Ca3, s1Cw3, s1V2, s1V3, s1b2, s1a2, s1b3, s1a3, s1T03, s1P03, s1w, s1rm \
    = meanCalcs(alpha1, Ca1, phi, psi, lambdam, iT01, iP01, omega)


# -------------------- Stage 2 Mean --------------------
# design params
lambdam = 0.5                           # degree of reaction @ rm
phi = 0.8107                            # flow coeff @ rm
psi = 2.99                              # blade loading coeff @ rm

# run meanline calcs
s2C1, s2Ca1, s2Cw1, s2C2, s2Ca2, s2Cw2, s2C3, s2Ca3, s2Cw3, s2V2, s2V3, s2b2, s2a2, s2b3, s2a3, s2T03, s2P03, s2w, s2rm \
    = meanCalcs(s1a3, s1Ca3, phi, psi, lambdam, s1T03, s1P03, omega)


# -------------------- Stage 3 Mean --------------------
wleft = wTurb - (s1w + s2w)
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Stage Work")
print(f"total turb w: {wTurb/1000:.2f} kJ/kg")
print(f"s1w: {s1w/1000:.2f} kJ/kg, s2w: {s2w/1000:.2f} kJ/kg")
print(f"s3w: {wleft/1000:.2f} kJ/kg")

# get root and tip radii for each
s1rr, s1rt = rootTipRadii(mdotT, s1C1, iP01, iT01, s1rm)
s2rr, s2rt = rootTipRadii(mdotT, s2C1, s1P03, s1T03, s2rm)
s1h = s1rt - s1rr
s2h = s2rt - s2rr
s3rm = 0.98*s2rm
s3C1, s3Ca1, s3Cw1, s3C2, s3Ca2, s3Cw2, s3C3, s3Ca3, s3Cw3, s3V2, s3V3, s3b2, s3a2, s3b3, s3a3, s3T03, s3P03, s3phi, \
    s3psi, s3lam = lastStageCalc(s2a3, s2Ca3, wleft, s2T03, s2P03, omega, s3rm)


# -------------------- Root and Tips --------------------
# get root and tip radii for each
s1rr, s1rt = rootTipRadii(mdotT, s1C1, iP01, iT01, s1rm)
s2rr, s2rt = rootTipRadii(mdotT, s2C1, s1P03, s1T03, s2rm)
s3rr, s3rt = rootTipRadii(mdotT, s3C1, s2P03, s2T03, s3rm)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Stage Geometries")
print(f"hs: {s1h:.4f} m, {s2h:.4f} m, {(s3rt-s3rr):.4f} m")
print(f"rms: {s1rm:.4f} m, {s2rm:.4f} m, {s3rm:.4f} m")
print(f"rts: {s1rt:.4f} m, {s2rt:.4f} m, {s3rt:.4f} m")
print(f"rrs: {s1rr:.4f} m, {s2rr:.4f} m, {s3rr:.4f} m")

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


# outputs
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Results')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Stage 1")
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Mean:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s1a2):.1f} deg, b2: {np.rad2deg(s1b2):.1f} deg")
print(f"a3: {np.rad2deg(s1a3):.1f} deg, b3: {np.rad2deg(s1b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s1C1:.1f} m/s, Ca1: {s1Ca1:.1f} m/s, Cw1: {s1Cw1:.1f} m/s")
print(f"C2: {s1C2:.1f} m/s, Ca2: {s1Ca2:.1f} m/s, Cw2: {s1Cw2:.1f} m/s")
print(f"C3: {s1C3:.1f} m/s, Ca3: {s1Ca3:.1f} m/s, Cw3: {s1Cw3:.1f} m/s")
print("")
print("Root:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s1a2r):.1f} deg, b2: {np.rad2deg(s1b2r):.1f} deg")
print(f"a3: {np.rad2deg(s1a3):.1f} deg, b3: {np.rad2deg(s1b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s1C1:.1f} m/s, Ca1: {s1Ca1:.1f} m/s, Cw1: {s1Cw1:.1f} m/s")
print(f"C2: {s1C2r:.1f} m/s, Ca2: {s1Ca2r:.1f} m/s, Cw2: {s1Cw2r:.1f} m/s")
print(f"C3: {s1C3:.1f} m/s, Ca3: {s1Ca3:.1f} m/s, Cw3: {s1Cw3:.1f} m/s")
print("")
print("Tip:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s1a2t):.1f} deg, b2: {np.rad2deg(s1b2t):.1f} deg")
print(f"a3: {np.rad2deg(s1a3):.1f} deg, b3: {np.rad2deg(s1b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s1C1:.1f} m/s, Ca1: {s1Ca1:.1f} m/s, Cw1: {s1Cw1:.1f} m/s")
print(f"C2: {s1C2t:.1f} m/s, Ca2: {s1Ca2t:.1f} m/s, Cw2: {s1Cw2t:.1f} m/s")
print(f"C3: {s1C3:.1f} m/s, Ca3: {s1Ca3:.1f} m/s, Cw3: {s1Cw3:.1f} m/s")
print("")
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Stage 2")
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Mean:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s2a2):.1f} deg, b2: {np.rad2deg(s2b2):.1f} deg")
print(f"a3: {np.rad2deg(s2a3):.1f} deg, b3: {np.rad2deg(s2b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s2C1:.1f} m/s, Ca1: {s2Ca1:.1f} m/s, Cw1: {s2Cw1:.1f} m/s")
print(f"C2: {s2C2:.1f} m/s, Ca2: {s2Ca2:.1f} m/s, Cw2: {s2Cw2:.1f} m/s")
print(f"C3: {s2C3:.1f} m/s, Ca3: {s2Ca3:.1f} m/s, Cw3: {s2Cw3:.1f} m/s")
print("")
print("Root:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s2a2r):.1f} deg, b2: {np.rad2deg(s2b2r):.1f} deg")
print(f"a3: {np.rad2deg(s2a3):.1f} deg, b3: {np.rad2deg(s2b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s2C1:.1f} m/s, Ca1: {s2Ca1:.1f} m/s, Cw1: {s2Cw1:.1f} m/s")
print(f"C2: {s2C2r:.1f} m/s, Ca2: {s2Ca2r:.1f} m/s, Cw2: {s2Cw2r:.1f} m/s")
print(f"C3: {s2C3:.1f} m/s, Ca3: {s2Ca3:.1f} m/s, Cw3: {s2Cw3:.1f} m/s")
print("")
print("Tip:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s2a2t):.1f} deg, b2: {np.rad2deg(s2b2t):.1f} deg")
print(f"a3: {np.rad2deg(s2a3):.1f} deg, b3: {np.rad2deg(s2b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s2C1:.1f} m/s, Ca1: {s2Ca1:.1f} m/s, Cw1: {s2Cw1:.1f} m/s")
print(f"C2: {s2C2t:.1f} m/s, Ca2: {s2Ca2t:.1f} m/s, Cw2: {s2Cw2t:.1f} m/s")
print(f"C3: {s2C3:.1f} m/s, Ca3: {s2Ca3:.1f} m/s, Cw3: {s2Cw3:.1f} m/s")
print("")
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Stage 3")
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Mean:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s3a2):.1f} deg, b2: {np.rad2deg(s3b2):.1f} deg")
print(f"a3: {np.rad2deg(s3a3):.1f} deg, b3: {np.rad2deg(s3b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s3C1:.1f} m/s, Ca1: {s3Ca1:.1f} m/s, Cw1: {s3Cw1:.1f} m/s")
print(f"C2: {s3C2:.1f} m/s, Ca2: {s3Ca2:.1f} m/s, Cw2: {s3Cw2:.1f} m/s")
print(f"C3: {s3C3:.1f} m/s, Ca3: {s3Ca3:.1f} m/s, Cw3: {s3Cw3:.1f} m/s")
print("")
print("Root:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s3a2r):.1f} deg, b2: {np.rad2deg(s3b2r):.1f} deg")
print(f"a3: {np.rad2deg(s3a3):.1f} deg, b3: {np.rad2deg(s3b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s3C1:.1f} m/s, Ca1: {s3Ca1:.1f} m/s, Cw1: {s3Cw1:.1f} m/s")
print(f"C2: {s3C2r:.1f} m/s, Ca2: {s3Ca2r:.1f} m/s, Cw2: {s3Cw2r:.1f} m/s")
print(f"C3: {s3C3:.1f} m/s, Ca3: {s3Ca3:.1f} m/s, Cw3: {s3Cw3:.1f} m/s")
print("")
print("Tip:")
print("")
print("Gas Angles:")
print(f"a1: {np.rad2deg(alpha1):.1f} deg")
print(f"a2: {np.rad2deg(s3a2t):.1f} deg, b2: {np.rad2deg(s3b2t):.1f} deg")
print(f"a3: {np.rad2deg(s3a3):.1f} deg, b3: {np.rad2deg(s3b3):.1f} deg")
print("")
print("Flow Velocities:")
print(f"C1: {s3C1:.1f} m/s, Ca1: {s3Ca1:.1f} m/s, Cw1: {s3Cw1:.1f} m/s")
print(f"C2: {s3C2t:.1f} m/s, Ca2: {s3Ca2t:.1f} m/s, Cw2: {s3Cw2t:.1f} m/s")
print(f"C3: {s3C3:.1f} m/s, Ca3: {s3Ca3:.1f} m/s, Cw3: {s3Cw3:.1f} m/s")