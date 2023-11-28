# AEEM4063 Project 2 Turbine Design Code
# Slade Brooks
# brooksl@mail.uc.edu
# this thing doin numbers (i hope)

import numpy as np
from numpy import sin, cos, tan, arcsin as asin, arccos as acos, arctan as atan


# given vals
cp = 1148                               # cp of hot air, J/kgK
y = 1.333                               # specific heat ratio of hot air
T01 = 2275.                             # turbine inlet temp, K
P01 = 0.                                # turbine inlet pressure, bar
etaInfT = 0.92                          # turbine polytropic eff
etaM = 0.99                             # mechanical efficiency
wComp = 0.                              # compressor work, J/kg
mdotA = 280.                            # air mass flow, kg/s
mdotF = 0.                              # fuel mass flow, kg/s
pTurb = wComp*mdotA*(1 - etaM)          # power required, W
wTurb = pTurb/(mdotA + mdotF)           # turbine specific work, J/kg
n = 0.                                  # rotational speed, rev/s
omega = 2*np.pi*n                       # rotational speed, rad/s
Ca1 = 0.                                # axial speed into turbine, m/s
alpha1 = np.deg2rad(0.)                 # inlet rel gas angle, rad

# design params
lambdam = 0.                            # degree of reaction @ rm
phi = 0.78                              # flow coeff @ rm
psi = 3.3                               # blade loading coeff @ rm
nstages = 2.                            # num stages
ws = wTurb/nstages                      # specific work per stage


# -------------------- Stage 1 Mean --------------------
# solve for inlet speeds
s1C1 = Ca1/cos(alpha1)
s1Cw1 = Ca1/sin(alpha1)
s1Ca1 = Ca1

# solve for gas angles from design params
s1b2 = atan((1/(2*phi))*(psi/2 - 2*lambdam))
s1a2 = atan(tan(s1b2) + 1/phi)
s1b3 = atan((1/(2*phi))*(psi/2 + 2*lambdam))
s1a3 = atan(tan(s1b3) - 1/phi)

# solve for delta T013
s1dT013 = ws/cp

# solve rest of velocity triangles
s1Ca3 = s1Ca1
s1C3 = s1Ca3/cos(s1a3)
s1C2 = np.sqrt(2*cp*(((s1C3)**2)/(2*cp) + s1dT013))
s1Cw2 = s1C2*sin(s1a2)
s1Ca2 = s1C2*cos(s1a2)
s1Cw3 = s1Ca3*tan(s1a3)

# solve for T03
s1T01 = T01
s1T03 = s1T01 - s1dT013

# solve for stage pressure ratio, then stag pres at exit
m1m = etaInfT*(y - 1)/y
s1T0ratio = s1T03/s1T01
s1P0ratio = s1T0ratio**(1/m1m)
s1P01 = P01
s1P03 = s1P01*s1P0ratio

# calculate mean radius from U
s1U = np.sqrt((2*cp*s1dT013)/psi)
s1rm = s1U/omega