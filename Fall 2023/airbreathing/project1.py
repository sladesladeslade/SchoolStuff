# AEEM4063 Project 1 Code
# Slade Brooks
# brooksl@mail.uc.edu
# jet emgime :o

import numpy as np
import matplotlib.pyplot as plt


# ---------- Ambient Conditions ----------
P0 = 0.22697        # bar
T0 = 216.78         # K
M0 = 0.85

# ---------- Given Vals ----------
yc = 1.4            # spec heat ratio cold
yh = 1.333          # spec heat ratio hot
T04 = 1560          # turbine inlet temp, K
hfuel = 43100       # fuel heat thing whatever kJ/kg
BPR = 10            # bypass ratio
pif = 1.5           # fan pres ratio
pic = 36            # compres pres ratio
pib = 0.96          # burner pres ratio
etai = 0.98         # inlet eff
etainff = 0.89      # fan polytropic eff
etainfc = 0.90      # compres polytropic eff
etab = 0.99         # 'bustion eff
etainft = 0.90      # turb polytropic eff
etam = 0.99         # mech eff
etaj = 0.99         # nozzle eff


# ---------- Inlet ----------
def inlet(P0, M0, etai, yc):
    # solve for stag P w/ isen relation
    P02 = P0*(1 + etai*(yc - 1)/2*M0**2)**(yc/(yc - 1))

    # solve for stag T - adiabatic so = T00
    T02 = T0*(1 + (yc - 1)/2*M0**2)
    
    return P02, T02


# ---------- Fan ----------
def fan(pif, etainff, P02, T02, yc):
    # solve for P w/ pressure ratio
    P025 = P02*pif
    
    # solve for n-1/n based on polytropic eff
    n1n = (1/etainff)*(yc - 1)/yc
    
    # solve for T from relation with n-1/n
    T025 = T02*pif**n1n
    
    return P025, T025


# ---------- Compressor ----------
def compressor(pic, etainfc, P025, T025, yc):
    # solve for P w/ pressure ratio
    P03 = P025*pic
    
    # solve for n-1/n based on polytropic eff
    n1n = (1/etainfc)*(yc - 1)/yc
    
    # solve for T from relation with n-1/n
    T03 = T025*pic**n1n
    
    return P03, T03


# ---------- Combustor ----------
def combustor(pib, P03, TIT):
    # solve for P from pressure loss
    P04 = P03*pib
    
    # given turbine inlet temp
    T04 = TIT
    
    return P04, T04


# ---------- Turbine ----------
def turbine(etainft, P04, T04, yh):
    # solve for m-1/m based on polytropic eff
    m1m = etainft*(yh - 1)/yh


# ---------- Testing ----------
P02, T02 = inlet(P0, M0, etai, yc)
P025, T025 = fan(pif, etainff, P02, T02, yc)
P03, T03 = compressor(pic, etainfc, P025, T025, yc)
P04, T04 = combustor(pib, P03, T04)
print(P02, T02)
print(P025, T025)
print(P03, T03)
print(P04, T04)