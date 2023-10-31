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
C0 = M0*np.sqrt(1.4*287*T0)     # m/s

# ---------- Given Vals ----------
yc = 1.4            # spec heat ratio cold
yh = 1.333          # spec heat ratio hot
T04 = 1560          # turbine inlet temp, K
hfuel = 43100       # fuel heat thing whatever, kJ/kg
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
cph = 1.148         # cp of hot scary hot air, kJ/kgK
cpc = 1.005         # cp of air, kJ/kgK
R = 287.            # R of air, kJ/kgK
TR = 70206.7        # thrust required, N


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
    P013 = P02*pif
    
    # solve for n-1/n based on polytropic eff
    n1n = (1/etainff)*(yc - 1)/yc
    
    # solve for T from relation with n-1/n
    T013 = T02*pif**n1n
    
    return P013, T013


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
def turbine(etainft, P04, T04, T03, T02, T013, BPR, hfuel, yh, cph, cpc, etab, etam):
    # solve for mdot thru compressor in terms of total
    mdotc = 1/(BPR + 1)
    
    # solve for f actual
    f = (cph*T04 - cpc*T03)/(etab*(hfuel - cph*T04))
    
    # solve for mdot thru hot section in terms of total
    mdoth = mdotc*(1 + f)
    
    # solve for T05
    T05 = T04 - (cpc*(T013 - T02) + mdotc*cpc*(T03 - T013))/(etam*mdoth*cph)
    
    # solve for m-1/m based on polytropic eff
    m1m = etainft*(yh - 1)/yh
    
    # solve for P05 from m-1/m but inverse because of how it works :)
    P05 = P04*(T05/T04)**(1/m1m)
    
    return P05, T05, f


# ---------- Core Nozzle ----------
def coreNozzle(etaj, P05, T05, P0, yh, R):
    # assume that nodzle is chonked and fully expanded
    P9 = P0
    
    # solve for T9
    T9 = T05 - etaj*T05*(1 - (P9/P05)**((yh - 1)/yh))
    
    # adiabatic nozzle
    T09 = T05
    
    # solve for smach number
    M9 = np.sqrt((T09/T9 - 1)/((yh - 1)/2))
    
    # solve for velocity
    C9 = M9*np.sqrt(yh*R*T9)
    
    return C9


# ---------- Fan Nozzle ----------
def fanNozzle(etaj, P013, T013, P0, yc, R):
    # assume choked/fully expanded
    P19 = P0
    
    # solve for T19
    T19 = T013 - etaj*T013*(1 - (P19/P013)**((yc - 1)/yc))
    
    # adiabatic nozzle
    T019 = T013
    
    # solve for M at exit
    M19 = np.sqrt((T019/T19 - 1)/((yc - 1)/2))
    
    # solve for v at exit
    C19 = M19*np.sqrt(yc*R*T19)
    
    return C19


# ---------- Finish Solving Engine ----------
def solve(TR, BPR, C9, C19, C0, f, idealFlag):
    # first solve for mdot total from thrust req
    mdotT = TR/(1/(BPR + 1)*(C9 - C0) + BPR/(BPR + 1)*(C19 - C0))
    
    # solve for other mdots
    mdotc = BPR/(BPR + 1)*mdotT
    mdoth = 1/(BPR + 1)*mdotT
    mdotf = f/(BPR + 1)*mdotT
    if idealFlag == False:
        mdoth += mdotf
    else:
        mdoth = mdoth
    
    # solve for fan diameter
    A2 = mdotT/(P0*100000/(R*T0)*C0)
    D = np.sqrt(4*A2/np.pi)
    
    return mdotT, D, mdotc, mdoth, mdotf


# ---------- Quick Solve ----------
def quickSolve(P0, M0, etai, yc, pif, etainff, pic, etainfc, pib, etainft, BPR, hfuel, yh, cph, cpc, etab, etam, etaj,
               R, TR, C0, T04, idealFlag):
    P02, T02 = inlet(P0, M0, etai, yc)
    P013, T013 = fan(pif, etainff, P02, T02, yc)
    P03, T03 = compressor(pic, etainfc, P013, T013, yc)
    P04, T04 = combustor(pib, P03, T04)
    P05, T05, f = turbine(etainft, P04, T04, T03, T02, T013, BPR, hfuel, yh, cph, cpc, etab, etam)
    C9 = coreNozzle(etaj, P05, T05, P0, yh, R)
    C19 = fanNozzle(etaj, P013, T013, P0, yc, R)
    mdotT, D, mdotc, mdoth, mdotf = solve(TR, BPR, C9, C19, C0, f, idealFlag)

    # calculate F/mdot
    fmdot = TR/mdotT

    # solve for specific fuel consumption
    TSFC = f*3600/((1 + BPR)*TR/mdotT)
    
    # calculate efficiencies
    etaT = 0.5*(mdoth*C9**2 + mdotc*C19**2 - mdotT*C0**2)/(mdotf*hfuel*1000)
    etaP = C0*(mdotc*(C19 - C0) + mdoth*(C9 - C0))/(0.5*(mdoth*C9**2 + mdotc*C19**2 - mdotT*C0**2))
    etaO = etaT*etaP
    
    return fmdot, TSFC, f, etaT, etaP, etaO
    

# ---------- Real Cycle ----------
# test and verify with hand calcs
fmdot, TSFC, f, etaT, etaP, etaO = quickSolve(P0, M0, etai, yc, pif, etainff, pic, etainfc, pib, etainft, BPR, hfuel,
                                              yh, cph, cpc, etab, etam, etaj, R, TR, C0, T04, False)
# print(fmdot, TSFC, f, etaT, etaP, etaO)


# ---------- Ideal Cycle ----------
# verify with cycle results posted
fmdot, TSFC, f, etaT, etaP, etaO = quickSolve(P0, M0, 1, yc, pif, 1, pic, 1, pib, 1, BPR, hfuel, yh, cph, cpc, 1, 1, 1,
                                              R, TR, C0, T04, True)
# print(fmdot, TSFC, f, etaT, etaP, etaO)


# ---------- Iterating Through Designs ----------
# set up vars to iterate through
BPRs = np.linspace(5, 20., 100)
pifs = np.linspace(1.2, 2., 100)
pics = np.linspace(20, 40., 100)

# iterate thru BPRs
# set up empty lists
fmdots = np.empty((2, 100))
TSFCs = np.empty((2, 100))
fs = np.empty((2, 100))
etaTs = np.empty((2, 100))
etaPs = np.empty((2, 100))
etaOs = np.empty((2, 100))

# loop through and get vals
for i, BPRo in enumerate(BPRs):
    fmdots[0][i], TSFCs[0][i], fs[0][i], etaTs[0][i], etaPs[0][i], etaOs[0][i] = quickSolve(P0, M0, etai, yc, pif,
        etainff, pic, etainfc, pib, etainft, BPRo, hfuel, yh, cph, cpc, etab, etam, etaj, R, TR, C0, T04, False)
    fmdots[1][i], TSFCs[1][i], fs[1][i], etaTs[1][i], etaPs[1][i], etaOs[1][i] = quickSolve(P0, M0, 1, yc, pif, 1, pic,
        1, pib, 1, BPRo, hfuel, yh, cph, cpc, 1, 1, 1, R, TR, C0, T04, False)

# plot specific thrust
plt.figure()
plt.plot(BPRs, fmdots[0], "k", label="Real", linewidth=2)
plt.plot(BPRs, fmdots[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("BPR"); plt.ylabel(r"$\frac{F}{\dot{m}}$ (Ns/kg)")
plt.xlim((5, 20)); plt.ylim((80, 200))
plt.legend(); plt.grid(); plt.title("Specific Thrust vs Bypass Ratio")
plt.savefig("Fall 2023/airbreathing/figs/BPRfmdot.png")

# plot specific fuel consumption
plt.figure()
plt.plot(BPRs, TSFCs[0], "k", label="Real", linewidth=2)
plt.plot(BPRs, TSFCs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("BPR"); plt.ylabel("TSFC (g/kNs)")
plt.xlim((5, 20)); plt.ylim((0.04, 0.075))
plt.legend(); plt.grid(); plt.title("Specific Fuel Consumption vs Bypass Ratio")
plt.savefig("Fall 2023/airbreathing/figs/BPRtsfc.png")

# plot fuel air ratio
plt.figure()
plt.plot(BPRs, fs[0], "k", label="Real", linewidth=2)
plt.plot(BPRs, fs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("BPR"); plt.ylabel("f")
plt.xlim((5, 20)); plt.ylim((0.02, 0.025))
plt.legend(); plt.grid(); plt.title("Fuel-Air Ratio vs Bypass Ratio")
plt.savefig("Fall 2023/airbreathing/figs/BPRf.png")

# plot thermal eff
plt.figure()
plt.plot(BPRs, etaTs[0], "k", label="Real", linewidth=2)
plt.plot(BPRs, etaTs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("BPR"); plt.ylabel(r"$\eta_T$")
plt.xlim((5, 20)); plt.ylim((0.4, 0.65))
plt.legend(); plt.grid(); plt.title("Thermal Efficiency vs Bypass Ratio")
plt.savefig("Fall 2023/airbreathing/figs/BPRetaT.png")

# plot propulsive eff
plt.figure()
plt.plot(BPRs, etaPs[0], "k", label="Real", linewidth=2)
plt.plot(BPRs, etaPs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("BPR"); plt.ylabel(r"$\eta_P$")
plt.xlim((5, 20)); plt.ylim((0.45, 0.85))
plt.legend(); plt.grid(); plt.title("Propulsive Efficiency vs Bypass Ratio")
plt.savefig("Fall 2023/airbreathing/figs/BPRetaP.png")

# plot overall eff
plt.figure()
plt.plot(BPRs, etaOs[0], "k", label="Real", linewidth=2)
plt.plot(BPRs, etaOs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("BPR"); plt.ylabel(r"$\eta_O$")
plt.xlim((5, 20)); plt.ylim((0.25, 0.55))
plt.legend(); plt.grid(); plt.title("Overall Efficiency vs Bypass Ratio")
plt.savefig("Fall 2023/airbreathing/figs/BPRetaO.png")

# loop through and get vals
for i, pifo in enumerate(pifs):
    fmdots[0][i], TSFCs[0][i], fs[0][i], etaTs[0][i], etaPs[0][i], etaOs[0][i] = quickSolve(P0, M0, etai, yc, pifo,
        etainff, pic, etainfc, pib, etainft, BPR, hfuel, yh, cph, cpc, etab, etam, etaj, R, TR, C0, T04, False)
    fmdots[1][i], TSFCs[1][i], fs[1][i], etaTs[1][i], etaPs[1][i], etaOs[1][i] = quickSolve(P0, M0, 1, yc, pifo, 1, pic,
        1, pib, 1, BPR, hfuel, yh, cph, cpc, 1, 1, 1, R, TR, C0, T04, False)

# plot specific thrust
plt.figure()
plt.plot(pifs, fmdots[0], "k", label="Real", linewidth=2)
plt.plot(pifs, fmdots[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("FPR"); plt.ylabel(r"$\frac{F}{\dot{m}}$ (Ns/kg)")
plt.xlim((1.2, 2.0)); plt.ylim((100, 180))
plt.legend(); plt.grid(); plt.title("Specific Thrust vs Fan Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/FPRfmdot.png")

# plot specific fuel consumption
plt.figure()
plt.plot(pifs, TSFCs[0], "k", label="Real", linewidth=2)
plt.plot(pifs, TSFCs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("FPR"); plt.ylabel("TSFC (g/kNs)")
plt.xlim((1.2, 2.0)); plt.ylim((0.04, 0.08))
plt.legend(); plt.grid(); plt.title("Specific Fuel Consumption vs Fan Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/FPRtsfc.png")

# plot fuel air ratio
plt.figure()
plt.plot(pifs, fs[0], "k", label="Real", linewidth=2)
plt.plot(pifs, fs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("FPR"); plt.ylabel("f")
plt.xlim((1.2, 2.0)); plt.ylim((0.02, 0.026))
plt.legend(); plt.grid(); plt.title("Fuel-Air Ratio vs Fan Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/FPRf.png")

# plot thermal eff
plt.figure()
plt.plot(pifs, etaTs[0], "k", label="Real", linewidth=2)
plt.plot(pifs, etaTs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("FPR"); plt.ylabel(r"$\eta_T$")
plt.xlim((1.2, 2.0)); plt.ylim((0.4, 0.7))
plt.legend(); plt.grid(); plt.title("Thermal Efficiency vs Fan Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/FPRetaT.png")

# plot propulsive eff
plt.figure()
plt.plot(pifs, etaPs[0], "k", label="Real", linewidth=2)
plt.plot(pifs, etaPs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("FPR"); plt.ylabel(r"$\eta_P$")
plt.xlim((1.2, 2.0)); plt.ylim((0.45, 0.85))
plt.legend(); plt.grid(); plt.title("Propulsive Efficiency vs Fan Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/FPRetaP.png")

# plot overall eff
plt.figure()
plt.plot(pifs, etaOs[0], "k", label="Real", linewidth=2)
plt.plot(pifs, etaOs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("FPR"); plt.ylabel(r"$\eta_O$")
plt.xlim((1.2, 2.0)); plt.ylim((0.25, 0.5))
plt.legend(); plt.grid(); plt.title("Overall Efficiency vs Fan Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/FPRetaO.png")

# loop through and get vals
for i, pico in enumerate(pics):
    fmdots[0][i], TSFCs[0][i], fs[0][i], etaTs[0][i], etaPs[0][i], etaOs[0][i] = quickSolve(P0, M0, etai, yc, pif,
        etainff, pico, etainfc, pib, etainft, BPR, hfuel, yh, cph, cpc, etab, etam, etaj, R, TR, C0, T04, False)
    fmdots[1][i], TSFCs[1][i], fs[1][i], etaTs[1][i], etaPs[1][i], etaOs[1][i] = quickSolve(P0, M0, 1, yc, pif, 1, pico,
        1, pib, 1, BPR, hfuel, yh, cph, cpc, 1, 1, 1, R, TR, C0, T04, False)

# plot specific thrust
plt.figure()
plt.plot(pics, fmdots[0], "k", label="Real", linewidth=2)
plt.plot(pics, fmdots[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("CPR"); plt.ylabel(r"$\frac{F}{\dot{m}}$ (Ns/kg)")
plt.xlim((20, 40)); plt.ylim((115, 155))
plt.legend(); plt.grid(); plt.title("Specific Thrust vs Compressor Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/CPRfmdot.png")

# plot specific fuel consumption
plt.figure()
plt.plot(pics, TSFCs[0], "k", label="Real", linewidth=2)
plt.plot(pics, TSFCs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("CPR"); plt.ylabel("TSFC (g/kNs)")
plt.xlim((20, 40)); plt.ylim((0.052, 0.066))
plt.legend(); plt.grid(); plt.title("Specific Fuel Consumption vs Compressor Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/CPRtsfc.png")

# plot fuel air ratio
plt.figure()
plt.plot(pics, fs[0], "k", label="Real", linewidth=2)
plt.plot(pics, fs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("CPR"); plt.ylabel("f")
plt.xlim((20, 40)); plt.ylim((0.021, 0.028))
plt.legend(); plt.grid(); plt.title("Fuel-Air Ratio vs Compressor Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/CPRf.png")

# plot thermal eff
plt.figure()
plt.plot(pics, etaTs[0], "k", label="Real", linewidth=2)
plt.plot(pics, etaTs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("CPR"); plt.ylabel(r"$\eta_T$")
plt.xlim((20, 40)); plt.ylim((0.45, 0.65))
plt.legend(); plt.grid(); plt.title("Thermal Efficiency vs Compressor Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/CPRetaT.png")

# plot propulsive eff
plt.figure()
plt.plot(pics, etaPs[0], "k", label="Real", linewidth=2)
plt.plot(pics, etaPs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("CPR"); plt.ylabel(r"$\eta_P$")
plt.xlim((20, 40)); plt.ylim((0.6, 0.8))
plt.legend(); plt.grid(); plt.title("Propulsive Efficiency vs Compressor Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/CPRetaP.png")

# plot overall eff
plt.figure()
plt.plot(pics, etaOs[0], "k", label="Real", linewidth=2)
plt.plot(pics, etaOs[1], "k--", label="Ideal", linewidth=2)
plt.xlabel("CPR"); plt.ylabel(r"$\eta_O$")
plt.xlim((20, 40)); plt.ylim((0.3, 0.4))
plt.legend(); plt.grid(); plt.title("Overall Efficiency vs Compressor Pressure Ratio")
plt.savefig("Fall 2023/airbreathing/figs/CPRetaO.png")

plt.show()