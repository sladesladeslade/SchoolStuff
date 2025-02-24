# Outputs for Comp. Flow and Norm. Shock Functions
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
from tabulate import tabulate
from compflow import compflow
from normshock import normshock


# range of mach numbers and gammas
Ms = np.arange(0., 4.001, 0.01)
ys = np.array([1.4, 1.33, 1.3])
T_Tts = np.empty((len(ys), len(Ms)))
P_Pts = np.empty((len(ys), len(Ms)))
rho_rhots = np.empty((len(ys), len(Ms)))
A_Astars = np.empty((len(ys), len(Ms)))
MFPsRgcs = np.empty((len(ys), len(Ms)))
mus = np.empty((len(ys), len(Ms)))
vs = np.empty((len(ys), len(Ms)))

# loop through and calculate parameters
for i, y in enumerate(ys):
    for k, M in enumerate(Ms):
        T_Tts[i, k], P_Pts[i, k], rho_rhots[i, k], A_Astars[i, k], MFPsRgcs[i, k], mus[i, k], vs[i, k] = compflow(M, y)

    # output tables
    data = list(zip(np.round(Ms, 2), T_Tts[i, :], P_Pts[i, :], rho_rhots[i, :], A_Astars[i, :], MFPsRgcs[i, :], mus[i, :], vs[i, :]))
    headers = ["M", "T/Tt", "P/Pt", "rho/rhot", "A/A*", "MFP√(R/gc)", "mu", "v"]
    units = ["-", "-", "-", "-", "-", "W-s/m-√T", "deg", "deg"]
    fulldata = [headers] + [units] + data
    colalign = ['center', 'center', 'center', 'center', 'center', 'center', 'center', 'center']
    with open(rf"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\utils\Compflow_y{ys[i]:.2f}_table.txt", "w", encoding="utf-8") as file:
        file.write(f"----------=============== Compressible Flow Tables (gamma={ys[i]:.2f}) ===============----------\n")
        file.write(tabulate(fulldata, headers="firstrow", colalign=colalign, tablefmt="fancy_grid"))  

# range of mach numbers and gammas
Mxs = np.arange(1., 4.001, 0.01)
yxs = np.array([1.4])
Mys = np.empty((len(ys), len(Ms)))
Pty_Ptxs = np.empty((len(ys), len(Ms)))
Py_Pxs = np.empty((len(ys), len(Ms)))
rhoy_rhoxs = np.empty((len(ys), len(Ms)))
Ty_Txs = np.empty((len(ys), len(Ms)))

# loop through and calculate parameters
for i, yx in enumerate(yxs):
    for k, Mx in enumerate(Mxs):
        Mys[i, k], Pty_Ptxs[i, k], Py_Pxs[i, k], rhoy_rhoxs[i, k], Ty_Txs[i, k] = normshock(Mx, yx)

    # output tables
    data = list(zip(np.round(Mxs, 2), Mys[i, :], Pty_Ptxs[i, :], Py_Pxs[i, :], rhoy_rhoxs[i, :], Ty_Txs[i, :]))
    headers = ["Mx", "My", "Pty/Ptx", "Py/Px", "rhoy/rhox", "Ty/Tx"]
    fulldata = [headers] + data
    colalign = ['center', 'center', 'center', 'center', 'center', 'center']
    with open(rf"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\utils\Normshock_y{yxs[i]:.2f}_table.txt", "w", encoding="utf-8") as file:
        file.write(f"----------=============== Normal Shock Tables (gamma={yxs[i]:.2f}) ===============----------\n")
        file.write(tabulate(fulldata, headers="firstrow", colalign=colalign, tablefmt="fancy_grid"))