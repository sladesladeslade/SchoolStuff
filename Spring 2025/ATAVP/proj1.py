# ATAVP Proj. 1 Data Plotting
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import matplotlib.pyplot as plt


# read in data function
def read_data_file(filename):
    with open(filename, 'r') as file:
        header = file.readline().strip().split(', ')
        data = np.loadtxt(file, delimiter=',')
    
    data_dict = {header[i]: data[:, i] for i in range(len(header))}
    return data_dict


# process files into dicts
fpicsweep = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\picsweep.txt"
fpicsweepr = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\picsweepr.txt"
falphasweep = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\alphasweep.txt"
falphasweepr = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\alphasweepr.txt"
ftt4sweep = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\tt4sweep.txt"
ftt4sweepr = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\tt4sweepr.txt"
fperf = r"C:\Users\spbro\Documents\SchoolStuff\Spring 2025\ATAVP\PARA\perf.txt"
picsweep = read_data_file(fpicsweep)
picsweepr = read_data_file(fpicsweepr)
alphasweep = read_data_file(falphasweep)
alphasweepr = read_data_file(falphasweepr)
tt4sweep = read_data_file(ftt4sweep)
tt4sweepr = read_data_file(ftt4sweepr)
perf = read_data_file(fperf)

# pic sweep plots
plt.figure("Pic Spec. T")
plt.plot(picsweep["Pic"], picsweep["F/m0"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(picsweepr["Pic"], picsweepr["F/m0"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$\pi_c$"); plt.ylabel("Specific Thrust (N-s/kg)")
plt.xlim(picsweep["Pic"][0], picsweep["Pic"][-1]); plt.ylim(188, 202)
plt.legend(); plt.title(r"Engine 1 Specific Thrust for Varying $\pi_c$")
plt.grid(); plt.tight_layout()
plt.savefig("results/piT.png", dpi=300)
plt.figure("Pic SFC")
plt.plot(picsweep["Pic"], picsweep["S"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(picsweepr["Pic"], picsweepr["S"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$\pi_c$"); plt.ylabel("Specific Fuel Consumption (g/kN-s)")
plt.xlim(picsweep["Pic"][0], picsweep["Pic"][-1]); plt.ylim(12, 24)
plt.legend(); plt.title(r"Engine 1 Specific Fuel Consumption for Varying $\pi_c$")
plt.grid(); plt.tight_layout()
plt.savefig("results/piS.png", dpi=300)
plt.figure("Pic f")
plt.plot(picsweep["Pic"], picsweep["f"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(picsweepr["Pic"], picsweepr["f"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$\pi_c$"); plt.ylabel("Fuel-Air Ratio")
plt.xlim(picsweep["Pic"][0], picsweep["Pic"][-1]); plt.ylim(0.02, 0.045)
plt.legend(); plt.title(r"Engine 1 Fuel-Air Ratio for Varying $\pi_c$")
plt.grid(); plt.tight_layout()
plt.savefig("results/pif.png", dpi=300)
plt.figure("Pic etas")
plt.plot(picsweep["Pic"], picsweep["OvAllEff"], color="tab:blue", linestyle="--", label=r"Ideal $\eta_O$")
plt.plot(picsweepr["Pic"], picsweepr["OvAllEff"], color="tab:blue", linestyle="-", label=r"Real $\eta_O$")
plt.plot(picsweep["Pic"], picsweep["ThermEff"], color="tab:orange", linestyle="--", label=r"Ideal $\eta_T$")
plt.plot(picsweepr["Pic"], picsweepr["ThermEff"], color="tab:orange", linestyle="-", label=r"Real $\eta_T$")
plt.plot(picsweep["Pic"], picsweep["PropEff"], color="tab:green", linestyle="--", label=r"Ideal $\eta_P$")
plt.plot(picsweepr["Pic"], picsweepr["PropEff"], color="tab:green", linestyle="-", label=r"Real $\eta_P$")
plt.xlabel(r"$\pi_c$"); plt.ylabel("Efficiency (%)")
plt.xlim(picsweep["Pic"][0], picsweep["Pic"][-1]); plt.ylim(0, 100)
plt.legend(loc="lower right"); plt.title(r"Engine 1 Efficiencies for Varying $\pi_c$")
plt.grid(); plt.tight_layout()
plt.savefig("results/pie.png", dpi=300)

# alpha sweep plots
plt.figure("alpha Spec. T")
plt.plot(alphasweep["Alpha"], alphasweep["F/m0"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(alphasweepr["Alpha"], alphasweepr["F/m0"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$\alpha$"); plt.ylabel("Specific Thrust (N-s/kg)")
plt.xlim(alphasweep["Alpha"][0], alphasweep["Alpha"][-1]); plt.ylim(150, 275)
plt.legend(); plt.title(r"Engine 1 Specific Thrust for Varying $\alpha$")
plt.grid(); plt.tight_layout()
plt.savefig("results/alphaT.png", dpi=300)
plt.figure("alpha SFC")
plt.plot(alphasweep["Alpha"], alphasweep["S"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(alphasweepr["Alpha"], alphasweepr["S"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$\alpha$"); plt.ylabel("Specific Fuel Consumption (g/kN-s)")
plt.xlim(alphasweep["Alpha"][0], alphasweep["Alpha"][-1]); plt.ylim(12, 24)
plt.legend(); plt.title(r"Engine 1 Specific Fuel Consumption for Varying $\alpha$")
plt.grid(); plt.tight_layout()
plt.savefig("results/alphaS.png", dpi=300)
plt.figure("alpha f")
plt.plot(alphasweep["Alpha"], alphasweep["f"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(alphasweepr["Alpha"], alphasweepr["f"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$\alpha$"); plt.ylabel("Fuel-Air Ratio")
plt.xlim(alphasweep["Alpha"][0], alphasweep["Alpha"][-1]); plt.ylim(0.02, 0.045)
plt.legend(); plt.title(r"Engine 1 Fuel-Air Ratio for Varying $\alpha$")
plt.grid(); plt.tight_layout()
plt.savefig("results/alphaf.png", dpi=300)
plt.figure("alpha etas")
plt.plot(alphasweep["Alpha"], alphasweep["OvAllEff"], color="tab:blue", linestyle="--", label=r"Ideal $\eta_O$")
plt.plot(alphasweepr["Alpha"], alphasweepr["OvAllEff"], color="tab:blue", linestyle="-", label=r"Real $\eta_O$")
plt.plot(alphasweep["Alpha"], alphasweep["ThermEff"], color="tab:orange", linestyle="--", label=r"Ideal $\eta_T$")
plt.plot(alphasweepr["Alpha"], alphasweepr["ThermEff"], color="tab:orange", linestyle="-", label=r"Real $\eta_T$")
plt.plot(alphasweep["Alpha"], alphasweep["PropEff"], color="tab:green", linestyle="--", label=r"Ideal $\eta_P$")
plt.plot(alphasweepr["Alpha"], alphasweepr["PropEff"], color="tab:green", linestyle="-", label=r"Real $\eta_P$")
plt.xlabel(r"$\alpha$"); plt.ylabel("Efficiency (%)")
plt.xlim(alphasweep["Alpha"][0], alphasweep["Alpha"][-1]); plt.ylim(0, 100)
plt.legend(loc="lower right"); plt.title(r"Engine 1 Efficiencies for Varying $\alpha$")
plt.grid(); plt.tight_layout()
plt.savefig("results/alphae.png", dpi=300)

# Tt4 sweep plots
plt.figure("tt4 Spec. T")
plt.plot(tt4sweep["Tt4R"], tt4sweep["F/m0"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(tt4sweepr["Tt4R"], tt4sweepr["F/m0"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$T_{t4}$ (K)"); plt.ylabel("Specific Thrust (N-s/kg)")
plt.xlim(tt4sweep["Tt4R"][0], tt4sweep["Tt4R"][-1]); plt.ylim(160, 260)
plt.legend(); plt.title(r"Engine 1 Specific Thrust for Varying $T_{t4}$")
plt.grid(); plt.tight_layout()
plt.savefig("results/tt4T.png", dpi=300)
plt.figure("tt4 SFC")
plt.plot(tt4sweep["Tt4R"], tt4sweep["S"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(tt4sweepr["Tt4R"], tt4sweepr["S"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$T_{t4}$ (K)"); plt.ylabel("Specific Fuel Consumption (g/kN-s)")
plt.xlim(tt4sweep["Tt4R"][0], tt4sweep["Tt4R"][-1]); plt.ylim(10, 30)
plt.legend(); plt.title(r"Engine 1 Specific Fuel Consumption for Varying $T_{t4}$")
plt.grid(); plt.tight_layout()
plt.savefig("results/tt4S.png", dpi=300)
plt.figure("tt4 f")
plt.plot(tt4sweep["Tt4R"], tt4sweep["f"], color="tab:blue", linestyle="--", label="Ideal")
plt.plot(tt4sweepr["Tt4R"], tt4sweepr["f"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$T_{t4}$ (K)"); plt.ylabel("Fuel-Air Ratio")
plt.xlim(tt4sweep["Tt4R"][0], tt4sweep["Tt4R"][-1]); plt.ylim(0.02, 0.06)
plt.legend(); plt.title(r"Engine 1 Fuel-Air Ratio for Varying $T_{t4}$")
plt.grid(); plt.tight_layout()
plt.savefig("results/tt4fpng", dpi=300)
plt.figure("tt4 etas")
plt.plot(tt4sweep["Tt4R"], tt4sweep["OvAllEff"], color="tab:blue", linestyle="--", label=r"Ideal $\eta_O$")
plt.plot(tt4sweepr["Tt4R"], tt4sweepr["OvAllEff"], color="tab:blue", linestyle="-", label=r"Real $\eta_O$")
plt.plot(tt4sweep["Tt4R"], tt4sweep["ThermEff"], color="tab:orange", linestyle="--", label=r"Ideal $\eta_T$")
plt.plot(tt4sweepr["Tt4R"], tt4sweepr["ThermEff"], color="tab:orange", linestyle="-", label=r"Real $\eta_T$")
plt.plot(tt4sweep["Tt4R"], tt4sweep["PropEff"], color="tab:green", linestyle="--", label=r"Ideal $\eta_P$")
plt.plot(tt4sweepr["Tt4R"], tt4sweepr["PropEff"], color="tab:green", linestyle="-", label=r"Real $\eta_P$")
plt.xlabel(r"$T_{t4}$ (K)"); plt.ylabel("Efficiency (%)")
plt.xlim(tt4sweep["Tt4R"][0], tt4sweep["Tt4R"][-1]); plt.ylim(0, 100)
plt.legend(loc="lower right"); plt.title(r"Engine 1 Efficiencies for Varying $T_{t4}$")
plt.grid(); plt.tight_layout()
plt.savefig("results/tt4e.png", dpi=300)

# perf stuff
perfspecT = perf["Thrust"]/perf["mdot"]
perff = perf["S"]/1000./1000.*(1 + perf["Alpha"])*perfspecT

# M0 sweep plots
plt.figure("M0 Spec. T")
plt.plot(perf["M0"], perfspecT, color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$M_0$"); plt.ylabel("Specific Thrust (N-s/kg)")
plt.xlim(perf["M0"][0], perf["M0"][-1]); plt.ylim(160, 300)
plt.title(r"Engine 1 Specific Thrust for Varying $M_0$")
plt.grid(); plt.tight_layout()
plt.savefig("results/M0t.png", dpi=300)
plt.figure("M0 SFC")
plt.plot(perf["M0"], perf["S"], color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$M_0$"); plt.ylabel("Specific Fuel Consumption (g/kN-s)")
plt.xlim(perf["M0"][0], perf["M0"][-1]); plt.ylim(10, 20)
plt.title(r"Engine 1 Specific Fuel Consumption for Varying $M_0$")
plt.grid(); plt.tight_layout()
plt.savefig("results/M0S.png", dpi=300)
plt.figure("M0 f")
plt.plot(perf["M0"], perff, color="tab:blue", linestyle="-", label="Real")
plt.xlabel(r"$M_0$"); plt.ylabel("Fuel-Air Ratio")
plt.xlim(perf["M0"][0], perf["M0"][-1]); plt.ylim(0.026, 0.032)
plt.title(r"Engine 1 Fuel-Air Ratio for Varying $M_0$")
plt.grid(); plt.tight_layout()
plt.savefig("results/M0f.png", dpi=300)
plt.figure("M0 rats")
plt.plot(perf["M0"], perf["Pic"]/10, color="tab:blue", linestyle="-", label=r"$\pi_c/10$")
plt.plot(perf["M0"], perf["pif"], color="tab:orange", linestyle="-", label=r"$\pi_f$")
plt.plot(perf["M0"], perf["Alpha"], color="tab:green", linestyle="-", label=r"$\alpha$")
plt.xlabel(r"$M_0$"); plt.ylabel("Ratio")
plt.xlim(perf["M0"][0], perf["M0"][-1]); plt.ylim(0, 10)
plt.title(r"Engine 1 $\pi$s for Varying $M_0$")
plt.legend()
plt.grid(); plt.tight_layout()
plt.savefig("results/M0r.png", dpi=300)

plt.show()