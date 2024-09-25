import numpy as np
import matplotlib.pyplot as plt


# constants
R = 200.
u = -102.759
v = 95.925
a = -3.3617e-4
b = -6.443e-5

# s ranges
npts = 100
s1 = np.linspace(0, 3*np.pi*R/2, npts)
s2 = np.linspace(0, R, npts)
s3 = np.linspace(0, np.pi*R, npts)
s4 = np.linspace(0, np.pi*R/2, npts)

# f and g eqtns
def f12(s): return -u*s - (R**2)*np.sin(s/R)
def f23(s): return f12(s1[-1]) - u*s
def f34(s): return f23(s2[-1]) - u*s + 4*(R**2)*np.cos(s/(2*R))-4*R**2
def f41(s): return f34(s3[-1]) - 3*R*s/2 - u*s - (R**2)/4*np.sin(2*s/R)
def g12(s): return -v*s + (R**2)*(np.cos(s/R) - 1)
def g23(s): return g12(s1[-1]) + R*s - v*s + (s**2)/2
def g34(s): return g23(s2[-1]) -v*s + 4*(R**2)*np.sin(s/(2*R))
def g41(s): return g34(s3[-1]) -v*s - (R**2)/4*(np.cos(2*s/R) - 1)

# get max and min
q12 = a*f12(s1) + b*g12(s1)
max = np.max(q12)
maxi = np.argmax(q12)
smax = s1[maxi]
print(f"Max q(s): {max:.3f} @ {smax:.1f} mm")
q34 = a*f34(s3) + b*g34(s3)
min = np.min(q34)
mini = np.argmin(q34)
smin = s3[mini] + s1[-1] + s2[-1]
print(f"Min q(s): {min:.3f} @ {smin:.1f} mm")

# plot it
plt.figure()
plt.plot(s1, a*f12(s1) + b*g12(s1), label="$q_{s_1}$", linewidth=2)
plt.plot(s1[-1] + s2, a*f23(s2) + b*g23(s2), label="$q_{s_2}$", linewidth=2)
plt.plot(s1[-1] + s2[-1] + s3, a*f34(s3) + b*g34(s3), label="$q_{s_3}$", linewidth=2)
plt.plot(s1[-1] + s2[-1] + s3[-1] + s4, a*f41(s4) + b*g41(s4), label="$q_{s_4}$", linewidth=2)
plt.ylim(-60, 10); plt.xlim(0, 2500)
plt.ylabel("Shear Flow $q(s)$ (N/mm)"); plt.xlabel("Arc Length $s$ (mm)")
plt.grid(); plt.legend(loc="lower right")
plt.show()