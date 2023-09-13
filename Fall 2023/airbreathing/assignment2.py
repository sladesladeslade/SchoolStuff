# AEEM4063 Assignment 2
# Slade Brooks
# brooksl@mail.uc.edu

import numpy as np
import matplotlib.pyplot as plt


# hopefully this is always true
gam = 1.4

# ----- Problem 2.1 -----
# set up list of Prs
rs = np.arange(1.01, 70.1, 1)

# calcd polytrop
etainfc = 0.876

# calc etac for each
etacs = np.array([(r**((gam - 1)/gam) - 1)/(r**((gam - 1)/(gam*etainfc)) - 1) for r in rs])

# plot and format
plt.figure()
plt.plot(rs, etacs, "k-")
plt.plot(2., 0.863, "r*")
plt.plot(10., 0.828, "r*")
plt.xlabel("Pressure Ratio")
plt.ylabel("$\eta_c$")
plt.title("Isentropic Efficiency vs Pressure Ratio of the Compressor")
plt.xlim((0, 70))
plt.ylim((0.76, 0.9))
plt.grid()
plt.savefig("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\airbreathing\\images\\a2af1.png")

# show it all
plt.waitforbuttonpress()