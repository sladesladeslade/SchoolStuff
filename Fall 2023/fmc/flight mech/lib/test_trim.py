from compute_trim import ComputeTrim
import numpy as np


# set test vars
Va = 35.
Y = np.deg2rad(0.)
R = np.inf

# compute trim and get outputs
xtrim, utrim = ComputeTrim().compute_trim(Va, Y, R)
deltae, deltat, deltaa, deltar = utrim.flatten()

# print out
print("--- States ---")
print(xtrim)
print("--- Trim Conditions ---")
print(f"E: {np.rad2deg(deltae):.2f} deg")
print(f"T: {deltat*100:.2f} %")
print(f"A: {np.rad2deg(deltaa):.2f} deg")
print(f"R: {np.rad2deg(deltar):.2f} deg")