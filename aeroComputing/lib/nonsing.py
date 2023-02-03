import numpy as np
import random as r

"""
thx josh <3

generates nonsingular matrices
"""

# Cycles stores N for N x N matrices
cycles = [2, 4, 6, 8, 12]
# Empty arrays size of cycles for data storage and then future plotting
runtime = np.empty((3, len(cycles))).astype('float64')
validation = np.copy(runtime)
# Array creation / loop for testing
det = 0
while det == 0:
    for h in range(len(cycles)):
        n = cycles[h]
        # Empty array size of current N value from cycle, and proper sized b vector
        a = np.empty([n, n])
        b = np.empty(n)
        # Nested for loop that puts a random int value at each location for new matrices
        for i in range(n):
            for j in range(n):
                a[i, j] = r.randint(-15, 15)
            b[i] = r.randint(-15, 15)
        # Flips from row to column vector
        b = b.reshape((n, 1))
        det = np.linalg.det(a)
        if det != 0:
            print(repr(a))
            print(repr(b))