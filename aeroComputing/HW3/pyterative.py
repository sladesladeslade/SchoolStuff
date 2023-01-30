# Slade Brooks
# brooksl@mail.uc.edu
# 01.30.2023

# Gauss-Seidel and conjGrad module

# this iteratively does things


import numpy as np
import math


def gaussSeidel(iterEqs, x, tol=1.0e-9):
    """
    Gauss-Seidel method for solving [A]{x}={b}.
    The matrix [A] should be sparse. User must supply the
    function iterEqs(x, omega) that returns the improved {x},
    given the current {x} ("omega" is the relaxation factor).
    """
    omega = 1.0
    k = 10
    p = 1
    for i in range(1, 501):
        xOld = x.copy()
        x = iterEqs(x, omega)
        dx = math.sqrt(np.dot(x-xOld, x-xOld))
        if dx < tol: return x, i, omega
        # compute relaxation factor after k+p iterations
        if i== k: dx1 = dx
        if i ==k + p:
            dx2 = dx
            omega = 2.0/(1.0 + math.sqrt(1.0 \
                - (dx2/dx1)**(1.0/p)))
    print("Gauss-Seidel failed to converge")


def conjGrad(Av, x, b, tol=1.0e-9):
    """
    Conjugate gradient method for solving [A]{x}={b}.
    The matrix [A] should be sparse. User must supply
    the function Av(v) that returns the vector [A]{v}.
    """
    n = len(b)
    r = b - Av(x)
    s = r.copy()
    for i in range(n):
        u = Av(s)
        alpha = np.dot(s, r)/np.dot(s, u)
        x = x + alpha*s
        r = b - Av(x)
        if(math.sqrt(np.dot(r, r))) < tol:
            break
        else:
            beta = -np.dot(r, u)/np.dot(s, u)
            s = r + beta*s
    return x, i