# Slade Brooks
# brooksl@mail.uc.edu
# 01.22.2023

# Gauss Pivot, LU Pivot, and Cramer's Rule Module

# this pivots stuff I guess


import numpy as np
import sys
import math


def err(string):
    ''' err(string).
    Prints 'string' and terminates program.
    '''
    print(string)
    input('Press return to exit')
    sys.exit()


def swapRows(v,i,j):
    ''' swapRows(v,i,j).
    Swaps rows i and j of a vector or matrix [v].
    '''
    if len(v.shape) == 1:
        v[i],v[j] = v[j],v[i]
    else:
        v[[i,j],:] = v[[j,i],:]
    return


def swapCols(v,i,j):
    """
    swapCols(v,i,j).
    Swaps columns of matrix [v].
    """
    v[:,[i,j]] = v[:,[j,i]]
    return


def swapCramer(a, b, i):
    """
    swapCramer(a, b, i).
    Swaps i-th column of matrix [a] with array [b].
    """
    import numpy as np
    ai = a.copy()
    ai[:, i] = np.transpose(b)
    return ai


def gaussPivot(a,b,tol=1.0e-12):
    """ 
    x = gaussPivot(a,b,tol=1.0e-12).
    Solves [a]{x} = {b} by Gauss elimination with
    scaled row pivoting
    """
    n = len(b)

    # Set up scale factors
    s = np.zeros(n)
    for i in range(n):
        s[i] = max(np.abs(a[i,:]))

    for k in range(0,n-1):

    # Row interchange, if needed
        p = np.argmax(np.abs(a[k:n,k])/s[k:n]) + k
        if abs(a[p,k]) < tol: err("Matrix is singular")
        if p != k:
            swapRows(b,k,p)
            swapRows(s,k,p)
            swapRows(a,k,p)

    # Elimination
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                b[i] = b[i] - lam*b[k]
    if abs(a[n-1,n-1]) < tol: err("Matrix is singular")

    # Back substitution
    b[n-1] = b[n-1]/a[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    
    return b


def LUdecomp(a, tol=1.0e-9):
    """
    a,seq = LUdecomp(a,tol=1.0e-9).
    LU decomposition of matrix [a] using scaled row pivoting.
    The returned matrix [a] = contains [U] in the upper
    triangle and the nondiagonal terms of [L] in the lower triangle.
    Note that [L][U] is a row-wise permutation of the original [a];
    the permutations are recorded in the vector {seq}.
    """
    n = len(a)
    seq = np.array(range(n))

    # Set up scale factors
    s = np.zeros((n))
    for i in range(n):
        s[i] = max(abs(a[i,:]))

    for k in range(0, n - 1):

        # Row interchange, if needed
        p = np.argmax(np.abs(a[k:n, k])/s[k:n]) + k
        if abs(a[p,k]) < tol: err("Matrix is singular")
        if p != k:
            swapRows(s, k, p)
            swapRows(a, k, p)
            swapRows(seq, k, p)

        # Elimination
        for i in range(k + 1, n):
            if a[i, k] != 0.0:
                lam = a[i, k]/a[k, k]
                a[i, k + 1:n] = a[i, k + 1:n] - lam*a[k, k + 1:n]
                a[i, k] = lam

    return a, seq


def LUsolve(a, b, seq):
    """
    x = LUsolve(a,b,seq).
    Solves [L][U]{x} = {b}, where the matrix [a] = and the
    permutation vector {seq} are returned from LUdecomp.
    """
    n = len(a)

    # Rearrange constant vector; store it in [x]
    x = b.copy()
    for i in range(n):
        x[i] = b[seq[i]]

    # Solution
    for k in range(1, n):
        x[k] = x[k] - np.dot(a[k, 0:k], x[0:k])
    x[n - 1] = x[n - 1]/a[n - 1, n - 1]
    for k in range(n - 2, -1, -1):
        x[k] = (x[k] - np.dot(a[k, k + 1:n], x[k + 1:n]))/a[k,k]

    return x


def cramer(a,b):
    """
    Solves [a]{x}={b} using the cramer method.

    Like seinfeld lol

    :type a: np.array
    :type b: np.array
    :returns: x
    :rtype: np.array
    """
    # get num of iterations
    num_rows, num_cols = a.shape

    # initialize x array
    x = np.zeros([num_rows, 1])
    ai = a

    # loop through
    for i in range(num_cols):
        # reset ai to a
        ai = np.array(a)

        # replace i column of A with b
        ai[:, i] = b[:, 0]

        # calculate xi
        x[i, 0] = (np.linalg.det(ai))/(np.linalg.det(a))

    return x


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


# testing
if __name__ == "__main__":

    test = cramer(np.array([[4., 2., -2.],[2., 4., 2.],[-2., 3., 1.]]), np.array([[2.],[16.],[7.]]))
    print(test)