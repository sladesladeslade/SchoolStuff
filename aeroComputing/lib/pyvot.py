# Slade Brooks
# brooksl@mail.uc.edu
# 01.22.2023

# Gauss Pivot, LU Pivot, and Cramer's Rule Module

# this pivots stuff I guess


import numpy as np
import sys
import math
import matplotlib.pyplot as plt


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

# LUdecomp3 for use in cubicspline

def LUdecomp3(c, d, e):
    """
    LU decomposition of tridiagonal matrix [c\d\e]. On output
    {c},{d} and {e} are the diagonals of the decomposed matrix.
    """
    n = len(d)
    for k in range(1, n):
        lam = c[k-1]/d[k-1]
        d[k] = d[k] - lam*e[k-1]
        c[k-1] = lam
    return c, d, e


def LUsolve3(c, d, e, b):
    """
    Solves [c\d\e]{x} = {b}, where {c}, {d} and {e} are the
    vectors returned from LUdecomp3.
    """
    n = len(d)
    for k in range(1, n):
        b[k] = b[k] - c[k-1]*b[k-1]
    b[n-1] = b[n-1]/d[n-1]
    for k in range(n-2, -1, -1):
        b[k] = (b[k] - e[k]*b[k+1])/d[k]
    return b

# newtonPoly

def evalPoly(a, xData, x):
    """
    Evaluates Newton's polynomial p at x. The coefficient
    vector {a} can be computed by the function "coeffts".
    """

    n = len(xData) - 1
    p = a[n]
    for k in range(1, n+1):
        p = a[n-k] + (x - xData[n-k])*p
    return p


def coeffts(xData, yData):
    """
    Computes the coefficients of Newton's polynomial.
    """
    m = len(xData)
    a = yData.copy()
    for k in range(1, m):
        a[k:m] = (a[k:m] - a[k-1])/(xData[k:m] - xData[k-1])
    return a

# rational

def rational(xData, yData, x):
    """
    Evaluates the diagonal rational function interpolant p(x)
    that passes through the data points
    """
    m = len(xData)
    r = yData.copy()
    rOld = np.zeros(m)
    for k in range(m-1):
        for i in range(m-k-1):
            if abs(x - (xData[i+k+1])) < 1.0e-9:
                return yData[i+k+1]
            else:
                c1 = r[i+1] - r[i]
                c2 = r[i+1] - rOld[i+1]
                c3 = (x - xData[i])/(x - xData[i+k+1])
                r[i] = r[i+1] + c1/(c3*(1.0 - c1/c2) - 1.0)
                rOld[i+1] = r[i+1]
    return r[0]

# Cubic Spline

def curvatures(xData, yData):
    """
    Returns the curvatures of cubic spline at its knots.
    """
    n = len(xData) - 1
    c = np.zeros(n)
    d = np.ones(n+1)
    e = np.zeros(n)
    k = np.zeros(n+1)
    c[0:n-1] = xData[0:n-1] - xData[1:n]
    d[1:n] = 2.0*(xData[0:n-1] - xData[2:n+1])
    e[1:n] = xData[1:n] - xData[2:n+1]
    k[1:n] = 6.0*(yData[0:n-1] - yData[1:n]) \
                /(xData[0:n-1] - xData[1:n]) \
                -6.0*(yData[1:n] - yData[2:n+1]) \
                /(xData[1:n] - xData[2:n+1])
    LUdecomp3(c, d, e)
    LUsolve3(c, d, e, k)
    return k


def evalSpline(xData, yData, k, x):
    """
    Evaluates cubic spline at x. The curvatures k can be
    computed with the function "curvatures".
    """
    def findSegment(xData, x):
        iLeft = 0
        iRight = len(xData) - 1
        while 1:
            if (iRight - iLeft) <= 1: return iLeft
            i = (iLeft + iRight)//2
            if x < xData[i]: iRight = i
            else: iLeft = i

    i = findSegment(xData, x)
    h = xData[i] - xData[i+1]
    y = ((x - xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0 \
         - ((x - xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0 \
         + (yData[i]*(x - xData[i+1]) - yData[i+1]*(x - xData[i]))/h
    return y

# polyFit

def polyFit(xData, yData, m):
    """
    Returns coefficients of the polynomial
    p(x) = c[0] + c[1]x + c[2]xˆ2 +...+ c[m]xˆm
    that fits the specified data in the least
    squares sense.
    """
    a = np.zeros((m+1, m+1))
    b = np.zeros(m+1)
    s = np.zeros(2*m+1)
    for i in range(len(xData)):
        temp = yData[i]
        for j in range(m+1):
            b[j] = b[j] + temp
            temp = temp*xData[i]
        temp = 1.0
        for j in range(2*m+1):
            s[j] = s[j] + temp
            temp = temp*xData[i]
    for i in range(m+1):
        for j in range(m+1):
            a[i, j] = s[i+j]
    return gaussPivot(a, b)


def stdDev(c, xData, yData):
    """
    Computes the std. deviation between p(x)
    and the data.
    """
    def evalPoly(c, x):
        m = len(c) - 1
        p = c[m]
        for j in range(m):
            p = p*x + c[m-j-1]
        return p

    n = len(xData) - 1
    m = len(c) - 1
    sigma = 0.0
    for i in range(n+1):
        p = evalPoly(c, xData[i])
        sigma = sigma + (yData[i] - p)**2
    sigma = math.sqrt(sigma/(n - m))
    return sigma

# plotPoly

def plotPoly(xData, yData, coeff, xlab="x", ylab="y"):
    """
    Plots data points and the fitting
    polynomial defined by its coefficient
    array coeff = [a0, a1. ...]
    xlab and ylab are optional axis labels
    """
    m = len(coeff)
    x1 = min(xData)
    x2 = max(xData)
    dx = (x2 - x1)/20.0
    x = np.arange(x1, x2 + dx/10.0, dx)
    y = np.zeros((len(x)))*1.0
    for i in range(m):
        y = y + coeff[i]*x**i
    plt.plot(xData, yData, "o", x, y, "-")
    plt.xlabel(xlab); plt.ylabel(ylab)
    plt.grid (True)
    plt.show()


# testing
if __name__ == "__main__":

    test = cramer(np.array([[4., 2., -2.],[2., 4., 2.],[-2., 3., 1.]]), np.array([[2.],[16.],[7.]]))
    print(test)