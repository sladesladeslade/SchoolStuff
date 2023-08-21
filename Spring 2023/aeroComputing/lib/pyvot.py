# Slade Brooks
# brooksl@mail.uc.edu
# 01.22.2023

# Gauss Pivot, LU Pivot, and Cramer's Rule Module

# this pivots stuff I guess


import numpy as np
from numpy import sign
import sys
import math
import matplotlib.pyplot as plt
import cmath
import random


def err(string):
    ''' err(string).
    Prints 'string' and terminates program.
    '''
    print(string)
    #input('Press return to exit')
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

# cramer

def swapCramer(a, b, i):
    """
    swapCramer(a, b, i).
    Swaps i-th column of matrix [a] with array [b].
    """
    import numpy as np
    ai = a.copy()
    ai[:, i] = np.transpose(b)
    return ai

# guassPivot

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

# LUdecomp

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

# cramer

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

# gaussSeidel

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

# conjGrad

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
    p(x) = c[0] + c[1]x + c[2]x^2 +...+ c[m]x^m
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

# rootsearch

def rootsearch(f, a, b, dx):
    """
    Searches the interval (a,b) in increments dx for
    the bounds (x1,x2) of the smallest root of f(x).
    Returns x1 = x2 = None if no roots were detected.
    """
    x1 = a; f1= f(a)
    x2 = a + dx; f2 = f(x2)
    while sign(f1) == sign(f2):
        if x1 >= b: return None, None
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2)
    else:
        return x1, x2
    
# bisection

def bisection(f, x1, x2, switch=1, tol=1.0e-9):
    """
    Finds a root of f(x) = 0 by bisection.
    The root must be bracketed in (x1,x2).
    Setting switch = 1 returns root = None if
    f(x) increases upon bisection.
    """
    f1 = f(x1)
    if f1 == 0.0: return x1
    f2 = f(x2)
    if f2 == 0.0: return x2
    if sign(f1) == sign(f2):
        err("Root is not bracketed")
    n = int(math.ceil(math.log(abs(x2 - x1)/tol)/math.log(2.0)))

    for i in range(n):
        x3 = 0.5*(x1 + x2); f3 = f(x3)
        if (switch == 1) and (abs(f3) > abs(f1)) \
                        and (abs(f3) > abs(f2)):
            return None
        if f3 == 0.0: return x3
        if sign(f2) != sign(f3): x1 = x3; f1 = f3
        else: x2 = x3; f2 = f3
    return (x1 + x2)/2.0

# newtonRahpson

def newtonRaphson(f, df, a, b, tol=1.0e-9):
    """
    Finds a root of f(x) = 0 by combining the Newton-Raphson
    method with bisection. The root must be bracketed in (a,b).
    Calls user-supplied functions f(x) and its derivative df(x).
    """
    fa = f(a)
    if fa == 0.0: return a
    fb = f(b)
    if fb == 0.0: return b
    if sign(fa) == sign(fb): err("Root is not bracketed")
    x = 0.5*(a+b)
    for i in range(30):
        fx = f(x)
        if fx == 0.0: return x
        if sign(fa) != sign(fx): b=x
        else: a=x
        dfx = df(x)
        try: dx = -fx/dfx
        except ZeroDivisionError: dx = b - a
        x = x + dx
        if (b-x)*(x-a) < 0.0:
            dx = 0.5*(b-a)
            x = a + dx
        if abs(dx) < tol*max(abs(b), 1.0): return x
    print("Too many iterations in Newton-Raphson")

# newtonraphson 2

def newtonRaphson2(f, x, tol=1.0e-9):
    """
    Solves the simultaneous equations f(x) = 0 by
    the Newton-Raphson method using {x} as the initial
    guess. Note that {f} and {x} are vectors.
    """
    def jacobian(f, x):
        h = 1.0e-4
        n = len(x)
        jac = np.zeros((n, n))
        f0 = f(x)
        for i in range(n):
            temp = x[i]
            x[i] = temp + h
            f1 = f(x)
            x[i] = temp
            jac[:,i] = (f1 - f0)/h
        return jac, f0

    for i in range(30):
        jac, f0 = jacobian(f, x)
        if math.sqrt(np.dot(f0,f0)/len(x)) < tol: return x
        dx = gaussPivot(jac, -f0)
        x = x + dx
        if math.sqrt(np.dot(dx, dx)) < tol*max(max(abs(x)), 1.0):
            return x
    print("Too many iterations")

# deflating polys

def deflPoly(a, root):
    """
    Deflates a poly
    """
    n = len(a) - 1
    b = [(0.0 + 0.0j)]*n
    b[n-1] = a[n]
    for i in range(n-2, -1, -1):
        b[i] = a[i+1] + root*b[i+1]
    return b


def evalPolydef(a, x):
    """
    Evaluates the polynomial
    p = a[0] + a[1]*x + a[2]*x^2 +...+ a[n]*x^n
    with its derivatives dp = p' and ddp = p''
    at x.
    """
    n = len(a) - 1
    p = a[n]
    dp = 0.0 + 0.0j
    ddp = 0.0 + 0.0j
    for i in range(1, n+1):
        ddp = ddp*x + 2.0*dp
        dp = dp*x + p
        p = p*x + a[n-i]
    return p, dp, ddp


def polyRoots(a, tol=1.0e-12):
    """
    Uses Laguerre's method to compute all the roots of
    a[0] + a[1]*x + a[2]*x^2 +...+ a[n]*x^n = 0.
    The roots are returned in the array 'roots',
    """
    def laguerre(a, tol):
        x = random.random()
        n = len(a) - 1
        for i in range(30):
            p, dp, ddp = evalPolydef(a, x)
            if abs(p) < tol: return x
            g = dp/p
            h = g*g - ddp/p
            f = cmath.sqrt((n-1)*(n*h - g*g))
            if abs(g+f) > abs(g-f): dx = n/(g+f)
            else: dx = n/(g-f)
            x = x - dx
            if abs(dx) < tol: return x
        print("too many iterations")

    n = len(a) - 1
    roots = np.zeros((n),dtype=complex)
    for i in range(n):
        x = laguerre(a, tol)
        if abs(x.imag) < tol: x=x.real
        roots[i] = x
        a = deflPoly(a, x)
    return roots


# testing
if __name__ == "__main__":

    test = cramer(np.array([[4., 2., -2.],[2., 4., 2.],[-2., 3., 1.]]), np.array([[2.],[16.],[7.]]))
    # print(test)