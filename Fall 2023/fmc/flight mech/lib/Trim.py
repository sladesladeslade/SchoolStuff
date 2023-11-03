import numpy as np
import control
from control.matlab import *
from rotations import Euler2Rotation
from scipy.optimize import minimize
from numpy import sin, cos, tan
from UAVdynamics import UAVdynamics
from UAVlinaero import UAVaero
from UAVparams import *


class Trim:
    def __init__(self):
        self.Ts = Ts
        self.fm = UAVaero()
        self.dyn = UAVdynamics()
        gamma = jx*jz - jxz**2
        gamma1 = (jxz*(jx - jy + jz))/gamma
        gamma2 = (jz*(jz - jy) + jxz**2)/gamma
        gamma3 = jz/gamma
        gamma4 = jxz/gamma
        gamma5 = (jz - jx)/jy
        gamma6 = jxz/jy
        gamma7 = ((jx - jy)*jx + jxz**2)/gamma
        gamma8 = jx/gamma
        self.gamma = np.array([gamma, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8])

    def costTrim(self, x, Va, Y, R):
        alpha = x[0]
        beta = x[1]
        phi = x[2]

        xDot = np.zeros((12, 1))
        xDot[2] = -Va * sin(Y)
        xDot[8] = Va / R

        xTrim, uTrim = self.compTrim(x, Va, Y, R)
        deltaa, deltae, deltar, deltat = uTrim.flatten()
        fx, fy, fz = self.fm.forces(xTrim, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
        l, m, n = self.fm.forces(xTrim, alpha, beta, deltaa, deltae, deltar, deltat, Va).flatten()
        f = np.array([[fx], [fy], [fz]])

        stateDot = self.der(xTrim, f, l, m, n)
        J = np.linalg.norm(xDot - stateDot) ** 2

        return J

    def minTrim(self, Va, Y, R):
        x0 = np.array([0, 0, 0])
        res = minimize(lambda x: self.costTrim(x, Va, Y, R), x0, method='nelder-mead',
                       options={'xatol': 1e-8, 'disp': True})
        xTrim, uTrim = self.compTrim(res.x, Va, Y, R)
        return xTrim, uTrim

    def compTrim(self, x, Va, Y, R):
        alpha = x[0]
        beta = x[1]
        phi = x[2]

        u = Va * cos(alpha) * cos(beta)
        v = Va * sin(beta)
        w = Va * sin(alpha) * cos(beta)

        theta = alpha + Y

        p = (-Va / R) * sin(theta)
        q = (Va / R) * sin(phi) * cos(theta)
        r = (Va / R) * cos(phi) * cos(theta)

        xTrim = np.array([0, 0, 0, u, v, w, phi, theta, 0, p, q, r])

        cl = C_L_0 + C_L_alpha * alpha
        cd = C_D_0 + C_D_alpha * alpha

        cx = -cd * cos(alpha) + cl * sin(alpha)
        cxq = -C_D_q * cos(alpha) + C_L_q * sin(alpha)
        cxDe = -C_D_delta_e * cos(alpha) + C_L_delta_e * sin(alpha)

        De = (((jxz * (p ** 2 - r ** 2) + (jx - jz) * p * r) / (
                0.5 * rho * (Va ** 2) * c * S)) - C_m_0 - C_m_alpha * alpha - C_m_q * (
                      (c * q) / (2 * Va))) / C_m_delta_e

        Dt = np.sqrt(((2 * m * (-r * v + q * w + g * sin(theta)) - rho * (Va ** 2) * S *
                       (cx + cxq * ((c * q) / (2 * Va)) + cxDe * De)) /
                      (rho * S_prop * C_prop * k_motor ** 2)) + ((Va ** 2) / (k_motor ** 2)))

        diffMat = np.linalg.inv(np.array([[C_ell_delta_a, C_ell_delta_r],
                                         [C_n_delta_a, C_n_delta_r]]))
        magMat = np.array([[((-self.gamma[1] * p * q + self.gamma[2] * q * r) / (
                0.5 * rho * (Va ** 2) * S * b)) - C_L_0 - C_ell_beta * beta - C_ell_p * (
                                    (b * p) / (2 * Va)) - C_ell_r * ((b * r) / (2 * Va))],
                           [((-self.gamma[7] * p * q + self.gamma[1] * q * r) / (
                                   0.5 * rho * (Va ** 2) * S * b)) - C_n_0 - C_n_beta * beta - C_n_p * (
                                    (b * p) / (2 * Va)) - C_n_r * ((b * r) / (2 * Va))]])
        DaDr = np.matmul(diffMat, magMat)

        Da, Dr = DaDr.flatten()

        uTrim = np.array([[Da], [De], [Dr], [Dt]])

        return xTrim, uTrim

    def der(self, state, f, l, m, n):
        pn, pe, pd, u, v, w, ph, th, ps, p, q, r = state.flatten()

        bv = np.array([[u], [v], [w]])
        iv = Euler2Rotation(ph, th, ps) @ bv
        fvb = 1 / m * f

        pnDot, peDot, pdDot = iv.flatten()

        uvwDot = np.array([[r * v - q * w],
                           [p * w - r * u],
                           [q * u - p * v]]) + fvb
        uDot, vDot, wDot = uvwDot.flatten()

        angVel = np.array([[p], [q], [r]])
        Rgb = np.array([[1, sin(ph) * tan(th), cos(ph) * tan(th)],
                        [0, cos(ph), -sin(ph)],
                        [0, sin(ph) / cos(th), cos(ph) / cos(th)]])
        ptpDot = Rgb @ angVel
        phDot, thDot, psDot = ptpDot.flatten()

        pqrDot = np.array([[self.gamma[1] * p * q - self.gamma[2] * q * r],
                          [self.gamma[5] * p * r - self.gamma[6] * (p ** 2 - r ** 2)],
                          [self.gamma[7] * p * q - self.gamma[1] * q * r]]) + np.array(
            [[self.gamma[3] * l + self.gamma[4] * n],
             [(1 / jy) * m],
             [self.gamma[4] * l + self.gamma[8] * n]])
        pDot, qDot, rDot = pqrDot.flatten()

        xDot = np.array([[pnDot],
                         [peDot],
                         [pdDot],
                         [uDot],
                         [vDot],
                         [wDot],
                         [phDot],
                         [thDot],
                         [psDot],
                         [pDot],
                         [qDot],
                         [rDot],
                         ])
        return xDot
