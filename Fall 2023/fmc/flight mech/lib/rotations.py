# Rotation Matrices
# Slade Brooks
# brooksl@mail.uc.edu
# stolen

import numpy as np


def Euler2Rotation(phi, theta, psi):
    """
    Converts euler angles to rotation matrix (R_b^i)
    """
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)

    R_roll = np.array([[1, 0, 0],
                       [0, c_phi, -s_phi],
                       [0, s_phi, c_phi]])
    R_pitch = np.array([[c_theta, 0, s_theta],
                        [0, 1, 0],
                        [-s_theta, 0, c_theta]])
    R_yaw = np.array([[c_psi, -s_psi, 0],
                      [s_psi, c_psi, 0],
                      [0, 0, 1]])
    R = R_yaw @ R_pitch @ R_roll

    return R


def Rvb(phi, theta, psi):
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    rvb = np.array([[ct*cpsi, ct*spsi, -st],
                    [sphi*st*cpsi - cphi*spsi, sphi*st*spsi + cphi*cpsi, sphi*ct],
                    [cphi*st*cpsi + sphi*spsi, cphi*st*spsi - sphi*cpsi, cphi*ct]])
    return rvb