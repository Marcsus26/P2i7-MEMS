from numba import njit
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt

# Define functions with numba acceleration
@njit
def calc_P(T, Vdc, Vac, OMEGA, t):
    return T * Vdc**2 + 2 * T * Vdc * Vac * np.cos(OMEGA * t)

@njit
def calc_Fnl(T, Vdc, y):
    return -(3 * T * Vdc**2) * y**2 - (4 * T * Vdc**2) * y**3

@njit
def calc_dFnl(T, Vdc, y):
    dFY = -6 * T * (Vdc**2) * (y + 2 * y**2)
    dFdY = 0
    return dFY, dFdY

@njit
def Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K):
    precNR = 1.e-12
    t = t_init
    Y, dY = Y0, dY0
    tt, Yt, dYt = zeros((NT, 1)), zeros((NT, 1)), zeros((NT, 1))
    
    P = calc_P(T, Vdc, Vac, OMEGA, t)
    Fnl = calc_Fnl(T, Vdc, Y)
    ddY = (P - C * dY - K * Y - Fnl) / M
    tt[0], Yt[0], dYt[0] = t, Y, dY

    for n in range(1, NT):
        t += dt
        iter = 0
        Y += dt * dY + (dt**2 / 2) * ddY
        dY += dt * ddY
        res = calc_P(T, Vdc, Vac, OMEGA, t) - M * ddY - C * dY - K * Y - calc_Fnl(T, Vdc, Y)
        normres = abs(res / P)

        while normres > precNR:
            iter += 1
            dFY, dFdY = calc_dFnl(T, Vdc, Y)
            J = (4 / dt**2) * M + (2 / dt) * (C + dFdY) + K + dFY
            deltaY = res / J
            Y += deltaY
            dY += (2 / dt) * deltaY
            ddY += (4 / dt**2) * deltaY
            res = calc_P(T, Vdc, Vac, OMEGA, t) - M * ddY - C * dY - K * Y - calc_Fnl(T, Vdc, Y)
            normres = abs(deltaY / Y)

        tt[n], Yt[n], dYt[n] = t, Y, dY

    return tt, Yt, dYt

@njit
def compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGAinit, nb_pts_per, nb_per, tolerance=0.00001):
    npas = int(abs((OMEGA_fin - OMEGA_debut) / dOMEGAinit) + 1)
    OME, AMPL = zeros((npas, 1)), zeros((npas, 1))
    Y0, dY0 = 0.25, 0
    k, OMEGA = 0, OMEGA_debut

    while (dOMEGAinit > 0 and OMEGA <= OMEGA_fin) or (dOMEGAinit < 0 and OMEGA >= OMEGA_fin):
        OME[k] = OMEGA
        periode = 2 * np.pi / OMEGA
        dt = periode / nb_pts_per
        NT = nb_per * nb_pts_per
        tt, Yt, dYt = Newmark(Y0, dY0, 0, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)
        AMPL[k] = max(Yt[-3 * nb_pts_per:])
        Y0, dY0 = Yt[-1, 0], dYt[-1, 0]
        OMEGA += dOMEGAinit
        k += 1
    
    val = -1000

    if dOMEGAinit < 0:
        i_max = 0
        for i in range(len(AMPL)):
            if AMPL[i] > AMPL[i_max]:
                i_max = i
        val = OME[i_max]

        npas = int(abs((OMEGA_fin - OMEGA_debut) / dOMEGAinit) + 1 + (0.002) / tolerance)
        OME, AMPL = zeros((npas, 1)), zeros((npas, 1))
        Y0, dY0 = 0.25, 0
        k, OMEGA = 0, OMEGA_debut

        while (dOMEGAinit > 0 and OMEGA <= OMEGA_fin) or (dOMEGAinit < 0 and OMEGA >= OMEGA_fin):
            OME[k] = OMEGA
            if (val - 0.001) < OME[k] and (val + 0.001) > OME[k]:
                dOMEGA = -tolerance
            else:
                dOMEGA = dOMEGAinit
            periode = 2 * np.pi / OMEGA
            dt = periode / nb_pts_per
            NT = nb_per * nb_pts_per
            tt, Yt, dYt = Newmark(Y0, dY0, 0, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)
            AMPL[k] = max(Yt[-3 * nb_pts_per:])
            Y0, dY0 = Yt[-1, 0], dYt[-1, 0]
            OMEGA += dOMEGA
            k += 1

    return OME[:k], AMPL[:k]

@njit
def init_params(deltam=0):
    rho, l, b, h, d = 2500, 250e-6, 40e-6, 1e-6, 0.03e-6
    Vdc, Vac = 5, 5 / 10
    epsilon0 = 8.5e-12
    f0 = 100e6
    omega0 = 2 * np.pi * f0
    Q = 1000
    xi = 1 / Q
    A = l * b
    m = rho * A * h
    T = (epsilon0 * A) / (2 * m * omega0**2 * d**3)
    M = 1 + (deltam / m)
    C = xi
    K = 1 - 2 * T * Vdc**2
    return T, Vdc, Vac, omega0, M, C, K