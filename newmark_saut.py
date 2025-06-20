from numba import njit
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt

# Calcule l'excitation extérieure p(t) et la dérivée de phase pour le balayage
@njit(fastmath=True)
def calc_P(T, Vdc, Vac, OMEGA_min, OMEGA_max, t, OMEGA_bal):
    OMEGA_c = (OMEGA_max + OMEGA_min)/2
    delta_OMEGA = (OMEGA_max - OMEGA_min)/2
    theta = OMEGA_c*t + (delta_OMEGA/OMEGA_bal)*np.cos(OMEGA_bal*t)
    dtheta = OMEGA_c - (delta_OMEGA)*(np.sin(OMEGA_bal*t))
    return T * Vdc**2 + 2 * T * Vdc * Vac * np.cos(theta), dtheta

# Calcule la force non linéaire interne
@njit(fastmath=True)
def calc_Fnl(T, Vdc, y):
    return -(3 * T * Vdc**2) * y**2 - (4 * T * Vdc**2) * y**3

# Calcule la dérivée de la force non linéaire pour Newton-Raphson
@njit(fastmath=True)
def calc_dFnl(T, Vdc, y):
    dFY = -6 * T * (Vdc**2) * (y + 2 * y**2)
    dFdY = 0
    return dFY, dFdY

# Schéma de Newmark adapté au balayage de fréquence (OMEGA variable)
@njit(fastmath=True)
def Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal):
    precNR = 1.e-12
    t = t_init
    Y, dY = Y0, dY0
    tt, Yt, dYt = zeros((NT, 1)), zeros((NT, 1)), zeros((NT, 1))
    P, dtheta = calc_P(T, Vdc, Vac, OMEGA_min, OMEGA_max, t, OMEGA_bal)
    Fnl = calc_Fnl(T, Vdc, Y)
    ddY = (P - C * dY - K * Y - Fnl) / M
    l_dtheta = []
    l_dtheta.append(dtheta)
    tt[0], Yt[0], dYt[0] = t, Y, dY
    
    for n in range(1, NT):
        # Intégration temporelle avec correction Newton-Raphson
        t += dt
        Y += dt * dY + (dt**2 / 2) * ddY
        dY += dt * ddY
        P, dtheta = calc_P(T, Vdc, Vac, OMEGA_min, OMEGA_max, t, OMEGA_bal)
        l_dtheta.append(dtheta)
        res = calc_P(T, Vdc, Vac, OMEGA_min, OMEGA_max, t, OMEGA_bal)[0] - M * ddY - C * dY - K * Y - calc_Fnl(T, Vdc, Y)
        normres = np.abs(res / P)

        while normres > precNR:
            dFY, dFdY = calc_dFnl(T, Vdc, Y)
            J = (4 / dt**2) * M + (2 / dt) * (C + dFdY) + K + dFY
            deltaY = res / J
            Y += deltaY
            dY += (2 / dt) * deltaY
            ddY += (4 / dt**2) * deltaY
            res = calc_P(T, Vdc, Vac, OMEGA_min, OMEGA_max, t, OMEGA_bal)[0] - M * ddY - C * dY - K * Y - calc_Fnl(T, Vdc, Y)
            normres = np.abs(deltaY / Y)

        tt[n], Yt[n], dYt[n] = t, Y, dY

    return tt, Yt, dYt, l_dtheta

# Calcule la courbe de réponse pour un balayage de la pulsation d'excitation
@njit(fastmath=True)
def compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGAinit, nb_pts_per, nb_per,OMEGA_bal, tolerance=0.00001):
    npas = int(abs((OMEGA_fin - OMEGA_debut) / dOMEGAinit) + 1)
    OME, AMPL = zeros((npas, 1)), zeros((npas, 1))
    Y0, dY0 = 0.25, 0.25
    k, OMEGA = 0, OMEGA_debut

    while (dOMEGAinit > 0 and OMEGA <= OMEGA_fin) or (dOMEGAinit < 0 and OMEGA >= OMEGA_fin):
        OME[k] = OMEGA
        periode = 2 * np.pi / OMEGA
        dt = periode / nb_pts_per
        NT = nb_per * nb_pts_per
        tt, Yt, dYt, l_theta = Newmark(Y0, dY0, 0, dt, NT, omega0, T, Vdc, Vac, OMEGA, OMEGA, M, C, K,OMEGA_bal)
        AMPL[k] = max(Yt[-3 * nb_pts_per:])
        Y0, dY0 = Yt[-1, 0], dYt[-1, 0]
        OMEGA += dOMEGAinit
        k += 1
    
    val = -1000
    
    return OME[:k], AMPL[:k]

# Trace la courbe de réponse sur un axe matplotlib
# Peut afficher la courbe de référence

def plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data, ax, deltam=0, tracer_data=False):

    if tracer_data:
        ax.plot(OMEGA_data, AMPL_data, color='green', marker='o', label='Courbe de référence')

    ax.plot(OME, AMPL, marker='>', label=f'montée en fréquence pour {deltam}')
    ax.plot(OME2, AMPL2, marker='<', label=f'descente en fréquence pour {deltam}')
    plt.xlabel(r"$\Omega$ pulsation de l'excitation")
    plt.ylabel("Amplitude adimensionnée de la réponse = $max(y(t))$")
    plt.title("Courbe de réponse du résonnateur")
    plt.legend()
    plt.xlim(0.988,0.996)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)


# Initialise les paramètres physiques et adimensionnés du système pour le cas avec balayage
@njit(fastmath=True)
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
    OMEGA_bal = (2*np.pi*250)/omega0
    return T, Vdc, Vac, omega0, M, C, K, OMEGA_bal, d

