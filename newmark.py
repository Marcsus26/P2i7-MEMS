# Import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import time

# Définition des fonctions de calcul
def calc_P(T, Vdc, Vac, OMEGA, t):
    return T * Vdc**2 + 2 * T * Vdc * Vac * np.cos(OMEGA * t)

def calc_Fnl(T, Vdc, y):
    return -(3 * T * Vdc**2) * y**2 - (4 * T * Vdc**2) * y**3

def calc_dFnl(T, Vdc, y):
    dFY = -6 * T * (Vdc**2) * (y + 2 * y**2)
    dFdY = 0
    return dFY, dFdY

# Fonction Newmark
def Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K):
    precNR = 1.e-9
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

# Initialisation des paramètres
def init_params():
    rho, l, b, h, d = 2500, 250e-6, 40e-6, 1e-6, 0.03e-6
    Vdc, Vac = 5, 5 / 10
    epsilon0 = 8.5e-12
    f0 = 100e6
    omega0 = 2 * np.pi * f0
    Q = 1000
    xi = 1 / Q
    A = l * b
    deltam = 0
    m = rho * A * h
    T = (epsilon0 * A) / (2 * m * omega0**2 * d**3)
    M = 1 + (deltam / m)
    C = xi
    K = 1 - 2 * T * Vdc**2
    return T, Vdc, Vac, omega0, M, C, K

# Calcul de la courbe de réponse
def compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per):
    npas = int(abs((OMEGA_fin - OMEGA_debut) / dOMEGA) + 1)
    OME, AMPL = zeros((npas, 1)), zeros((npas, 1))
    Y0, dY0 = 0.25, 0
    k, OMEGA = 0, OMEGA_debut

    while (dOMEGA > 0 and OMEGA <= OMEGA_fin) or (dOMEGA < 0 and OMEGA >= OMEGA_fin):
        OME[k] = OMEGA
        periode = 2 * np.pi / OMEGA
        dt = periode / nb_pts_per
        NT = nb_per * nb_pts_per
        tt, Yt, dYt = Newmark(Y0, dY0, 0, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)
        AMPL[k] = max(Yt[-3 * nb_pts_per:])
        print(f'ome= {OME[k, 0]:0.5f}  y= {AMPL[k, 0]:0.5g}', end="\r", flush=True)
        Y0, dY0 = Yt[-1, 0], dYt[-1, 0]
        OMEGA += dOMEGA
        k += 1

    return OME[:k], AMPL[:k]

# Affichage des résultats
def plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data):
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot()
    ax.plot(OME, AMPL, color='red', marker='>', label='montée en fréquence')
    ax.plot(OME2, AMPL2, color='blue', marker='<', label='descente en fréquence')
    ax.plot(OMEGA_data, AMPL_data, color='green', marker='o', label='données fichier')
    plt.xlabel("$\Omega$ pulsation de l'excitation")
    plt.ylabel("Amplitude de la réponse = $max(y(t))$")
    plt.title("Courbe de réponse")
    plt.legend()
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.show()

# Main
if __name__ == "__main__":
    T, Vdc, Vac, omega0, M, C, K = init_params()
    OMEGA_debut, OMEGA_fin, dOMEGA = 0.9, 1.10, 0.005
    nb_pts_per, nb_per = 50, 500

    OME, AMPL = compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per)

    # Descente en fréquence
    OME2, AMPL2 = compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -dOMEGA, nb_pts_per, nb_per)

    # Chargement des données de la courbe de réponse
    data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
    OMEGA_data, AMPL_data = data[:, 0], data[:, 1]

    # Affichage
    plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data)