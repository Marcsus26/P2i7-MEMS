# Import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import time

import newmark as nk
from numba import njit

# OMEGA_debut, OMEGA_fin, dOMEGA = 0.985, 1, 0.0003


# T, Vdc, Vac, omega0, M, C, K = nk.init_params()
# nb_pts_per, nb_per = 50, 500
# fig,ax = plt.subplots(1,1,figsize=(10,4))

# OME, AMPL = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per)

    
# # Descente en fréquence
# OME2, AMPL2 = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -dOMEGA, nb_pts_per, nb_per)

# val1 = OME2[np.argmax(AMPL2)]
# Chargement des données de la courbe de réponse
# data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
# OMEGA_data, AMPL_data = data[:, 0], data[:, 1]

# # Affichage
# nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data,ax)
@njit
def courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA, deltam=0, montee=False):
    T, Vdc, Vac, omega0, M, C, K = nk.init_params(deltam)
    nb_pts_per, nb_per = 50, 500
    if montee:
        OME, AMPL = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per)

    # Descente en fréquence
    OME2, AMPL2 = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -dOMEGA, nb_pts_per, nb_per)

    OMEGA_MAX = OME2[np.argmax(AMPL2)]
    # Chargement des données de la courbe de réponse
    # data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
    # OMEGA_data, AMPL_data = data[:, 0], data[:, 1]

    # Affichage
    # nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data,ax,deltam, True)
    return OMEGA_MAX, OME, AMPL, OME2, AMPL2

def affichage_courbe_rep():
    Omega_max, OME, AMPL, OME2, AMPL2 = courbe_reponse_delta_m(0.985, 1, 0.0003, 0, True)
    Omega_max, OME13, AMPL13, OME213, AMPL213 = courbe_reponse_delta_m(0.985, 1, 0.0003, 1e-13, True)
    Omega_max, OME14, AMPL14, OME214, AMPL214 = courbe_reponse_delta_m(0.985, 1, 0.0003, 1e-14, True)
    data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
    OMEGA_data, AMPL_data = data[:, 0], data[:, 1]
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(OME, AMPL, label='Montee')
    ax.plot(OME2, AMPL2, label='Descente')
    # ax.plot(OME13, AMPL13, label='Montee avec deltam=1e-13')
    # ax.plot(OME213, AMPL213, label='Descente avec deltam=1e-13')
    ax.plot(OME14, AMPL14, label='Montee avec deltam=1e-14')
    ax.plot(OME214, AMPL214, label='Descente avec deltam=1e-14')
    ax.plot(OMEGA_data, AMPL_data, label='Données expérimentales', linestyle='--')
    ax.legend()
    ax.set_xlabel('OMEGA')
    ax.set_ylabel('AMPL')
    ax.set_title('Courbe de réponse avec et sans masse ajoutée')
    ax.grid()



affichage_courbe_rep()
plt.show()
