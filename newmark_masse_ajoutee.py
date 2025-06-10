# Import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import time

import newmark as nk
from numba import njit

# Calcule la courbe de réponse (amplitude vs pulsation) pour un MEMS avec masse ajoutée
@njit(fastmath=True)
def courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA, deltam=0, montee=False, data = False):
    T, Vdc, Vac, omega0, M, C, K, d = nk.init_params(deltam)  # Initialisation des paramètres avec masse ajoutée
    nb_pts_per, nb_per = 50, 500  # Paramètres numériques
    if montee:
        OME, AMPL = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per)  # Courbe en montée

    # Descente en fréquence
    OME2, AMPL2 = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -dOMEGA, nb_pts_per, nb_per)  # Courbe en descente

    OMEGA_MAX = OME2[np.argmax(AMPL2)]  # Extraction du maximum d'amplitude
    return OMEGA_MAX, OME, AMPL, OME2, AMPL2

