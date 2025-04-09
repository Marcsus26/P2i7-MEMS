# Import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import time

import newmark as nk

T, Vdc, Vac, omega0, M, C, K = nk.init_params()
OMEGA_debut, OMEGA_fin, dOMEGA = 0.9, 1.10, 0.005
nb_pts_per, nb_per = 50, 500
fig,ax = plt.subplots(1,1,figsize=(10,4))

OME, AMPL = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per)

# Descente en fréquence
OME2, AMPL2 = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -dOMEGA, nb_pts_per, nb_per)

# Chargement des données de la courbe de réponse
data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
OMEGA_data, AMPL_data = data[:, 0], data[:, 1]

# Affichage
nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data,ax)
deltam = 1e-13
T, Vdc, Vac, omega0, M, C, K = nk.init_params(deltam)
OMEGA_debut, OMEGA_fin, dOMEGA = 0.9, 1.10, 0.005
nb_pts_per, nb_per = 50, 500

OME, AMPL = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, dOMEGA, nb_pts_per, nb_per)

# Descente en fréquence
OME2, AMPL2 = nk.compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -dOMEGA, nb_pts_per, nb_per)

# Chargement des données de la courbe de réponse
data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
OMEGA_data, AMPL_data = data[:, 0], data[:, 1]

# Affichage
nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data,ax,deltam)

plt.show()
