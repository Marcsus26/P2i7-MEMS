import numpy as np
from newmark import compute_response_curve, init_params

# Plage de masses à tester
masses = [1e-13,1e-14,1e-15,1e-16,1e-17,1e-18,1e-19,1e-20,1e-21,1e-22,1e-23,1e-24,1e-25]  # 1e-12 à 1e-25 inclus

# Paramètres de balayage
OMEGA_debut = 0.97
OMEGA_fin = 1.0
pas_OMEGA = 0.00005
nb_pts_per = 50
nb_per = 500

def save_response_curve_for_mass(mass):
    T, Vdc, Vac, omega0, M, C, K = init_params(mass)
    # Montée en fréquence
    OME_up, AMPL_up = compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_debut, OMEGA_fin, pas_OMEGA, nb_pts_per, nb_per)
    fname_up = f"Courbes de réponse/courbe_reponse_masse_{mass:.0e}_up.txt"
    data_up = np.hstack((OME_up, AMPL_up))
    np.savetxt(fname_up, data_up, delimiter=",", header="OMEGA,AMPL", comments="")
    print(f"Fichier créé : {fname_up}")
    print(f"Courbe montée sauvegardée pour masse {mass:.0e} dans {fname_up}")

    # Descente en fréquence
    OME_down, AMPL_down = compute_response_curve(T, Vdc, Vac, omega0, M, C, K, OMEGA_fin, OMEGA_debut, -pas_OMEGA, nb_pts_per, nb_per)
    fname_down = f"Courbes de réponse/courbe_reponse_masse_{mass:.0e}_down.txt"
    data_down = np.hstack((OME_down, AMPL_down))
    np.savetxt(fname_down, data_down, delimiter=",", header="OMEGA,AMPL", comments="")
    print(f"Fichier créé : {fname_down}")
    print(f"Courbe descente sauvegardée pour masse {mass:.0e} dans {fname_down}")

for mass in masses:
    save_response_curve_for_mass(mass)