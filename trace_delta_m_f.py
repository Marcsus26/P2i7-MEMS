import newmark_masse_ajoutee as nkm
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt

# Calcule le décalage de fréquence (delta_f) pour une liste de masses ajoutées
@njit
def compute_delta_f(liste_delta_m, OMEGA_debut, OMEGA_fin, dOMEGA, f0):
    OMEGA_ref = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA)  # Référence sans masse ajoutée
    liste_delta_f = np.zeros(len(liste_delta_m))
    for i in prange(len(liste_delta_m)):
        masse = liste_delta_m[i]
        OMEGA = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA, masse)  # Avec masse ajoutée
        A = OMEGA[0]
        B = OMEGA_ref[0]
        liste_delta_f[i] = f0 * np.abs((A - B)[0])  # Calcul du décalage de fréquence
    return liste_delta_f

# Trace la courbe delta_m vs delta_f et la régression linéaire
def tracer_delta_m_f(liste_delta_m, OMEGA_debut,OMEGA_fin,dOMEGA, f0):
    liste_delta_f = compute_delta_f(np.array(liste_delta_m), OMEGA_debut, OMEGA_fin, dOMEGA, f0)
    plt.plot(liste_delta_m, liste_delta_f, label='Courbe réelle')
    
    # Régression linéaire sur les points calculés
    m, b = np.polyfit(liste_delta_m, liste_delta_f, 1)
    y_line = m * liste_delta_m + b
    plt.plot(liste_delta_m, y_line, '--r', label=f"Regression valide pour l'échelle linéaire : y = {m:.2e}x + {b:.2e}")

# Trace la courbe delta_m vs delta_f (version log si besoin)
def tracer_delta_m_f_log(liste_delta_m, OMEGA_debut,OMEGA_fin,dOMEGA, f0):
    liste_delta_f = compute_delta_f(np.array(liste_delta_m), OMEGA_debut, OMEGA_fin, dOMEGA, f0)
    plt.plot(liste_delta_m, liste_delta_f)