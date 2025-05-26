import newmark_masse_ajoutee as nkm
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt

@njit
def compute_delta_f(liste_delta_m, OMEGA_debut, OMEGA_fin, dOMEGA, f0):
    OMEGA_ref = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA)
    liste_delta_f = np.zeros(len(liste_delta_m))
    for i in prange(len(liste_delta_m)):
        masse = liste_delta_m[i]
        OMEGA = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA, masse)
        A = OMEGA[0]
        B = OMEGA_ref[0]
        liste_delta_f[i] = f0 * np.abs((A - B)[0])
    return liste_delta_f

def tracer_delta_m_f(liste_delta_m, OMEGA_debut,OMEGA_fin,dOMEGA, f0):
    liste_delta_f = compute_delta_f(np.array(liste_delta_m), OMEGA_debut, OMEGA_fin, dOMEGA, f0)
    plt.plot(liste_delta_m, liste_delta_f, label='Courbe réelle')
    
    # Régression linéaire
    m, b = np.polyfit(liste_delta_m, liste_delta_f, 1)
    y_line = m * liste_delta_m + b
    plt.plot(liste_delta_m, y_line, '--r', label=f'Regression linéaire : y = {m:.2e}x + {b:.2e}')

def tracer_delta_m_f_log(liste_delta_m, OMEGA_debut,OMEGA_fin,dOMEGA, f0):
    liste_delta_f = compute_delta_f(np.array(liste_delta_m), OMEGA_debut, OMEGA_fin, dOMEGA, f0)
    plt.plot(liste_delta_m, liste_delta_f)


# tracer_delta_m_f(np.logspace(-13, -15, 100, base=10), 0.985, 1, 0.0003, 100e6)
# plt.xscale('log')
# plt.xlabel('Masse ajoutée (kg)')
# plt.ylabel('Delta f (Hz)')
# plt.title('Variation de la fréquence de résonance en fonction de la masse ajoutée')
# plt.grid()
# plt.show()