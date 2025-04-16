import newmark_masse_ajoutee as nkm
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
OMEGA_debut, OMEGA_fin, dOMEGA = 0.985, 1, 0.0003
f0 = 100

@njit
def compute_delta_f(liste_delta_m, OMEGA_debut, OMEGA_fin, dOMEGA, f0):
    OMEGA_ref = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA)
    liste_delta_f = np.zeros(len(liste_delta_m))
    for i in prange(len(liste_delta_m)):
        masse = liste_delta_m[i]
        OMEGA = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA, masse)
        liste_delta_f[i] = f0 * np.abs(OMEGA - OMEGA_ref)[0]
    return liste_delta_f

def tracer_delta_m_f(liste_delta_m):
    liste_delta_f = compute_delta_f(np.array(liste_delta_m), OMEGA_debut, OMEGA_fin, dOMEGA, f0)
    plt.plot(liste_delta_m, liste_delta_f)

tracer_delta_m_f([1e-13, 1e-14, 1e-15, 1e-16])
plt.show()