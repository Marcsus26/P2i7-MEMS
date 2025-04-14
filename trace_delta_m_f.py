import newmark_masse_ajoutee as nkm
import numpy as np
import matplotlib.pyplot as plt

OMEGA_debut, OMEGA_fin, dOMEGA = 0.985, 1, 0.0003
f0 = 100

def tracer_delta_m_f(liste_delta_m):
    OMEGA_ref = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA)
    liste_delta_f = []
    for masse in liste_delta_m:
        OMEGA = nkm.courbe_reponse_delta_m(OMEGA_debut, OMEGA_fin, dOMEGA, masse)
        liste_delta_f.append(f0*(abs(OMEGA-OMEGA_ref)))
    plt.plot(liste_delta_m,liste_delta_f)

tracer_delta_m_f([1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15])
plt.show()