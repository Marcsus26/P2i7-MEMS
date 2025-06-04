from newmark_saut import *
import numpy as np
import newmark as nk
from numba import prange

T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params()

# Initialize parameters

OMEGA_min = 0.99
OMEGA_max = 0.992877
OMEGA_c = (OMEGA_max+OMEGA_min)/2

periode = 2*np.pi / OMEGA_c  # Period of excitation and response
nb_pts_per = 50     # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 500    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                # Initial time
NT = nb_per * nb_pts_per
nb_newmark = 750

Y0 = 0
dY0 = 0
liste_dtheta = np.zeros((nb_per*nb_pts_per*nb_newmark))
liste_Yt = np.zeros((nb_per*nb_pts_per*nb_newmark))

for i in range(nb_newmark):
    tt, Yt, dYt, dtheta = Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
    t_init = tt[-1]
    Y0 = Yt[-1]
    dY0 = dYt[-1]
    for j in range(len(Yt)):
        liste_dtheta[j+i*(len(Yt))] = dtheta[j].item() if hasattr(dtheta[j], 'item') else dtheta[j]
        liste_Yt[j+i*(len(Yt))] = Yt[j].item() if hasattr(Yt[j], 'item') else Yt[j]

print('fini')
OME, AMPL = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_0e+00_up.txt')
print('fini')
OME2, AMPL2 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_0e+00_down.txt')
print('fini')

OME_data, AMPL_data = nk.recuperer_courbe_data('Courbes de réponse/courbe_reference.txt')
t_init = tt[-1]
Y0 = Yt[-1]
dY0 = dYt[-1]
T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params(1e-17)

liste_dtheta_new = np.zeros((nb_per*nb_pts_per*nb_newmark*2))
liste_Yt_new = np.zeros((nb_per*nb_pts_per*nb_newmark*2))

for i in range(nb_newmark*2):
    tt_new, Yt_new, dYt_new, dtheta_new = Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
    t_init = tt_new[-1]
    Y0 = Yt_new[-1]
    dY0 = dYt_new[-1]
    for j in range(len(Yt)):
        liste_dtheta_new[j+i*(len(Yt))] = dtheta_new[j].item() if hasattr(dtheta_new[j], 'item') else dtheta_new[j]
        liste_Yt_new[j+i*(len(Yt))] = Yt_new[j].item() if hasattr(Yt_new[j], 'item') else Yt_new[j]

# tt_new, Yt_new, dYt_new, dtheta_new = Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
print('fini')
OME3, AMPL3 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_1e-17_up.txt')
print('fini')
OME4, AMPL4 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_1e-17_down.txt')
print('fini')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(liste_dtheta_new, np.abs(liste_Yt_new), label='Déplacement après ajout de masse (Yt_new)')
ax.plot(liste_dtheta, np.abs(liste_Yt), label='Déplacement sans ajout de masse (Yt)')
nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OME_data, AMPL_data, ax, 0, True)
plot_response_curve(OME3, AMPL3, OME4, AMPL4, 0, 0, ax,1e-17)
ax.set_xlabel('OMEGA adimensionné')
ax.set_ylabel('AMPLITUDE adimensionnée')
ax.grid()
plt.show()
