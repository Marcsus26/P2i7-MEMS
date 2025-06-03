from newmark_saut import *
import numpy as np
import newmark as nk
from numba import prange

T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params()

# Initialize parameters

OMEGA_min = 0.988
OMEGA_max = 0.994
OMEGA_c = (OMEGA_max+OMEGA_min)/2

periode = 2*np.pi / OMEGA_c  # Period of excitation and response
nb_pts_per = 500     # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 100000    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                # Initial time
NT = nb_per * nb_pts_per

tt, Yt, dYt, dtheta = Newmark(0, 0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
print('fini')
OME, AMPL = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_0e+00_up.txt')
print('fini')
OME2, AMPL2 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_0e+00_down.txt')
print('fini')

OME_data, AMPL_data = nk.recuperer_courbe_data('Courbes de réponse/courbe_reference.txt')
t_init = tt[-1]
Y0 = Yt[-1]
dY0 = dYt[-1]
T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params(1e-13)
tt_new, Yt_new, dYt_new, dtheta_new = Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
print('fini')
OME3, AMPL3 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_1e-13_up.txt')
print('fini')
OME4, AMPL4 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_1e-13_down.txt')
print('fini')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(dtheta_new, np.abs(Yt_new), label='Déplacement après ajout de masse (Yt_new)')
ax.plot(dtheta, np.abs(Yt), label='Déplacement (Yt)')
nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OME_data, AMPL_data, ax, 0, True)
plot_response_curve(OME3, AMPL3, OME4, AMPL4, 0, 0, ax,1e-13)
ax.set_xlabel('Temps (s)')
ax.set_ylabel('Déplacement')
ax.grid()
plt.show()
