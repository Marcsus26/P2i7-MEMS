from newmark_saut import *
import numpy as np
import newmark as nk
from numba import prange, njit


@njit
def saut_ampl(T, Vdc, Vac, omega0, M, C, K, OMEGA_bal, OMEGA_min, OMEGA_max, nb_pts_per, dt, nb_per, t_init, NT, nb_newmark, delta_m):
    Y0 = np.array([0.0])
    dY0 = np.array([0.0])
    liste_dtheta = np.zeros((nb_per*nb_pts_per*nb_newmark, 1))
    liste_Yt = np.zeros((nb_per*nb_pts_per*nb_newmark, 1))

    for i in range(nb_newmark):
        tt, Yt, dYt, dtheta = Newmark(Y0[0], dY0[0], t_init[0], dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
        t_init = tt[-1]
        Y0 = Yt[-1]
        dY0 = dYt[-1]
        for j in prange(len(Yt)):
            liste_dtheta[j+i*(len(Yt)),0] = dtheta[j]
            liste_Yt[j+i*(len(Yt)),0] = Yt[j,0]


    t_init = tt[-1]
    Y0 = Yt[-1]
    dY0 = dYt[-1]
    T, Vdc, Vac, omega0, M, C, K, OMEGA_bal, d = init_params(delta_m)

    liste_dtheta_new = np.zeros((nb_per*nb_pts_per*nb_newmark*2,1))
    liste_Yt_new = np.zeros((nb_per*nb_pts_per*nb_newmark*2,1))

    for i in range(nb_newmark*2):
        tt_new, Yt_new, dYt_new, dtheta_new = Newmark(Y0[0], dY0[0], t_init[0], dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
        t_init = tt_new[-1]
        Y0 = Yt_new[-1]
        dY0 = dYt_new[-1]
        for j in prange(len(Yt)):
            liste_dtheta_new[j+i*(len(Yt)),0] = dtheta_new[j]
            liste_Yt_new[j+i*(len(Yt)),0] = Yt_new[j,0]

    return liste_dtheta,liste_Yt,liste_dtheta_new,liste_Yt_new

def data_for_saut(delta_m):
    OME, AMPL = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_0e+00_up.txt')
    OME2, AMPL2 = nk.recuperer_courbe_data('Courbes de réponse/courbe_reponse_masse_0e+00_down.txt')
    OME_data, AMPL_data = nk.recuperer_courbe_data('Courbes de réponse/courbe_reference.txt')
    OME3, AMPL3 = nk.recuperer_courbe_data(f'Courbes de réponse/courbe_reponse_masse_{delta_m}_up.txt')
    OME4, AMPL4 = nk.recuperer_courbe_data(f'Courbes de réponse/courbe_reponse_masse_{delta_m}_down.txt')
    return OME, AMPL, OME2, AMPL2, OME_data, AMPL_data, OME3, AMPL3, OME4, AMPL4

T, Vdc, Vac, omega0, M, C, K, OMEGA_bal, d = init_params()


OMEGA_min = 0.99
OMEGA_max = 0.992877
OMEGA_c = (OMEGA_max+OMEGA_min)/2

periode = 2*np.pi / OMEGA_c  # Period of excitation and response
nb_pts_per = 50     # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 500    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = np.array([0.0])               # Initial time
NT = nb_per * nb_pts_per
nb_newmark = 750
deltam = 1e-17

liste_dtheta, liste_Yt, liste_dtheta_new, liste_Yt_new = saut_ampl(T, Vdc, Vac, omega0, M, C, K, OMEGA_bal, OMEGA_min, OMEGA_max, nb_pts_per, dt, nb_per, t_init, NT, nb_newmark,deltam)
OME, AMPL, OME2, AMPL2, OME_data, AMPL_data, OME3, AMPL3, OME4, AMPL4 = data_for_saut(deltam)

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(liste_dtheta_new, np.abs(liste_Yt_new), label='Déplacement après ajout de masse (Yt_new)')
ax.plot(liste_dtheta, np.abs(liste_Yt), label='Déplacement sans ajout de masse (Yt)')
nk.plot_response_curve(OME, AMPL, OME2, AMPL2, OME_data, AMPL_data, ax, 0, True)
plot_response_curve(OME3, AMPL3, OME4, AMPL4, 0, 0, ax,deltam)
ax.set_xlabel('OMEGA adimensionné')
ax.set_ylabel('AMPLITUDE adimensionnée')
ax.grid()
plt.show()
