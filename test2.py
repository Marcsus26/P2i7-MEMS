from newmark_saut import *
import numpy as np
import newmark as nk
from numba import prange

T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params()

# Initialize parameters

OMEGA_min = 0.985
OMEGA_max = 0.993
OMEGA_c = (OMEGA_max+OMEGA_min)/2

periode = 2*np.pi / OMEGA_c  # Period of excitation and response
nb_pts_per = 100     # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 100000    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                # Initial time
NT = nb_per * nb_pts_per

tt, Yt, dYt, theta = Newmark(0, 0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)

t_init = tt[-1]
Y0 = Yt[-1]
dY0 = dYt[-1]
T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params(1e-13)
tt_new, Yt_new, dYt_new, theta_new = Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
# OME, AMPL = compute_response_curve(T, Vdc, Vac, omega0, M, C, K, 0.97, 1, 0.00001, 60, nb_per, OMEGA_bal, tolerance=0.00001)
# OME2, AMPL2 = compute_response_curve(T, Vdc, Vac, omega0, M, C, K, 1, 0.97, 0.00001, 60, nb_per, OMEGA_bal, tolerance=0.00001)

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(tt/omega0, (Yt), label='Déplacement (Yt)')
ax.plot(tt_new / omega0, (Yt_new), label='Déplacement après ajout de masse (Yt_new)')
# plot_response_curve(OME, AMPL, OME2, AMPL2, 0, 0, ax, deltam=0, tracer_data=False)
ax.set_xlabel('Temps (s)')
ax.set_ylabel('Déplacement')
ax.grid()
plt.show()
