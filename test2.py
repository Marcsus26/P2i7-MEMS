from newmark_saut import *
import newmark as nk

T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params()

# Initialize parameters

OMEGA_min = 0.98
OMEGA_max = 0.9928
OMEGA_c = (OMEGA_max+OMEGA_min)/2

periode = 2*np.pi / OMEGA_c  # Period of excitation and response
nb_pts_per = 500     # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 10000    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                # Initial time
NT = nb_per * nb_pts_per

tt, Yt, dYt = Newmark(0.5, 0.5, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)
t_init = tt[-1]
Y0 = Yt[-1]
dY0 = dYt[-1]
T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params(1e-15)
tt_new, Yt_new, dYt_new = Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA_min, OMEGA_max, M, C, K, OMEGA_bal)

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(tt / omega0, (Yt), label='Déplacement (Yt)')
ax.plot(tt_new / omega0, (Yt_new), label='Déplacement après ajout de masse (Yt_new)')
ax.set_xlabel('Temps (s)')
ax.set_ylabel('Déplacement')
ax.grid()
plt.show()
