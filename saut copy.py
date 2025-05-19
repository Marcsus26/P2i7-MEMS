import newmark as nk
import numpy as np
import matplotlib.pyplot as plt
from numba import njit, prange

@njit
def saut_amp(T, Vdc, Vac, omega0, M, C, K, l_masses):
    OMEGA = 0.9923
    periode = 2 * np.pi / OMEGA  # Period of excitation and response
    nb_pts_per = 60             # Number of points per period for time integration
    dt = periode / nb_pts_per   # Time step size
    nb_per = 5000               # Number of periods for time integration
    t_tot = nb_per * periode    # Final time
    t_init = 0                  # Initial time
    NT = nb_per * nb_pts_per
    M0 = M

    # Updated initial conditions
    Y0 = 0.3  # Initial displacement
    dY0 = 0.3  # Initial velocity
    fig,ax = plt.subplots(figsize = (16,9))

    for i in prange(len(l_masses)):
        saut = [-float('inf')]
        T, Vdc, Vac, omega0, M_new, C, K = nk.init_params(l_masses[i])
        for omega in prange(-50, 50, 1):
            tt, Yt, dYt = nk.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA+(omega*1e-3), M0, C, K)
            ampl_0 = np.max(Yt[-3 * nb_pts_per:])
            Y0 = Yt[-1]
            dY0 = dYt[-1]
            tt_new, Yt_new, dYt_new = nk.Newmark(Y0, dY0, tt[-1], dt, NT, omega0, T, Vdc, Vac, OMEGA+(omega*1e-3), M_new, C, K)
            ampl = np.max(Yt_new[-3 * nb_pts_per:])
            saut.append(ampl-ampl_0)
            Y0 = Yt[0]
            dY0 = dYt[0]
            print('hey')
        omega_max = np.argmax(saut)
        tt, Yt, dYt = nk.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA+((omega_max-50)*1e-3), M0, C, K)
        t_init = tt[-1]
        Y0 = Yt[-1]
        dY0 = dYt[-1]
        tt_new, Yt_new, dYt_new = nk.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA+((omega_max-50)*1e-3), M_new, C, K) 
        Y0 = Yt_new[-1]
        dY0 = dYt_new[-1]
        t_init = tt_new[-1]
        ax.plot(tt/(omega0), Yt, label='Displacement (Yt) - deltam=0')
        ax.plot(tt_new/(omega0), Yt_new, label=f'Displacement (Yt) - deltam={l_masses[i]}')

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Displacement')
    ax.set_title('Displacement Over Time with and without Added Mass')
    ax.legend()
    ax.grid()
    plt.show()

# Initialize parameters
T, Vdc, Vac, omega0, M, C, K = nk.init_params()
liste_masse = [1e-13,1e-14]
saut_amp(T, Vdc, Vac, omega0, M, C, K, liste_masse)