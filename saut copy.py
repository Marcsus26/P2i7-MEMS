import newmark as nk
import numpy as np
import matplotlib.pyplot as plt

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

    for i in range (len(l_masses)):
        for omega in range(-50, 50, 1):
            tt, Yt, dYt = nk.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA+(omega*1e-3), M0, C, K)
            if np.max(Yt[-3 * nb_pts_per:]) < 0.2:
                break
        T, Vdc, Vac, omega0, M_new, C, K = nk.init_params(l_masses[i])
        ax.plot(tt/(omega0), Yt, label='Displacement (Yt) - deltam=0')
        t_init = tt[-1]
        Y0 = Yt[-1]
        dY0 = dYt[-1]
        tt, Yt, dYt = nk.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M_new, C, K)
        print(np.max(Yt[-3 * nb_pts_per:]))
        ax.plot(tt/(omega0), Yt, label=f'Displacement (Yt) - deltam={l_masses[i]}')
        Y0 = Yt[-1]
        dY0 = dYt[-1]
        t_init = tt[-1]

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Displacement')
    ax.set_title('Displacement Over Time with and without Added Mass')
    ax.legend()
    ax.grid()
    plt.show()

# Initialize parameters
T, Vdc, Vac, omega0, M, C, K = nk.init_params()
liste_masse = [1e-13,1e-14, 1e-15, 1e-16, 1e-17, 1e-18]
saut_amp(T, Vdc, Vac, omega0, M, C, K, liste_masse)