import numpy as np
import matplotlib.pyplot as plt
import newmark_saut as nks

def saut_amp(T, Vdc, Vac, omega0, M, C, K, l_masses):
    OMEGA = 0.9928
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

    for i in range(len(l_masses)):
        T, Vdc, Vac, omega0, M_new, C, K, OMEGA_bal = nks.init_params(l_masses[i])
        tt, Yt, dYt = nks.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M0, C, K, OMEGA_bal)
        ax.plot(tt/(omega0), Yt, label='Displacement (Yt) - deltam=0')
        t_init = tt[-1]
        Y0 = Yt[-1]
        dY0 = dYt[-1]
        tt, Yt, dYt = nks.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M_new, C, K, OMEGA_bal)
        print(np.max(Yt[-3 * nb_pts_per:]))
        ax.plot(tt/(omega0), Yt, label=f'Displacement (Yt) - deltam={l_masses[i]}')
        Y0 = Yt[-1]
        dY0 = dYt[-1]
        t_init = tt[-1]
    
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Displacement')
    ax.set_title('Displacement Over Time with and without Added Mass')
    ax.grid()
    plt.show()

# Initialize parameters
T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = nks.init_params()
liste_masse = [1e-13,0]
# liste_masse = [1e-25, 1e-14]
saut_amp(T, Vdc, Vac, omega0, M, C, K, liste_masse)