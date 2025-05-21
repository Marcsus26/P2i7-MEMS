import numpy as np
import matplotlib.pyplot as plt
import newmark_saut as nks
import newmark as nk
import test as tst

def saut_amp(T, Vdc, Vac, omega0, M, C, K, l_masses):
    OMEGA_min = 0.99
    OMEGA_max = 0.993
    T, Vdc, Vac, omega0, M_new, C, K, OMEGA_bal = nks.init_params()
    T_bal = 2 * np.pi / OMEGA_bal  # Period of excitation and response
    nb_pts_per = 60             # Number of points per period for time integration
    dt = T_bal / nb_pts_per   # Time step size
    nb_per = 50000          # Number of periods for time integration
    t_tot = nb_per * T_bal   # Final time
    t_init = 0                  # Initial time
    NT = nb_per * nb_pts_per
    M0 = M

    # Updated initial conditions
    Y0 = 0.3  # Initial displacement
    dY0 = 0.3  # Initial velocity
    fig,ax = plt.subplots(figsize = (16,9))
    OME_montee = np.linspace(OMEGA_min, OMEGA_max, NT)
    OME_descente = np.linspace(OMEGA_max, OMEGA_min, NT)
    OME = np.concatenate((OME_montee,OME_descente))
    AMPL = np.zeros(2*NT)

    T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = nks.init_params()
    data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
    OMEGA_data, AMPL_data = data[:, 0], data[:, 1]
    tst.plot_response_curve(0,0,0,0,OMEGA_data, AMPL_data, ax, 0, True)
    tt, Yt, dYt = nks.Newmark(Y0,dY0,t_init,dt,NT,omega0,T,Vdc,Vac,OMEGA_min,OMEGA_max,M,C,K,OMEGA_bal)
    AMPL[0:NT] = max(max(Yt[-3 * nb_pts_per:]))
    t_init_new = tt[-1]
    Y0_new = Yt[-1]
    dY0_new = dYt[-1]
    tt_new, Yt_new, dYt_new = nks.Newmark(Y0_new,dY0_new,t_init_new,dt,NT,omega0,T,Vdc,Vac,OMEGA_max,OMEGA_min,M,C,K,OMEGA_bal)
    AMPL[NT:2*NT] = max(max(Yt_new[-3 * nb_pts_per:]))
    ax.plot(OME,AMPL)

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Displacement')
    ax.set_title('Displacement Over Time with and without Added Mass')
    ax.legend()
    ax.grid()
    plt.show()

# Initialize parameters
T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = nks.init_params()
liste_masse = [1e-13,0]
# liste_masse = [1e-25, 1e-14]
saut_amp(T, Vdc, Vac, omega0, M, C, K, liste_masse)