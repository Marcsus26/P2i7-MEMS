import newmark as nk
import newmark_masse_ajoutee as nkm
import numpy as np
import newmark_saut as nks
import matplotlib.pyplot as plt



# # Initialize parameters
# T, Vdc, Vac, omega0, M, C, K = nk.init_params()
# OMEGA = 1
# periode = 2 * np.pi / OMEGA  # Period of excitation and response
# nb_pts_per = 50         # Number of points per period for time integration
# dt = periode / nb_pts_per   # Time step size
# nb_per = 5000            # Number of periods for time integration
# t_tot = nb_per * periode    # Final time
# t_init = 0                  # Initial time
# NT = nb_per * nb_pts_per

# tt, Yt, dYt = nk.Newmark(0.3, 0.3, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)


# # Plot the displacement over time
# plt.figure(figsize=(50, 6))
# plt.plot(tt / omega0, Yt, label='Déplacement (Yt)')
# plt.xlabel('Temps (s)')
# plt.ylabel('Déplacement')
# plt.title('Déplacement au cours du temps pour deltam = 0')
# plt.legend()
# plt.xlim(0, 2e-5)
# plt.grid()
# plt.show()


def plot_response_curve(OME, AMPL, OME2, AMPL2, OMEGA_data, AMPL_data, ax, deltam=0, tracer_data=False):


    ax.plot(OME, AMPL, marker='>', label=f'montée en fréquence pour {deltam}')


    ax.plot(OME2, AMPL2, marker='<', label=f'descente en fréquence pour {deltam}')

    if tracer_data:

        ax.plot(OMEGA_data, AMPL_data, color='green', marker='o', label='données fichier')


    ax.plot(OME, AMPL, marker='>', label=f'montée en fréquence pour {deltam}')


    ax.plot(OME2, AMPL2, marker='<', label=f'descente en fréquence pour {deltam}')

    plt.xlabel(r"$\Omega$ pulsation de l'excitation")

    plt.ylabel("Amplitude de la réponse = $max(y(t))$")

    plt.title("Courbe de réponse")

    plt.legend()

    plt.xlim(0.988,0.996)

    plt.grid(color='gray', linestyle='--', linewidth=0.5)

T, Vdc, Vac, omega0, M, C, K = nk.init_params()
OMEGA = 1
periode = 2*np.pi / OMEGA  # Period of excitation and response
nb_pts_per = 100      # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 2500    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                  # Initial time
NT = nb_per * nb_pts_per

tt, Yt, dYt  = nk.Newmark(0.3,0.3,t_init,dt, NT, omega0, T, Vdc, Vac,OMEGA,M,C,K)

plt.plot(tt/omega0, Yt)
plt.show()