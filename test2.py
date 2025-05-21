from newmark_saut import *

T, Vdc, Vac, omega0, M, C, K, OMEGA_bal = init_params()

# Initialize parameters

periode = 1 / OMEGA_bal  # Period of excitation and response
nb_pts_per = 100      # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 50      # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                  # Initial time
NT = nb_per * nb_pts_per

tt, Yt, dYt = Newmark(0.3, 0.3, t_init, dt, NT, omega0, T, Vdc, Vac, 0.99, 0.993, M, C, K, OMEGA_bal)

# Plot the displacement over time
plt.figure(figsize=(50, 6))
plt.plot(tt / omega0, Yt, label='Déplacement (Yt)')
plt.xlabel('Temps (s)')
plt.ylabel('Déplacement')
plt.title('Déplacement au cours du temps pour deltam = 0')
plt.legend()
plt.grid()
plt.show()