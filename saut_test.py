import newmark as nk
import numpy as np
import matplotlib.pyplot as plt

# Initialize parameters
T, Vdc, Vac, omega0, M, C, K = nk.init_params()

OMEGA = 0.9923
periode = 2 * np.pi / OMEGA  # Period of excitation and response
nb_pts_per = 60             # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 5000               # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                  # Initial time
NT = nb_per * nb_pts_per

# Updated initial conditions
Y0 = 0.3  # Initial displacement
dY0 = 0.3  # Initial velocity

# Perform Newmark integration
tt, Yt, dYt = nk.Newmark(Y0, dY0, t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)

# # Plot the displacement over time
# plt.figure(figsize=(10, 6))
# plt.plot(tt, Yt, label='Displacement (Yt)')
# plt.xlabel('Time (s)')
# plt.ylabel('Displacement')
# plt.title('Displacement Over Time for deltam=0')
# plt.legend()
# plt.grid()


# Verify the amplitude of the periodic solution
amplitude_periodic_solution = np.max(Yt[-3 * nb_pts_per:])
print(f"Amplitude of the periodic solution: {amplitude_periodic_solution}")
#assert np.isclose(amplitude_periodic_solution, 0.169, atol=0.01), "The amplitude does not match the expected value of 0.169."

# Add a small mass and continue integration

# Update parameters with the added mass
T, Vdc, Vac, omega0, M_new, C, K = nk.init_params(1e-13)

# Use final conditions from the previous integration as initial conditions
Y0_new = Yt[-1]
dY0_new = dYt[-1]

# Perform Newmark integration for the extended time period
tt_new, Yt_new, dYt_new = nk.Newmark(Y0_new, dY0_new, t_tot, dt, NT, omega0, T, Vdc, Vac, OMEGA, M_new, C, K)

amplitude_periodic_solution = np.max(Yt_new[-3 * nb_pts_per:])
print(f"Amplitude of the periodic solution: {amplitude_periodic_solution}")

# Plot the displacement over time for both integrations
# plt.figure(figsize=(10, 6))
# plt.plot(tt, Yt, label='Displacement (Yt) - deltam=0')
# plt.plot(tt_new, Yt_new, label='Displacement (Yt) - deltam=1e-13', color='orange')
# plt.xlabel('Time (s)')
# plt.ylabel('Displacement')
# plt.title('Displacement Over Time with and without Added Mass')
# plt.legend()
# plt.grid()
# plt.show()

T, Vdc, Vac, omega0, M_new, C, K = nk.init_params()
Y0_new2 = Yt_new[-1]
dY0_new2 = dYt_new[-1]

tt_new2, Yt_new2, dYt_new2 = nk.Newmark(Y0_new2, dY0_new2, 2*t_tot, dt, NT, omega0, T, Vdc, Vac, OMEGA, M_new, C, K)

amplitude_periodic_solution = np.max(Yt_new2[-3 * nb_pts_per:])
print(f"Amplitude of the periodic solution: {amplitude_periodic_solution}")

# plt.figure(figsize=(10, 6))
# plt.plot(tt, Yt, label='Displacement (Yt) - deltam=0')
# plt.plot(tt_new, Yt_new, label='Displacement (Yt) - deltam=1e-13', color='orange')
# plt.plot(tt_new2, Yt_new2, label='Displacement (Yt) - deltam=0')
# plt.xlabel('Time (s)')
# plt.ylabel('Displacement')
# plt.title('Displacement Over Time with and without Added Mass')
# plt.legend()
# plt.grid()
# plt.show()
T, Vdc, Vac, omega0, M_new, C, K = nk.init_params()
Y0_new3 = Yt_new2[-1]
dY0_new3 = dYt_new2[-1]

tt_new3, Yt_new3, dYt_new3 = nk.Newmark(Y0_new3, dY0_new3, 3*t_tot, dt, NT, omega0, T, Vdc, Vac, OMEGA, M_new, C, K)

amplitude_periodic_solution = np.max(Yt_new3[-3 * nb_pts_per:])
print(f"Amplitude of the periodic solution: {amplitude_periodic_solution}")

plt.figure(figsize=(10, 6))
plt.plot(tt/(omega0), Yt, label='Displacement (Yt) - deltam=0')
plt.plot(tt_new/omega0, Yt_new, label='Displacement (Yt) - deltam=1e-13', color='orange')
plt.plot(tt_new2/omega0, Yt_new2, label='Displacement (Yt) - deltam=0')
plt.plot(tt_new3/omega0, Yt_new3, label='Displacement (Yt) - deltam=1e-14')
plt.xlabel('Time (s)')
plt.ylabel('Displacement')
plt.title('Displacement Over Time with and without Added Mass')
plt.legend()
plt.grid()
plt.show()