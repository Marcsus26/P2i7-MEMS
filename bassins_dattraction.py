import newmark as nk
import numpy as np
import matplotlib.pyplot as plt

T, Vdc, Vac, omega0, M, C, K, d = nk.init_params(1e-13)

OMEGA = 0.99
periode = 2 * np.pi / OMEGA  # periode de l'excitation et de la reponse
nb_pts_per = 60             # nb de points par periode pour l integration temporelle
dt = periode / nb_pts_per   # taille du pas de temps
nb_per = 2500               # nb de periodes pour le calcul temporel
t_tot = nb_per * periode    # temps final
t_init = 0                  # temps initial
NT = nb_per * nb_pts_per

# Conditions initiales
# pas de la grille
# pasY = 0.5; pasdY = 1     # grille grossiere (837 cases)
# pasY = 0.2; pasdY = .4    # grille moyenne (5016 cases)
pasY = 0.025; pasdY = 0.05    # grille fine (19781 cases)
# bornes de la grille
Y0 = np.arange(-1, 1 + pasY, pasY)
dY0 = np.arange(-2, 2 + pasdY, pasdY)

AMPL = np.zeros((len(dY0), len(Y0)))

fig1 = plt.figure(figsize=(6, 5))  # Crée une figure de dimensions données.
ax1 = fig1.add_subplot()
ax1.set_xlabel('Y0')
ax1.set_ylabel('dY0')
im = ax1.imshow(AMPL, interpolation='none', aspect='auto', extent=[-2, 2, -2, 2], vmin=0, vmax=1, cmap='jet')
fig1.colorbar(im, ax=ax1)

plt.ion()  # Turn on interactive mode for live updates
plt.show()

for i in range(len(Y0)):  # boucle sur les CI
    for j in range(len(dY0)):
        tt, Yt, dYt = nk.Newmark(Y0[i], dY0[j], t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)
        AMPL[j, i] = np.max(Yt[-2 * nb_pts_per:])
        im.set_data(AMPL)
    fig1.canvas.draw()
    fig1.canvas.flush_events()

plt.ioff()  # Turn off interactive mode
plt.show()