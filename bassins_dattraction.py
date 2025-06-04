import newmark as nk
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

def bassin_attraction(deltam):
    T, Vdc, Vac, omega0, M, C, K, d = nk.init_params(deltam)
    OMEGA = 0.99
    nb_pts_per = 60
    nb_per = 2500
    periode = 2 * np.pi / OMEGA
    dt = periode / nb_pts_per
    t_init = 0
    NT = nb_per * nb_pts_per
    pasY = 0.025
    pasdY = 0.05
    Y0 = np.arange(-1, 1 + pasY, pasY)
    dY0 = np.arange(-2, 2 + pasdY, pasdY)
    AMPL = np.zeros((len(dY0), len(Y0)))

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.set_xlabel('Y0')
    ax.set_ylabel('dY0')
    im = ax.imshow(AMPL, interpolation='none', aspect='auto', extent=[-1, 1, -2, 2], vmin=0, vmax=1, cmap='jet')
    fig.colorbar(im, ax=ax)

    for i in range(len(Y0)):
        for j in range(len(dY0)):
            tt, Yt, dYt = nk.Newmark(Y0[i], dY0[j], t_init, dt, NT, omega0, T, Vdc, Vac, OMEGA, M, C, K)
            AMPL[j, i] = np.max(Yt[-2 * nb_pts_per:])
            im.set_data(AMPL)
        clear_output(wait=True)
        display(fig)
        plt.pause(0.01)

    plt.close(fig)