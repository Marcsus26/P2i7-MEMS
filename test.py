import newmark as nk
import newmark_masse_ajoutee as nkm
import numpy as np
import newmark_saut as nks
import matplotlib.pyplot as plt



T, Vdc, Vac, omega0, M, C, K = nk.init_params()
OMEGA = 1
periode = 2*np.pi / OMEGA  # Period of excitation and response
nb_pts_per = 100      # Number of points per period for time integration
dt = periode / nb_pts_per   # Time step size
nb_per = 2500    # Number of periods for time integration
t_tot = nb_per * periode    # Final time
t_init = 0                  # Initial time
NT = nb_per * nb_pts_per
print(dt)

tt, Yt, dYt  = nk.Newmark(0.3,0.3,t_init,dt, NT, omega0, T, Vdc, Vac,OMEGA,M,C,K)

plt.plot(tt/omega0, Yt)
plt.show()