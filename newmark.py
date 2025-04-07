# import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint
import time

#definition des fonctions

def calc_P(T,Vdc,Vac,OMEGA,tau):
    P= T*Vdc**2 + 2*T*Vdc*Vac*np.cos(OMEGA*tau)
    return P

def calc_Fnl(T,Vdc,y):
    Fnl= -(3*T*Vdc**2)*y**2 - (4*T*Vdc**2)*y**3
    return Fnl

def calc_dFnl(T,Vdc,y,dY=0):
    dFY= -6*T*Vdc**2*(y+2*y**2)
    dFdY= 0
    return dFY,dFdY

#Fonction newmark
def Newmark(X0,dX0,t_init,dt,NT,omega0):
    precNR=1.e-9
    # C.I. déplacement et vitesse
    t = t_init
    tau = omega0*t
    Y=Y0;dY=dY0
    # initialisation des variables
    Fnl=0
    dFY=0
    dFdY=0
    P=0
    # tableaux de stockage des résultats à chaque pas de temps
    tt=zeros((NT,1))      # temps
    Yt=zeros((NT,1))      # déplacement
    dYt=zeros((NT,1))     # vitesse
    
    P = calc_P(T,Vdc,Vac,OMEGA,tau)     # calcul de l'effort extérieur
    Fnl=calc_Fnl(T,Vdc,Y)   # calcul de l'effort non-linéaire
    ddY=(P-C*dY-K*Y-Fnl)/M   # accélération initiale
    tt[0]=t
    Yt[0]=Y
    dYt[0]=dY
    
    # Boucle sur les pas de temps
    for n in range(1,NT):
        t=t+dt
        tau = omega0*t
        # prediction
        iter=0
        Y = Y+ dt*dY+(dt**2/2)*ddY
        dY = dY + dt*ddY
        ddY = ddY
        # Calcul du residu
        P=calc_P(T,Vdc,Vac,OMEGA,tau)
        Fnl=calc_Fnl(T,Vdc,Y)
        res=P-M*ddY-C*dY-K*Y-Fnl
        normres=abs(res/P)
        # Iterations de Newton Raphson
        while (normres>precNR):
            iter=iter+1
            # Calcul de la Jacobienne
            dFY,dFdY=calc_dFnl(T,Vdc,Y)
            J=(4/dt**2)*M+(2/dt)*(C+dFdY)+K+dFY
            # Calcul de la correction
            deltaY=res/J
            Y = Y + deltaY
            # mise à jour vitesse et accélération
            dY = dY + (2/dt)*deltaY
            ddY = ddY+(4/dt**2)*deltaY
            # Calcul du residu
            Fnl=calc_Fnl(T,Vdc,Y)
            res=P-M*ddY-C*dY-K*Y-Fnl
            #normres=abs(res/P)
            normres=abs(deltaY/Y)
        tt[n]=t     # stockage t(n)
        Yt[n]=Y     # stockage X(n)
        dYt[n]=dY  # stockage dX(n)
    return tt,Yt,dYt

#Déclaration des parametres du système
rho = 2500
l = 250e-6
b = 40e-6
h = 1e-6
d = 0.03e-6
Vdc = 5
Vac = Vdc/10
epsilon0 = 8.5e-12

f0 = 100e6
omega0 = 2*np.pi*f0

Q = 1000
xi = 1/Q

A = l*b
deltam = 0
m = rho*A*h
T = (epsilon0*A)/(2*m*omega0**2*d**3)


M = 1 + deltam/m
C = xi
K = 1 - 2*T*Vdc**2
OMEGA = 1

periode=2*np.pi/OMEGA  # periode de l'excitation et de la reponse
nb_pts_per = 50          # nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per  # taille du pas de temps
nb_per = 2000              # nb de periodes pour le calcul temporel
t_tot=nb_per*periode   # temps final
t_init=0               # temps initial
NT=nb_per*nb_pts_per+1 # Nb total de pas de temps

# Résolution equation mvt

OMEGA=1

Y0=1 ; dY0=0             # conditions initiales
start = time.process_time()
tt,Yt,dYt=Newmark(Y0,dY0,t_init,dt,NT,omega0)   # Integration par Newmark
print('ellapsed time: ',time.process_time() - start)

fig = plt.figure(figsize=(10, 4))  # Crée une figure de dimensions données.
ax = fig.add_subplot()
ax.plot(tt, Yt, color='red',marker = 'o', label='Newmark')  # Trace le deplacement en fonction du temps.
plt.xlabel("Temps $t$")        # titre axe horizontal
plt.ylabel("Déplacement $y(t)$")  # titre axe vertical
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.show()

#