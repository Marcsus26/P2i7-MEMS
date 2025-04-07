# import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

#definition des fonctions

def calc_P(p0,OMEGA,t):
    P= p0*np.cos(OMEGA*t)
    return P

def calc_Fnl(knl,X,dX):
    Fnl= knl*X**3
    return Fnl

def calc_dFnl(knl,X,dX):
    dFX= 3*knl*X**2
    dFdX= 0
    return dFX,dFdX

#Fonction newmark
def Newmark(X0,dX0,t_init,dt,NT):
    precNR=1.e-9
    # C.I. déplacement et vitesse
    t=t_init
    X=X0;dX=dX0
    # initialisation des variables
    Fnl=0
    dFX=0
    dFdX=0
    P=0
    # tableaux de stockage des résultats à chaque pas de temps
    tt=zeros((NT,1))      # temps
    Xt=zeros((NT,1))      # déplacement
    dXt=zeros((NT,1))     # vitesse
    
    P=calc_P(p0,OMEGA,t)     # calcul de l'effort extérieur
    Fnl=calc_Fnl(knl,X,dX)   # calcul de l'effort non-linéaire
    ddX=(P-C*dX-K*X-Fnl)/M   # accélération initiale
    tt[0]=t
    Xt[0]=X
    dXt[0]=dX
    
    # Boucle sur les pas de temps
    for n in range(1,NT):
        t=t+dt;
        # prediction
        iter=0;
        X= X+ dt*dX+(dt**2/2)*ddX
        dX= dX + dt*ddX
        ddX= ddX
        # Calcul du residu
        P=calc_P(p0,OMEGA,t)
        Fnl=calc_Fnl(knl,X,dX)
        res=P-M*ddX-C*dX-K*X-Fnl
        normres=abs(res/P)
        # Iterations de Newton Raphson
        while (normres>precNR):
            iter=iter+1
            # Calcul de la Jacobienne
            dFX,dFdX=calc_dFnl(knl,X,dX)
            J=(4/dt**2)*M+(2/dt)*(C+dFdX)+K+dFX
            # Calcul de la correction
            deltaX=res/J
            X= X + deltaX
            # mise à jour vitesse et accélération
            dX= dX + (2/dt)*deltaX
            ddX= ddX+(4/dt**2)*deltaX
            # Calcul du residu
            Fnl=calc_Fnl(knl,X,dX)
            res=P-M*ddX-C*dX-K*X-Fnl
            #normres=abs(res/P)
            normres=abs(deltaX/X)
        tt[n]=t     # stockage t(n)
        Xt[n]=X     # stockage X(n)
        dXt[n]=dX   # stockage dX(n)
    return tt,Xt,dXt

#Déclaration des parametres du système

M=1
C=0.1
K=1
knl=0.25
p0=.5
OMEGA=1.

periode=2*np.pi/OMEGA  # periode de l'excitation et de la reponse
nb_pts_per=30          # nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per  # taille du pas de temps
nb_per=30              # nb de periodes pour le calcul temporel
t_tot=nb_per*periode   # temps final
t_init=0               # temps initial
NT=nb_per*nb_pts_per+1 # Nb total de pas de temps

# Résolution equation mvt

OMEGA=1.
X0=1;dX0=0             # conditions initiales
start = time.process_time()
tt,Xt,dXt=Newmark(X0,dX0,t_init,dt,NT)   # Integration par Newmark
print('ellapsed time: ',time.process_time() - start)

fig = plt.figure(figsize=(10, 4))  # Crée une figure de dimensions données.
ax = fig.add_subplot()
ax.plot(tt, Xt, color='red',marker = 'o', label='Newmark')  # Trace le deplacement en fonction du temps.
plt.xlabel("Temps $t$")        # titre axe horizontal
plt.ylabel("Déplacement $x(t)$")  # titre axe vertical
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.show()

#