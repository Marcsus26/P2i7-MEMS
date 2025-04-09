# import des librairies additionnelles
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint
import time

#definition des fonctions

def calc_P(T,Vdc,Vac,OMEGA,t):
    P= T*Vdc**2 + 2*T*Vdc*Vac*np.cos(OMEGA*t)
    return P

def calc_Fnl(T,Vdc,y):
    Fnl= -(3*T*Vdc**2)*y**2 - (4*T*Vdc**2)*y**3
    return Fnl

def calc_dFnl(T,Vdc,y,dY=0):
    dFY= -6*T*(Vdc**2)*(y+2*y**2)
    dFdY= 0
    return dFY,dFdY

#Fonction newmark
def Newmark(Y0,dY0,t_init,dt,NT,omega0):
    precNR=1.e-9
    # C.I. déplacement et vitesse
    t = t_init
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
    
    P = calc_P(T,Vdc,Vac,OMEGA,t)     # calcul de l'effort extérieur
    Fnl=calc_Fnl(T,Vdc,Y)   # calcul de l'effort non-linéaire
    ddY=(P-C*dY-K*Y-Fnl)/M   # accélération initiale
    tt[0]=t
    Yt[0]=Y
    dYt[0]=dY
    
    # Boucle sur les pas de temps
    for n in range(1,NT):
        t=t+dt
        # prediction
        iter=0
        Y = Y+ dt*dY+(dt**2/2)*ddY
        dY = dY + dt*ddY
        ddY = ddY
        # Calcul du residu
        P=calc_P(T,Vdc,Vac,OMEGA,t)
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
        Yt[n]=Y     # stockage Y(n)
        dYt[n]=dY  # stockage dY(n)
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


M = 1 + (deltam/m)
C = xi
K = 1 - 2*T*Vdc**2
OMEGA = 1

periode=2*np.pi/OMEGA  # periode de l'excitation et de la reponse
nb_pts_per = 50          # nb de points par periode pour l integration temporelle
dt=periode/nb_pts_per  # taille du pas de temps
nb_per = 5000              # nb de periodes pour le calcul temporel
t_tot=nb_per*periode   # temps final
t_init=0               # temps initial
NT=nb_per*nb_pts_per+1 # Nb total de pas de temps

# Résolution equation mvt


'''
Y0=1 ; dY0=0             # conditions initiales
start = time.process_time()
tt,Yt,dYt=Newmark(Y0,dY0,t_init,dt,NT,omega0)   # Integration par Newmark
print('ellapsed time: ',time.process_time() - start)

fig = plt.figure(figsize=(10, 4))  # Crée une figure de dimensions données.
ax = fig.add_subplot()
ax.plot(tt, Yt, color='red',marker = '.', label='Newmark')  # Trace le deplacement en fonction du temps.
plt.xlabel("Temps $t$")        # titre axe horizontal
plt.ylabel("Déplacement $y(t)$")  # titre axe vertical
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.show()'''

# Courbe de réponse

OMEGA_debut=0.9;OMEGA_fin=1.10;dOMEGA=0.005
nb_pts_per=50         # nb de points par periode pour l integration temporelle
nb_per=500              # nb de periodes pour le calcul temporel
t_init=0               # temps initial
NT=nb_per*nb_pts_per   # Nb total de pas de temps
npas=int((OMEGA_fin-OMEGA_debut)/dOMEGA+1)
# conditions initiales
Y0=0.25;dY0=0
k=0; OMEGA=OMEGA_debut

OME=zeros((npas,1))
AMPL=zeros((npas,1))

# boucle sur Omega 
while (OMEGA<OMEGA_fin):
    OME[k]=OMEGA
    periode=2*np.pi/OMEGA   # periode de l'excitation et de la reponse
    dt=periode/nb_pts_per   # mise a jour du pas de temps
    t_tot=nb_per*periode    # mise a jour du temps final
    tt,Yt,dYt=Newmark(Y0,dY0,t_init,dt,NT,omega0)   # Integration par Newmark
    AMPL[k]=max(Yt[-3*nb_pts_per:])   # recherche de l’amplitude max
    print(f'ome= {OME[k,0]:0.5f}  y= {AMPL[k,0]:0.5g}', end="\r", flush=True)   # affichage resultat
    Y0=Yt[-1,0];dY0=dYt[-1,0]   # nouvelles CI
    OMEGA=OMEGA+dOMEGA
    k=k+1

fig = plt.figure(figsize=(10, 4))  # Crée une figure de dimensions données.
ax = fig.add_subplot()
ax.plot(OME, AMPL, color='red',marker = '>', label='montée en fréquence')  # Trace la courbe de réponse.
plt.xlabel(r"$\Omega$ pulsation de l'excitation")
plt.ylabel("Amplitude de la réponse = $max(y(t))$")
plt.title("Courbe de réponse")
plt.legend()
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)

k=0; OMEGA=OMEGA_fin

npas=int((OMEGA_fin-OMEGA_debut)/dOMEGA+1)
#OMEGA=np.linspace(OMEGA_debut, OMEGA_fin, num=npas)
OME2=zeros((npas,1))
AMPL2=zeros((npas,1))
# boucle sur Omega 
start = time.process_time()
while (OMEGA>OMEGA_debut):
    OME2[k]=OMEGA
    periode=2*np.pi/OMEGA      # periode de l'excitation et de la reponse
    dt=periode/nb_pts_per  # taille du pas de temps
    t_tot=nb_per*periode   # temps final
    tt,Yt,dYt=Newmark(Y0,dY0,t_init,dt,NT,omega0)   # Integration par Newmark
    AMPL2[k]=max(Yt[-3*nb_pts_per:])  # recherche de l’amplitude max
    print(f'ome= {OME2[k,0]:0.5f}  x= {AMPL2[k,0]:0.5g}', end="\r", flush=True)   # affichage resultat
    Y0=Yt[-1,0];dY0=dYt[-1,0]   # nouvelles CI
    OMEGA=OMEGA-dOMEGA
    k=k+1

ax.plot(OME2, AMPL2, color='blue',marker = '<', label='descente en fréquence')  # Trace la courbe de réponse.

#Récupération des données de la courbe de reponse prof
data = np.loadtxt('courbe_reponse_modified.txt', delimiter=',')
OMEGA_data = data[:, 0]
AMPL_data = data[:, 1]

# Ajout des données du fichier au graphe existant
ax.plot(OMEGA_data, AMPL_data, color='green', marker='o', label='données fichier')  # Trace les données du fichier.

plt.legend()

plt.show()
