import numpy as np

def newmark_beta(m, c, k, f, dt, u0, v0, beta=0.25, gamma=0.5):
    """
    Implémentation de l'algorithme de Newmark pour résoudre des systèmes dynamiques.

    Arguments :
    m -- masse
    c -- coefficient d'amortissement
    k -- raideur
    f -- force appliquée (tableau)
    dt -- pas de temps
    u0 -- condition initiale de déplacement
    v0 -- condition initiale de vitesse
    beta -- paramètre de Newmark (par défaut 0.25 pour méthode implicite)
    gamma -- paramètre de Newmark (par défaut 0.5 pour méthode implicite)

    Retourne :
    u -- tableau des déplacements
    v -- tableau des vitesses
    a -- tableau des accélérations
    """
    n = len(f)
    u = np.zeros(n)
    v = np.zeros(n)
    a = np.zeros(n)

    # Conditions initiales
    u[0] = u0
    v[0] = v0
    a[0] = (f[0] - c * v0 - k * u0) / m

    # Coefficients constants
    a0 = 1 / (beta * dt**2)
    a1 = gamma / (beta * dt)
    a2 = 1 / (beta * dt)
    a3 = 1 / (2 * beta) - 1
    a4 = gamma / beta - 1
    a5 = dt * (gamma / (2 * beta) - 1)

    # Matrice effective
    keff = k + a0 * m + a1 * c

    for i in range(1, n):
        # Force effective
        feff = f[i] + m * (a0 * u[i-1] + a2 * v[i-1] + a3 * a[i-1]) + c * (a1 * u[i-1] + a4 * v[i-1] + a5 * a[i-1])

        # Déplacement
        u[i] = feff / keff

        # Accélération
        a[i] = a0 * (u[i] - u[i-1]) - a2 * v[i-1] - a3 * a[i-1]

        # Vitesse
        v[i] = v[i-1] + dt * ((1 - gamma) * a[i-1] + gamma * a[i])

    return u, v, a

print('merci mon coeur <3' \
'oueeee')