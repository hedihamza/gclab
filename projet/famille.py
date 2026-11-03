"""La famille x^3 - 2x + c, et ce que devient Newton parti de 0 quand c varie."""

import numpy as np

from bassins import iterer


def orbite_zero(C, tours):
    """Position après tours pas de Newton depuis 0, pour chaque valeur de c du tableau C."""
    return iterer([C, -2, 0, 1], np.zeros_like(C), tours)


def converge_depuis_zero(C, tours=4000, tol=1e-9):
    """Vrai pour les c dont l'orbite de 0 finit sur une racine."""
    Z = orbite_zero(C, tours)
    with np.errstate(all="ignore"):
        return np.abs(Z ** 3 - 2 * Z + C) < tol


def tours_depuis_zero(C, tours=200, tol=1e-10):
    """Nombre de tours pour que l'orbite de 0 arrive sur une racine (tours + 1 si jamais)."""
    Z = np.zeros_like(C)
    nb = np.full(C.shape, tours + 1)
    with np.errstate(all="ignore"):
        for k in range(tours):
            Z = Z - (Z ** 3 - 2 * Z + C) / (3 * Z ** 2 - 2)
            arrive = (nb > tours) & (np.abs(Z ** 3 - 2 * Z + C) < tol)
            nb[arrive] = k + 1
    return nb


def periode(c, tours=4000, tol=1e-9, qmax=64):
    """Période du cycle sur lequel finit l'orbite de 0, ou None si on n'en trouve pas."""
    z = 0.0
    traj = []
    for _ in range(tours):
        z = z - (z * z * z - 2 * z + c) / (3 * z * z - 2)
        traj.append(z)
    for q in range(1, qmax + 1):
        if abs(traj[-1 - q] - traj[-1]) < tol:
            return q
    return None


def bifurcation(C, tours=3000, garder=64):
    """Les garder positions qui suivent tours pas de Newton depuis 0, pour chaque c."""
    Z = orbite_zero(C, tours)
    res = []
    with np.errstate(all="ignore"):
        for _ in range(garder):
            Z = Z - (Z ** 3 - 2 * Z + C) / (3 * Z ** 2 - 2)
            res.append(Z.real.copy())
    return np.array(res)


def retour(c, k):
    """Position de l'orbite de 0 après 2^k pas. Elle vaut 0 si 0 est sur un cycle de période 2^k.

    Le pas est écrit en clair plutôt qu'avec evaluer : c'est trois fois plus rapide, et
    pour k = 16 on fait 65 536 pas par appel.
    """
    z = 0.0
    for _ in range(2 ** k):
        z = z - (z * z * z - 2 * z + c) / (3 * z * z - 2)
    return z


def retour_bis(c, k):
    """Comme retour, avec le pas de Newton mis sous la forme (2z^3 - c) / (3z^2 - 2)."""
    z = 0.0
    for _ in range(2 ** k):
        z = (2 * z * z * z - c) / (3 * z * z - 2)
    return z


def dichotomie(f, a, b):
    """Zéro de f entre a et b (signes opposés), en coupant tant que les flottants le permettent."""
    fa = f(a)
    while True:
        m = (a + b) / 2
        if m == a or m == b:
            return m
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm


def superstables(kmax, f=retour):
    """Dictionnaire k -> s_k : le c, juste sous le précédent, où 0 a la période 2^k.

    s_1 = 2. On cherche s_2 dans [1.95 ; 1.965], repéré sur le diagramme de bifurcation.
    Ensuite on prédit s_k avec le dernier rapport des écarts, on cherche un changement de
    signe de f(c, k) autour de la prédiction, et on termine par dichotomie.
    """
    s = {1: 2.0}
    delta = 4.0
    for k in range(2, kmax + 1):
        if k == 2:
            cs = np.linspace(1.95, 1.965, 301)
        else:
            d = s[k - 2] - s[k - 1]
            centre = s[k - 1] - d / delta
            cs = np.linspace(centre - 0.3 * d / delta, centre + 0.3 * d / delta, 61)
        cs = [float(c) for c in cs]  # les float de Python vont plus vite que ceux de numpy ici
        vals = [f(c, k) for c in cs]
        trouve = None
        for i in range(len(cs) - 1, 0, -1):
            # |f| < 0.5 écarte les sauts de signe dus à un passage près d'un zéro de p'
            if vals[i] * vals[i - 1] < 0 and abs(vals[i]) < 0.5 and abs(vals[i - 1]) < 0.5:
                trouve = dichotomie(lambda c: f(c, k), cs[i - 1], cs[i])
                break
        if trouve is None:
            break
        s[k] = trouve
        if k >= 3:
            delta = (s[k - 2] - s[k - 1]) / (s[k - 1] - s[k])
    return s
