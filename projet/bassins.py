"""Méthode de Newton appliquée à beaucoup de points de départ du plan complexe à la fois."""

import numpy as np

from polynomes import evaluer, deriver


def grille(xmin, xmax, ymin, ymax, n):
    """Tableau n x n des points x + iy d'un rectangle du plan."""
    X, Y = np.meshgrid(np.linspace(xmin, xmax, n), np.linspace(ymin, ymax, n))
    return X + 1j * Y


def iterer(P, Z, tours):
    """Applique tours fois le pas de Newton à chaque point de Z."""
    dP = deriver(P)
    # un point peut tomber pile sur un zéro de P' : numpy donne inf ou nan, on laisse faire
    with np.errstate(all="ignore"):
        for _ in range(tours):
            Z = Z - evaluer(P, Z) / evaluer(dP, Z)
    return Z


def classer(P, racines, Z, tours=100, tol=1e-8):
    """Pour chaque point de départ : numéro de la racine atteinte et nombre de tours.

    Un point qui n'atteint aucune racine a le numéro -1 et le nombre de tours tours + 1.
    """
    classe = np.full(Z.shape, -1)
    nb = np.full(Z.shape, tours + 1)
    dP = deriver(P)
    with np.errstate(all="ignore"):
        for k in range(tours):
            Z = Z - evaluer(P, Z) / evaluer(dP, Z)
            for i, r in enumerate(racines):
                arrive = (classe == -1) & (np.abs(Z - r) < tol)
                classe[arrive] = i
                nb[arrive] = k + 1
    return classe, nb


def finit_sur(P, Z, cycle, tours=150, tol=1e-6):
    """Vrai pour les points de départ dont l'orbite finit sur un des points de cycle."""
    Z = iterer(P, Z, tours)
    res = np.zeros(Z.shape, dtype=bool)
    for z in cycle:
        res |= np.abs(Z - z) < tol
    return res


def couronne(r, n, rng, log=True):
    """n points au hasard dans la couronne r <= |z| < 1.5 r.

    log=True : rayon uniforme en échelle logarithmique. log=False : uniforme en surface.
    """
    if log:
        rho = r * 1.5 ** rng.random(n)
    else:
        rho = np.sqrt(r ** 2 + rng.random(n) * ((1.5 * r) ** 2 - r ** 2))
    return rho * np.exp(2j * np.pi * rng.random(n))


def bord(classe):
    """Pixels dont le voisin de droite ou du dessous n'a pas la même classe."""
    b = np.zeros(classe.shape, dtype=bool)
    b[:-1, :] |= classe[:-1, :] != classe[1:, :]
    b[:, :-1] |= classe[:, :-1] != classe[:, 1:]
    return b


def compter_boites(b, t):
    """Nombre de boîtes de t x t pixels qui contiennent au moins un pixel vrai de b."""
    n = b.shape[0] // t * t
    return int(b[:n, :n].reshape(n // t, t, n // t, t).any(axis=(1, 3)).sum())
