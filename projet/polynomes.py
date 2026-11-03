"""Polynômes donnés par la liste de leurs coefficients, du degré 0 au degré n (comme au TP3)."""


def evaluer(P, x):
    """Évalue P en x par le schéma de Horner. x peut être complexe, ou un tableau numpy."""
    s = 0
    for c in reversed(P):
        s = s * x + c
    return s


def deriver(P):
    """Coefficients du polynôme dérivé de P."""
    return [i * P[i] for i in range(1, len(P))]


def diviser(P, r):
    """Quotient de P par (x - r), par Horner. Le reste, qui vaut P(r), est jeté."""
    Q = [0] * (len(P) - 1)
    s = 0
    for i in range(len(P) - 1, 0, -1):
        s = s * r + P[i]
        Q[i - 1] = s
    return Q


def newton(P, x, tol=1e-14, tours=100):
    """Racine de P par Newton depuis x. Renvoie (racine, nombre de tours)."""
    dP = deriver(P)
    for k in range(tours):
        d = evaluer(dP, x)
        if d == 0:
            raise ZeroDivisionError("tangente horizontale")
        y = x - evaluer(P, x) / d
        if abs(y - x) <= tol * max(1.0, abs(y)):
            return y, k + 1
        x = y
    raise RuntimeError(f"pas de convergence en {tours} tours")


def racines(P, depart=0.5 + 0.5j):
    """Toutes les racines complexes de P : Newton, division par (x - r), et on recommence.

    Le départ n'est pas réel : avec des coefficients réels, Newton parti d'un réel reste
    sur l'axe réel et ne trouverait jamais une racine complexe. Chaque racine est repolie
    sur P lui-même, pour effacer l'erreur accumulée par les divisions.
    """
    Q = list(P)
    res = []
    while len(Q) > 2:
        r, _ = newton(Q, depart)
        res.append(newton(P, r)[0])
        Q = diviser(Q, r)
    res.append(newton(P, -Q[0] / Q[1])[0])
    return res


def newton_securise(P, a, b, x, tol=1e-12, tours=200):
    """Racine réelle de P dans [a, b], où P(a) et P(b) sont de signes opposés.

    Pas de Newton s'il tombe dans l'intervalle, pas de dichotomie sinon. L'intervalle
    qui encadre la racine rétrécit à chaque tour, donc on ne peut pas tourner en rond.
    """
    dP = deriver(P)
    fa = evaluer(P, a)
    if fa * evaluer(P, b) > 0:
        raise ValueError("P(a) et P(b) doivent être de signes opposés")
    for k in range(tours):
        fx = evaluer(P, x)
        if fx == 0:
            return x, k
        if fa * fx < 0:
            b = x
        else:
            a, fa = x, fx
        d = evaluer(dP, x)
        y = x - fx / d if d != 0 else a
        if not a < y < b:
            y = (a + b) / 2
        if abs(y - x) <= tol * max(1.0, abs(y)):
            return y, k + 1
        x = y
    raise RuntimeError(f"pas de convergence en {tours} tours")
