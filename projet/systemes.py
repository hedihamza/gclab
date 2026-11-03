"""Matrices en listes de listes, pivot de Gauss (TP4) et droite des moindres carrés."""


def transposee(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def produit(A, B):
    """Produit matriciel A x B."""
    n, p, m = len(A), len(B), len(B[0])
    if len(A[0]) != p:
        raise ValueError("dimensions incompatibles")
    C = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(p):
                C[i][j] += A[i][k] * B[k][j]
    return C


def resoudre(A, b, tol=1e-14):
    """Résout A x = b par le pivot de Gauss avec pivot partiel. A et b ne sont pas modifiés."""
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[p][k]) <= tol:
            raise ValueError("système singulier")
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]
    return x


def moindres_carres(xs, ys):
    """Pente a et ordonnée à l'origine b de la droite y = a x + b la plus proche des points.

    Le système A (a, b) = y a une équation par point, donc trop d'équations : on résout
    à la place les équations normales tA A (a, b) = tA y, un système 2 x 2.
    """
    A = [[x, 1] for x in xs]
    tA = transposee(A)
    second = [ligne[0] for ligne in produit(tA, [[y] for y in ys])]
    a, b = resoudre(produit(tA, A), second)
    return a, b
