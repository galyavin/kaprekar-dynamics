"""Exact rational minimization by two-phase simplex with Bland's rule."""

from fractions import Fraction as F


def lp_min(c, A, b, Aeq=(), beq=()):
    """Minimize c.x subject to A x <= b and Aeq x = beq with free x.

    Return ('opt', value, x), ('infeasible',) or ('unbounded',).
    """
    n = len(c)
    rows = list(A) + list(Aeq)
    rhs = list(map(F, b)) + list(map(F, beq))
    m_ub = len(A)
    m = len(rows)
    # Split x into positive and negative parts, with slacks and artificials.
    nv = 2 * n + m_ub
    tableau = []
    for i, row in enumerate(rows):
        r = [F(0)] * (nv + m + 1)
        for j in range(n):
            r[j] = F(row[j])
            r[n + j] = -r[j]
        if i < m_ub:
            r[2 * n + i] = F(1)
        r[-1] = rhs[i]
        if r[-1] < 0:
            r = [-v for v in r]
        r[nv + i] = F(1)
        tableau.append(r)
    basis = [nv + i for i in range(m)]

    def pivot(pr, pc):
        value = tableau[pr][pc]
        tableau[pr] = [v / value for v in tableau[pr]]
        for i, row in enumerate(tableau):
            if i != pr and row[pc]:
                factor = row[pc]
                tableau[i] = [a - factor * b for a, b in zip(row, tableau[pr])]
        basis[pr] = pc

    def run(obj, allowed):
        while True:
            cb = [obj[i] for i in basis]
            enter = None
            for j in allowed:
                if j in basis:
                    continue
                cost = obj[j] - sum(cb[i] * tableau[i][j] for i in range(m) if cb[i])
                if cost < 0:
                    enter = j
                    break
            if enter is None:
                return "opt"
            best = None
            pr = None
            for i, row in enumerate(tableau):
                if row[enter] > 0:
                    ratio = row[-1] / row[enter]
                    if best is None or ratio < best or (ratio == best and basis[i] < basis[pr]):
                        best, pr = ratio, i
            if pr is None:
                return "unbounded"
            pivot(pr, enter)

    obj1 = [F(0)] * nv + [F(1)] * m
    run(obj1, range(nv + m))
    value = sum(tableau[i][-1] for i in range(m) if basis[i] >= nv)
    if value > 0:
        return ("infeasible",)
    for i in range(m):
        if basis[i] >= nv:
            for j in range(nv):
                if tableau[i][j]:
                    pivot(i, j)
                    break
    obj2 = [F(v) for v in c] + [-F(v) for v in c] + [F(0)] * (m_ub + m)
    if run(obj2, range(nv)) == "unbounded":
        return ("unbounded",)
    parts = [F(0)] * nv
    for i in range(m):
        if basis[i] < nv:
            parts[basis[i]] = tableau[i][-1]
    x = [parts[j] - parts[n + j] for j in range(n)]
    return ("opt", sum(F(c[j]) * x[j] for j in range(n)), x)
