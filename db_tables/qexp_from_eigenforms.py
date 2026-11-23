"""
Fourier coefficients of paramodular eigenforms for smf_qexp_coeffs and smf_qexp_reduction,
from Brandon Williams' eigenform files (Eigenforms_Weight{k}/level_{N}_weight_{k}_label_{i}.txt).

A file gives the Hecke field Q(b) of the eigenform and its Fourier expansion
    sum a(n, r, m) q^n r^r s^m + O(q, s)^prec,
where a(n, r, m) is the coefficient at the index [[N n, r/2], [r/2, m]], i.e. at the binary
quadratic form m x^2 + r x y + N n y^2.  For each file this script

  1. computes T(p) for the good p in (2, 3, 5) on the expansion (the action of
     Paramodular_Qexp/HeckeOperator), at every index where it is determined, checks that the
     expansion is an eigenvector there and reads off lambda_p in Q(b).  Not every file is an
     eigenform (some are a basis of an eigenspace of T(2), or a list of expansions spanning an
     orbit): space_eigenforms takes the span of all the files of a space instead and
     diagonalizes T(2), T(3), T(5) on it;
  2. finds the newform orbit of the db with the same dimension and an embedding of Q(b) into its
     Hecke field taking these lambda_p to the stored ones (exactly one orbit must match);
  3. writes every coefficient in the integral basis of the Hecke ring H of that orbit,
     normalized: scaled by a rational number so that the coordinates are coprime integers, then
     so that the coefficients generate c O', where O' = {x : x I in I} is the multiplier ring of
     the H-module I they generate (an order containing H) and c the least positive integer with
     c O' in H (normalize_order).  If O' = H, this is unit content (the coefficients generate H;
     the convention of Poor-Yuen and Poor-Shurman-Yuen), which most orbits reach; in some O' is
     larger than H at 2, and then unit content is impossible.  The scaling is unique up to a
     unit of O'; the first nonzero coordinate at the index of least discriminant is made
     positive.  [H : I], [O' : H] and c are reported;
  4. reduces the indices with fourier_tables.legendre_reduce to (qf_legendre, qf_tag, det), and
     checks that det^k a(T) is constant on each class (a(u^T T u) = det(u)^k a(T)).

For each level it writes the classes (smf_qexp_reduction: the representative is the index of
least n + m in the class, then the least (m, r, n); classes sorted by discriminant) and for each
orbit the coefficient at the representative (smf_qexp_coeffs).  Nothing is written to the db.

    sage -python -m smf_lmfdb.db_tables.qexp_from_eigenforms OUTPREFIX FILE [FILE ...]
"""
import json
import os

import numpy
import re
import sys
from collections import Counter, defaultdict

from sage.all import (QQ, ZZ, prod, PolynomialRing, NumberField, LaurentPolynomialRing, Matrix, vector,
                      sage_eval, lcm, gcd, prime_range, MatrixSpace, BinaryQF, pari)
from lmfdb import db

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fourier_tables import P1N_lex_minimize, legendre_reduce  # noqa: E402

CHECK_TAGS = os.environ.get('CHECK_TAGS') == '1'

INT_MAX = 2**63 - 1      # smf_qexp_coeffs.coeff is bigint[]


# ---------------------------------------------------------------- the eigenform files

def parse(fname):
    N, k = map(int, re.search(r'level_(\d+)_weight_(\d+)_label_\d+(?:_wrong_label)?\.txt$', fname).groups())
    lines = open(fname).read().split("\n")
    head = lines[0].strip()
    if head == "Hecke field: QQ":
        K = QQ
    else:
        K = NumberField(PolynomialRing(QQ, 'x')(head.replace("Hecke field: b = root of ", "")), 'b')
    ser = max(lines, key=len).strip()
    # one expansion "f + O(q, s)^prec", or a list "[f1 + O(q, s)^prec, f2 + ...]" of expansions
    # spanning the orbit (with rational coefficients; the eigenform is a combination of them)
    if ser.startswith('[') and ser.endswith(']'):
        ser = ser[1:-1]
    parts = re.split(r'\s*\+\s*O\(q, s\)\^(\d+)\s*(?:,\s*|$)', ser)
    assert len(parts) >= 3 and parts[-1] == '', "no precision in %s" % fname
    precs = {int(x) for x in parts[1::2]}
    assert len(precs) == 1, "several precisions in %s" % fname
    prec = precs.pop()
    R = LaurentPolynomialRing(K, 'q,r,s')
    loc = dict(zip(['q', 'r', 's'], R.gens()))
    if K is not QQ:
        loc['b'] = K.gen()
    ds = []
    for src in parts[0::2][:-1]:
        f = sage_eval(src, locals=loc)
        # (A, B, C) = (N n, r, m): the index [[A, B/2], [B/2, C]]
        ds.append({(N * int(e[0]), int(e[1]), int(e[2])): K(c) for e, c in f.dict().items()})
    return N, k, K, ds, prec


# ---------------------------------------------------------------- T(p) (Paramodular_Qexp, fast_hecke.py)

def _div(t, q):
    if t[0] % q or t[1] % q or t[2] % q:
        return None
    return (t[0] // q, t[1] // q, t[2] // q)


def _images(m, p, N):
    A, B, C = m
    out = [((p * A, p * B, p * C), 'one'), (_div(m, p), 'p2k3')]
    out.append(((p * A, B, C // p) if C % p == 0 else None, 'pk2'))
    for j in range(p):
        jN = j * N
        a2 = A + B * jN + C * jN * jN
        out.append(((a2 // p, B + 2 * C * jN, C * p) if a2 % p == 0 else None, 'pk2'))
    return out


def box(N, prec):
    """The indices (N n, r, m) of the expansion: n, m >= 1, n + m < prec, 4 N n m > r^2."""
    out = []
    for n in range(1, prec):
        for mm in range(1, prec - n):
            R = int((4 * N * n * mm) ** 0.5) + 1
            out += [(N * n, r, mm) for r in range(-R, R + 1) if 4 * N * n * mm - r * r > 0]
    return out


def hecke(p, d, N, k, prec):
    """{index: (T(p) f)(index)} at every index of the box where it is determined by the known
    coefficients (those of the box; d holds the nonzero ones)."""
    fac = {'one': 1, 'p2k3': p**(2 * k - 3), 'pk2': p**(k - 2)}
    out = {}
    for m in box(N, prec):
        A, B, C = m
        val, ok = 0, True
        for img, f in _images(m, p, N):
            if img is None:
                continue
            if img[0] + N * img[2] >= N * prec:
                ok = False
                break
            val += fac[f] * d.get(img, 0)
        if ok:
            out[m] = val
    return out


def eigenvalue(p, d, N, k, prec):
    """lambda with T(p) f = lambda f at every index where T(p) f is determined, and the number of
    those indices; None if f vanishes at all of them (and so does T(p) f: no information)."""
    Tp = hecke(p, d, N, k, prec)
    nz = [m for m in Tp if d.get(m, 0) != 0]
    if not nz:
        assert all(v == 0 for v in Tp.values()), "T(%d): not an eigenvector (f = 0 where T(%d) f is known, T(%d) f != 0)" % (p, p, p)
        return None, len(Tp)
    lam = Tp[nz[0]] / d[nz[0]]
    bad = sum(1 for m in Tp if Tp[m] != lam * d.get(m, 0))
    assert not bad, "T(%d): not an eigenvector (%d of %d indices)" % (p, bad, len(Tp))
    return lam, len(Tp)


def eigen_from_span(ds, K, N, k, prec):
    """The eigenform in the span of the expansions ds (with coefficients in a subfield of K):
    T(p) on the span, for the first good p in (2, 3, 5) at which it has an eigenvalue in K of
    multiplicity one, then the eigenvector, with coordinates in K."""
    n = len(ds)
    for p in (2, 3, 5):
        if N % p == 0:
            continue
        Ts = [hecke(p, d, N, k, prec) for d in ds]
        S = sorted(Ts[0])
        A = Matrix(K, [[d.get(m, 0) for m in S] for d in ds])          # rows: the f_i on S
        if A.rank() < n:
            continue
        B = Matrix(K, [[T[m] for m in S] for T in Ts])                   # rows: the T(p) f_i on S
        M = A.solve_left(B)                                              # T(p) f_i = sum_j M_ij f_j
        assert M * A == B, "T(%d) does not preserve the span" % p
        cp = M.charpoly()
        for th, e in cp.roots(K):
            if e == 1:
                v = (M - th).left_kernel().basis()[0]
                return {m: c for m in set().union(*ds) for c in [sum(v[i] * ds[i].get(m, 0) for i in range(n))] if c != 0}
    raise AssertionError("no good p in (2, 3, 5) separates an eigenform in the span of %d expansions" % n)


def space_eigenforms(fnames):
    """The eigenforms of the space spanned by all the files of one level and weight, without
    trusting the files to be eigenforms: the rational span V of their expansions (each expansion
    over Q(b) gives its d coordinates in the power basis), the matrix of T(p) on V for the good p
    in (2, 3, 5) (with the rank of V on the indices where T(p) is determined: if it is dim V,
    T(p) f is determined by f in V, and V is checked to be T(p)-stable there), then a combination
    T of these with as many irreducible factors of multiplicity one in its characteristic
    polynomial as possible; each such factor gives an eigenform, the T-eigenvector with
    coordinates in the field of the factor.  Returns (diagnostics, [(K, d)])."""
    N = k = prec = None
    exps = []
    for fn in fnames:
        N, k, K, ds, prec1 = parse(fn)
        assert prec is None or prec1 == prec, "files of different precisions"
        prec = prec1
        for d in ds:
            if K is QQ:
                exps.append({m: QQ(c) for m, c in d.items()})
            else:
                for i in range(K.degree()):
                    e = {m: list(c)[i] for m, c in d.items() if list(c)[i] != 0}
                    if e:
                        exps.append(e)
    S0 = box(N, prec)
    A0 = Matrix(QQ, [[e.get(m, 0) for m in S0] for e in exps])
    basis = [exps[i] for i in A0.pivot_rows()]
    n = len(basis)
    diag = {'level': N, 'weight': k, 'prec': prec, 'files': len(fnames), 'expansions': len(exps), 'dim_span': n, 'rank': {}, 'stable': {}}
    Ms = {}
    for p in (2, 3, 5):
        if N % p == 0:
            continue
        Ts = [hecke(p, e, N, k, prec) for e in basis]
        S = sorted(Ts[0])
        A = Matrix(QQ, [[e.get(m, 0) for m in S] for e in basis])
        diag['rank'][p] = (int(A.rank()), len(S))
        if A.rank() < n:
            continue
        B = Matrix(QQ, [[T[m] for m in S] for T in Ts])
        M = A.solve_left(B)
        diag['stable'][p] = bool(M * A == B)
        if diag['stable'][p]:
            Ms[p] = M
    assert Ms, "T(p) is determined on the span for no good p in (2, 3, 5)"
    ps = sorted(Ms)
    best = None
    for cs in [(1,), (1, 1), (1, 2), (1, -1), (1, 3), (1, 1, 1), (1, 2, 3), (1, -2, 5)]:
        if len(cs) > len(ps):
            continue
        T = sum(c * Ms[p] for c, p in zip(cs, ps))
        fac = T.charpoly().factor()
        simple = sum(g.degree() for g, e in fac if e == 1)
        if best is None or simple > best[0]:
            best = (simple, cs, T, fac)
        if simple == n:
            break
    simple, cs, T, fac = best
    diag['combination'] = dict(zip(map(str, ps), cs))
    diag['charpoly_factors'] = sorted([int(g.degree()), int(e)] for g, e in fac)
    diag['charpoly_T2_squarefree'] = bool(2 in Ms and Ms[2].charpoly().is_squarefree())
    out = []
    for g, e in fac:
        if e != 1:
            continue
        Kg = QQ if g.degree() == 1 else NumberField(g, 'b')
        th = -g[0] / g[1] if g.degree() == 1 else Kg.gen()
        v = (T.change_ring(Kg) - th).left_kernel().basis()
        assert len(v) == 1
        v = v[0]
        f = {}
        for i, e2 in enumerate(basis):
            if v[i] != 0:
                for m, c in e2.items():
                    f[m] = f.get(m, 0) + v[i] * c
        out.append((Kg, {m: c for m, c in f.items() if c != 0}))
    return diag, out


# ---------------------------------------------------------------- the db orbit and the embedding

def ring_basis(h, Kd):
    if Kd is QQ:
        return Matrix(QQ, [[1]])
    return Matrix(QQ, [[QQ(c) / den for c in num] for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])])


def embedding(K, Kd, lam, lamd):
    """The embedding K -> Kd taking lam[p] to lamd[p] for every p, or None."""
    ps = sorted(lam)
    if K is QQ:
        return (lambda x: QQ(x)) if all(QQ(lam[p]) == lamd[p] for p in ps) else None
    if Kd is QQ or Kd.degree() != K.degree():
        return None
    d = K.degree()
    for t in range(-3, 4):
        x = lam[ps[0]] + (t * lam[ps[1]] if len(ps) > 1 else 0)
        y = lamd[ps[0]] + (t * lamd[ps[1]] if len(ps) > 1 else 0)
        M = Matrix(QQ, [list(x**i) for i in range(d)])
        if M.rank() < d:
            continue
        g = vector(QQ, [0, 1] + [0] * (d - 2)) * M.inverse()
        img = sum(c * y**i for i, c in enumerate(g))
        if K.polynomial()(img) != 0:
            return None
        s = K.hom([img], Kd)
        return s if all(s(lam[p]) == lamd[p] for p in ps) else None
    raise RuntimeError("no combination of the lambda_p generates the field")


def match_orbit(N, k, K, lam):
    deg = 1 if K is QQ else K.degree()
    hits = []
    for f in db.smf_newforms.search({'space_label': '2.K.%d.%d.0.a' % (N, k), 'dim': deg}, ['label', 'hecke_orbit_code']):
        h = db.smf_hecke_nf.lucky({'hecke_orbit_code': f['hecke_orbit_code']})
        if h is None:
            continue
        x = PolynomialRing(QQ, 'x')([QQ(c) for c in h['field_poly']])
        Kd = QQ if x.degree() == 1 else NumberField(x, 'a')
        B = ring_basis(h, Kd)
        el = (lambda v: QQ(v[0])) if Kd is QQ else (lambda v: Kd(list(vector(QQ, v) * B)))
        P = prime_range(100)
        lamd = {p: el(h['lambda_p'][P.index(p)]) for p in lam}
        s = embedding(K, Kd, lam, lamd)
        if s is not None:
            hits.append((f, h, Kd, B, s))
    assert len(hits) == 1, "level %d weight %d: %d db orbits match" % (N, k, len(hits))
    return hits[0]


# ---------------------------------------------------------------- reduction of the indices

def reduce_signed(abc, N):
    """fourier_tables.legendre_reduce for a positive definite index, with the determinant of the
    whole transformation: legendre_reduce returns det(u) for the reduction u, but the tag is then
    minimized over the orbit of the automorphisms of the reduced form, some of determinant -1
    (the TODO there).  Here the orbit is walked with the determinant of each automorphism, so
    that a(T) = sign^k a(class); sign is 0 when the tag is reached with both determinants (then
    a(T) = 0 in odd weight).  qf_legendre and qf_tag are those of legendre_reduce."""
    a, b, c = abc
    S = MatrixSpace(ZZ, 2, 2)
    qf = BinaryQF(a, b, N * c)
    pol = qf.reduced_form().polynomial()
    x, y = pol.parent().gens()
    ap, bp, cp = pol.coefficient(x**2), pol.coefficient(x * y), pol.coefficient(y**2)
    m1 = S([2 * a, b, b, 2 * N * c])
    m2 = S([2 * ap, bp, bp, 2 * cp])
    u = S(pari.qfisom(m2, m1))
    auts = [S(v) for v in pari.qfauto(m2)[1]]
    key = lambda v: (int(v[0, 0]), int(v[1, 0]))
    v0 = P1N_lex_minimize(S(u ** (-1)) * Matrix([[0], [1]]), N)
    dets = {key(v0): {1}}
    todo = [(v0, 1)]
    while todo:
        v, dv = todo.pop()
        for m in auts:
            w = P1N_lex_minimize(m * v, N)
            dw = dv * int(m.det())
            if dw not in dets.setdefault(key(w), set()):
                dets[key(w)].add(dw)
                todo.append((w, dw))
    tag = min(dets)
    ds = dets[tag]
    sign = 0 if len(ds) == 2 else int(u.det()) * ds.pop()
    return (int(ap), int(bp), int(cp)), tag, sign


# ---------------------------------------------------------------- coefficients

def coordinates(c, Kd, B, Bi, s):
    if Kd is QQ:
        return [QQ(s(c))]
    return list(vector(QQ, list(Kd(s(c)))) * Bi)


def zspan(vs):
    M = Matrix(QQ, vs)
    den = lcm([x.denominator() for x in M.list()] + [1])
    return [r / den for r in (M * den).change_ring(ZZ).echelon_form().rows() if r]


def colon(I, J):
    """A Z-basis of (J : I) = {x in K : x I in J}, for lattices I, J of full rank given by Z-bases
    (elements of K)."""
    K = I[0].parent()
    d = K.degree()
    MJi = Matrix(QQ, [list(y) for y in J]).inverse()
    W = [K.gen()**j for j in range(d)]
    A = Matrix(QQ, [sum([list(vector(QQ, list(w * y)) * MJi) for y in I], []) for w in W])   # x = t.W: x I in J iff t A integral
    den = lcm([c.denominator() for c in A.list()] + [1])
    Lam = (A * den).change_ring(ZZ).saturation()        # Z^n meet the row space of A
    return [sum(c * w for c, w in zip(A.solve_left(vector(QQ, l)), W)) for l in Lam.rows()]


def multiplier_ring(I, Bi=None):
    """A Z-basis of {x in K : x I in I}, for I given by a Z-basis (elements of K)."""
    return colon(I, I)


def order_basis(O, B):
    """The Z-basis of the order O (elements of K) that the coefficients are written in: 1 first, then
    the Hermite basis of O in the coordinates of the basis B of the Hecke ring H (rows: power-basis
    coordinates), taken from the last coordinate down, so that it is determined by O.  Since O
    contains H with small index, these are close to the basis of H, and coordinates stay small."""
    K = O[0].parent()
    d = K.degree()
    Bi = B.inverse()
    M = Matrix(QQ, [vector(QQ, list(y)) * Bi for y in O])
    den = lcm([x.denominator() for x in M.list()] + [1])
    E = (M * den).change_ring(ZZ).matrix_from_columns(list(range(d - 1, -1, -1))).echelon_form()
    rows = [K(list(vector(QQ, list(reversed(list(r)))) * B)) / den for r in E.rows()][::-1]
    assert rows[0] == 1, "1 is not the last Hermite vector: %s" % rows[0]
    return rows


def coords_in(basis):
    """y -> its coordinates in the Z-basis basis (elements of K)."""
    Mi = Matrix(QQ, [list(y) for y in basis]).inverse()
    return lambda y: vector(QQ, list(y)) * Mi


def least_norm_ideals(I, O, limit=2**20):
    """The integral ideals J = beta I (beta in K^*) of the order O of least index [O : J], with a beta
    for each: [(index, J as Z-basis, beta), ...] (one per J), or [] if none is found.  beta I lies in O
    iff beta lies in L = (O : I), and then a = beta I O_K is an integral ideal of O_K, in the class of
    I O_K, with [O : J] = N(a) [O : I] / N(I O_K).  So the ideals a of O_K of norm 1, 2, ... are run
    through; for each with a (I O_K)^{-1} = beta_0 O_K principal (PARI under GRH), the beta are beta_0 u,
    u a unit of O_K, with beta_0 u in L.  That depends on u modulo e O_K, e the least integer with
    e beta_0 O_K in L, and J on u modulo m O_K, m the least integer with m O_K in O (then u = u' mod m
    gives u/u' in O^*); the image of the units in (O_K / E)^*, E = lcm(e, m), is walked breadth first.
    Only the completeness of the list depends on GRH: each J found is checked exactly."""
    K = I[0].parent()
    d = K.degree()
    OK = K.maximal_order()
    ob = [K(x) for x in OK.basis()]
    okco = coords_in(ob)
    oco = coords_in(O)
    L = colon(I, O)
    lco = coords_in(L)
    inL = lambda y: all(x.denominator() == 1 for x in lco(y))
    IOK = K.ideal(I)
    idxI = abs(Matrix(QQ, [oco(y) for y in I]).det())            # [O : I] (I need not lie in O)
    m = lcm([x.denominator() for y in ob for x in oco(y)] + [1])
    U = K.unit_group(proof=False)
    gens = [K(g) for g in U.gens_values()]
    walks = {}

    def walk(E):
        """The image of the units in (O_K / E)^*: an integer array of residues (O_K-coordinates)
        and the exponents of a unit in each, breadth first."""
        if E not in walks:
            gm = [Matrix(ZZ, [okco(g * b) for b in ob]) for g in gens]     # multiplication by g on O_K-coordinates
            start = tuple(int(x) % E for x in okco(K(1)))
            seen = {start: (0,) * len(gens)}
            todo = [start]
            while todo:
                assert len(seen) < limit, "too many unit residues mod %d" % E
                r = todo.pop(0)
                v = vector(ZZ, r)
                for i, G in enumerate(gm):
                    r2 = tuple(int(x) % E for x in v * G)
                    if r2 not in seen:
                        ex = seen[r]
                        seen[r2] = ex[:i] + (ex[i] + 1,) + ex[i + 1:]
                        todo.append(r2)
            keys = list(seen)
            walks[E] = (numpy.array(keys, dtype=object), [seen[k] for k in keys])
        return walks[E]

    for n in range(1, 10**6):
        found = {}
        for a in K.ideals_of_bdd_norm(n)[n]:
            Q = a / IOK
            if not Q.is_principal(proof=False):
                continue
            b0 = Q.gens_reduced(proof=False)[0]
            assert K.ideal(b0) == Q
            e = 1
            while not all(inL(e * b0 * y) for y in ob):
                e += 1
            E = lcm(e, m)
            R, EX = walk(E)
            # b0 r in L, for r with O_K-coordinates x: x . C integral, C the matrix of y -> b0 y from
            # O_K-coordinates to L-coordinates; as integers: x . (D C) = 0 mod D
            C = Matrix(QQ, [lco(b0 * b) for b in ob])
            D = lcm([c.denominator() for c in C.list()] + [1])
            Ci = numpy.array([[int(c) for c in row] for row in (C * D).change_ring(ZZ).rows()], dtype=object)
            ok = numpy.all((R.dot(Ci)) % int(D) == 0, axis=1)
            for j in numpy.nonzero(ok)[0]:
                u = prod([g**x for g, x in zip(gens, EX[j])], K(1))
                assert inL(b0 * u)
                J = zspan([oco(b0 * u * y) for y in I])      # in O-coordinates
                assert all(x.denominator() == 1 for v in J for x in v)
                found.setdefault(tuple(tuple(v) for v in J), b0 * u)
        if found:
            idx = n * idxI / IOK.norm()
            out = []
            for J, beta in sorted(found.items()):
                assert abs(Matrix(QQ, J).det()) == idx
                out.append((idx, [vector(QQ, v) for v in J], beta))
            return out
    return []


def unit_search(I, O, beta, hco, HI, limit=2**20):
    """A unit u of O_K with I = beta u O', given I O_K = beta O_K, or None.  With L = beta^{-1} I
    (covolume that of O'), I = beta u O' iff u lies in L, and L must contain m O_K, m the exponent of
    O_K / O'; so only u mod m O_K matters, and the image of the unit group in (O_K / m)^* is walked
    breadth first (units of K under GRH; a u found is checked exactly, so only None depends on GRH)."""
    K = beta.parent()
    OK = K.maximal_order()
    ob = [K(x) for x in OK.basis()]
    Mob = Matrix(QQ, [list(x) for x in ob]); Mobi = Mob.inverse()
    okco = lambda y: vector(QQ, list(y)) * Mobi
    m = lcm([x.denominator() for y in ob for x in okco_rel(y, O)] + [1])
    L = [y / beta for y in I]
    ML = Matrix(QQ, [okco(y) for y in L])
    if any(x.denominator() != 1 for y in ob for x in (okco(m * y) * ML.inverse())):
        return None                       # L does not contain m O_K: I is not beta u O' for any unit u
    MLi = ML.inverse()
    inL = lambda r: all(x.denominator() == 1 for x in vector(QQ, r) * MLi)    # r: O_K-coordinates
    U = K.unit_group(proof=False)
    gens = [K(g) for g in U.gens_values()]
    key = lambda y: tuple(int(x) % m for x in okco(y))
    red = lambda y: sum(QQ(int(x) % m) * b for x, b in zip(okco(y), ob))
    start = K(1)
    seen = {key(start): []}
    todo = [(start, [])]
    while todo and len(seen) < limit:
        r, ex = todo.pop(0)
        if inL(okco(r)):
            u = prod([g**e for g, e in zip(gens, ex + [0] * (len(gens) - len(ex)))], K(1)) if ex else K(1)
            if zspan([hco(beta * u * y) for y in O]) == HI:
                return u
        for i, g in enumerate(gens):
            r2 = red(r * g)
            k2 = key(r2)
            if k2 not in seen:
                ex2 = (ex + [0] * len(gens))[:len(gens)]
                ex2[i] += 1
                seen[k2] = ex2
                todo.append((r2, ex2))
    return None


def okco_rel(y, O):
    """Coordinates of y in the Z-basis O (to find m with m O_K in O')."""
    return vector(QQ, list(y)) * Matrix(QQ, [list(z) for z in O]).inverse()


def basis_rows(O):
    """[numerators, denominators] of the Z-basis O (elements of K), in the power basis, in the format
    of hecke_ring_numerators / hecke_ring_denominators."""
    nums, dens = [], []
    for y in O:
        den = lcm([x.denominator() for x in list(y)] + [1])
        nums.append([int(x * den) for x in list(y)])
        dens.append(int(den))
    return [nums, dens]


def normalize_order(co, Kd, B, Bi):
    """The coefficients (H-coordinates co, coprime integers) generate a Z-module I, which is a module
    over H and over its multiplier ring O' = {x : x I in I}, an order containing H.  Write them in a
    basis of O' (that of H if O' = H, order_basis otherwise) and scale them:
    * if I = alpha O' (alpha searched among short elements of I, the generator of I O_K and its unit
      multiples), to a(T) / alpha, which generate O' (unit content in O'), unique up to a unit of O';
    * otherwise to beta a(T), which generate the integral ideal J = beta I of O' of least index in its
      class (least_norm_ideals; the least in Hermite form if there are several), unique up to a unit
      of O'.
    Returns the new coordinates and what was done."""
    elt = lambda v: Kd(list(vector(QQ, v) * B))
    hco = lambda y: vector(QQ, list(y)) * Bi
    S = zspan(list(co.values()))
    I = [elt(v) for v in zspan([hco(elt(v) * elt(w)) for v in S for w in Matrix.identity(QQ, B.nrows()).rows()])]   # the H-module
    O = multiplier_ring(I)
    mult = abs(1 / Matrix(QQ, [hco(y) for y in O]).det())       # [O' : H]
    if mult == 1:
        O = [elt(r) for r in Matrix.identity(QQ, B.nrows()).rows()]
        oco = hco
        order = None
    else:
        O = order_basis(O, B)
        oco = coords_in(O)
        order = basis_rows(O)
    covol = abs(Matrix(QQ, [oco(y) for y in I]).det())          # [O' : I] = |N(alpha)| if I = alpha O'
    OI = zspan([oco(y) for y in I])
    cands = list(I)
    LL = Matrix(QQ, [hco(y) for y in I])
    den = lcm([x.denominator() for x in LL.list()] + [1])
    short = [elt(r / den) for r in (LL * den).change_ring(ZZ).LLL().rows()]
    cands += short + [x + s_ * y for i, x in enumerate(short) for y in short[i + 1:] for s_ in (1, -1)]
    gen = ideal_generator(Kd.ideal(I), I)
    if gen is not None:
        cands.append(gen)
        u = unit_search(I, O, gen, oco, OI)
        if u is not None:
            cands.append(gen * u)
    info = {'multiplier_index': int(mult), 'order': order}
    for al in cands:
        if al == 0 or abs(al.norm()) != covol:
            continue
        if zspan([oco(al * y) for y in O]) == OI:
            new = {m: [ZZ(x) for x in oco(elt(v) / al)] for m, v in co.items()}
            info.update({'generator_found': True, 'ideal_index': 1, 'ideal': None, 'least_ideals': 1})
            return new, info
    res = least_norm_ideals(I, O)
    assert res, "no integral ideal found in the class of I"
    idx, J, beta = res[0]
    new = {m: [ZZ(x) for x in oco(beta * elt(v))] for m, v in co.items()}
    assert zspan(list(new.values()) + [list(oco(y * z)) for y in O for z in [beta * elt(v) for v in S]]) == J
    info.update({'generator_found': False, 'ideal_index': int(idx), 'ideal': [[int(x) for x in v] for v in J], 'least_ideals': len(res)})
    return new, info


def content_index(vecs, Kd, B):
    """[H : I] for the ideal I of the Hecke ring H generated by the coefficients (coordinates vecs)."""
    if Kd is QQ:
        return 1
    S = zspan(vecs)
    Bi = B.inverse()
    elt = lambda v: Kd(list(vector(QQ, v) * B))
    prods = [list(vector(QQ, list(elt(v) * elt(w))) * Bi) for v in S for w in Matrix.identity(QQ, B.nrows()).rows()]
    I = zspan(prods)
    assert len(I) == B.nrows()
    return abs(Matrix(QQ, I).det())


def ideal_generator(J, els):
    """A generator of the ideal J of O_K, or None if none is found.  First an element of norm N(J)
    among the given generators of J and their sums and differences of two; then PARI's class group
    under GRH (Sage's default certifies it, which takes hours in degree 10-15).  A generator found
    either way is checked to generate J, so only a None depends on GRH."""
    nJ = J.norm()
    cands = list(els) + [x + s * y for i, x in enumerate(els) for y in els[i + 1:] for s in (1, -1)]
    for x in cands:
        if x != 0 and abs(x.norm()) == nJ and J.number_field().ideal(x) == J:
            return x
    if not J.is_principal(proof=False):
        return None
    al = J.gens_reduced(proof=False)[0]
    assert J.number_field().ideal(al) == J
    return al


def process_file(fname, stats):
    N, k, K, ds, prec = parse(fname)
    d = ds[0] if len(ds) == 1 else eigen_from_span(ds, K, N, k, prec)
    return normalize(N, k, K, d, prec, os.path.basename(fname))


def normalize(N, k, K, d, prec, name):
    """Check that the expansion d (coefficients in K) is an eigenvector, match it to a db orbit,
    and write it in that orbit's Hecke ring basis, normalized; see the module docstring."""
    lam, nchk = {}, {}
    for p in (2, 3, 5):
        if N % p:
            x, nchk[p] = eigenvalue(p, d, N, k, prec)
            if x is not None:
                lam[p] = x
    assert lam, "no good p in (2, 3, 5) gives an eigenvalue"
    f, h, Kd, B, s = match_orbit(N, k, K, lam)
    Bi = B.inverse()
    raw = {m: coordinates(c, Kd, B, Bi, s) for m, c in d.items() if c != 0}
    L = lcm([x.denominator() for v in raw.values() for x in v])
    G = gcd([ZZ(x * L) for v in raw.values() for x in v])
    scale = QQ(L) / G
    co = {m: [ZZ(x * scale) for x in v] for m, v in raw.items()}
    first = min(co, key=lambda m: (4 * m[0] * m[2] - m[1]**2, m))
    if next(x for x in co[first] if x) < 0:
        co = {m: [-x for x in v] for m, v in co.items()}
    idx = content_index(list(co.values()), Kd, B)
    order = {'multiplier_index': 1, 'order': None, 'generator_found': True, 'ideal_index': 1, 'ideal': None, 'least_ideals': 1}
    if idx != 1:
        co, order = normalize_order(co, Kd, B, Bi)       # now in O'-coordinates
        if next(x for x in co[first] if x) < 0:
            co = {m: [-x for x in v] for m, v in co.items()}
    info = {'file': name, 'label': f['label'], 'level': N, 'weight': k, 'dim': B.nrows(),
            'prec': prec, 'nonzero': len(co), 'lambda_checked_at': {str(p): n for p, n in nchk.items()},
            'content_index': int(idx), 'max_abs': int(max(abs(x) for v in co.values() for x in v))}
    info.update(order)
    return N, k, f, co, prec, info


def classes(N, prec, k, coeffs_by_orbit):
    """The classes of the indices in the precision box, their representatives, and the
    coefficient of each orbit at each representative (checking det^k a(T) on each class)."""
    cls = {}
    for n in range(1, prec):
        for mm in range(1, prec - n):
            rmax = int((4 * N * n * mm) ** 0.5) + 1
            for r in range(-rmax, rmax + 1):
                if 4 * N * n * mm - r * r <= 0:
                    continue
                qf, tag, det = reduce_signed((mm, r, n), N)
                if CHECK_TAGS:
                    q0, t0, _ = legendre_reduce((mm, r, n), N)
                    assert (tuple(int(x) for x in q0), (int(t0[0]), int(t0[1]))) == (qf, tag), ((mm, r, n), q0, t0, qf, tag)
                if det == 0 and k % 2 == 0:
                    det = 1
                key = (qf, tag)
                cls.setdefault(key, []).append(((N * n, r, mm), int(det)))
    reps, values = {}, defaultdict(dict)
    bad = Counter()
    for key, members in cls.items():
        rep, rdet = min(members, key=lambda t: (t[0][0] // N + t[0][2], t[0][2], t[0][1], t[0][0]))
        reps[key] = rep
        for code, co in coeffs_by_orbit.items():
            dim = len(next(iter(co.values())))
            # a(T) = det^k a(class); det 0: a(T) = 0 if k is odd
            seen = {tuple((det**k) * x for x in co.get(m, [0] * dim)) for m, det in members if det}
            zero = all(not any(co.get(m, [])) for m, det in members if not det) or k % 2 == 0
            if len(seen) > 1 or not zero:
                bad[code] += 1
                continue
            v = seen.pop() if seen else (0,) * dim
            values[code][key] = [(rdet**k if rdet else 1) * x for x in v]
    return cls, reps, values, bad


def main(prefix, files):
    stats, infos = Counter(), []
    by_level = defaultdict(dict)
    for fn in files:
        N, k, f, co, prec, info = process_file(fn, stats)
        print(json.dumps(info), flush=True)
        infos.append(info)
        by_level[(N, k, prec)][f['hecke_orbit_code']] = co
    red_lines, coeff_lines, summary = [], [], []
    for (N, k, prec), orbits in sorted(by_level.items()):
        cls, reps, values, bad = classes(N, prec, k, orbits)
        disc = lambda key: 4 * reps[key][0] * reps[key][2] - reps[key][1]**2   # -disc of m x^2 + r x y + N n y^2
        keys = sorted(reps, key=lambda key: (disc(key), key))
        for i, key in enumerate(keys):
            A, r, m = reps[key]
            red_lines.append({'level': N, 'family': 'K', 'index': i, 'is_minimal': True, 'qf_legendre': list(key[0]),
                              'qf_tag': list(key[1]), 'qf_orbit_rep': [m, r, A]})
        for code, vals in values.items():
            for key in keys:
                if key in vals:
                    coeff_lines.append({'hecke_orbit_code': code, 'qf_legendre': list(key[0]), 'qf_tag': list(key[1]), 'coeff': [int(x) for x in vals[key]]})
        summary.append({'level': N, 'weight': k, 'classes': len(keys), 'sign_inconsistent': {str(c): n for c, n in bad.items()},
                        'overflow': sum(1 for l in coeff_lines if any(abs(x) > INT_MAX for x in l['coeff']))})
        print(json.dumps(summary[-1]), flush=True)
    json.dump({'files': infos, 'levels': summary, 'reduction': red_lines, 'coeffs': coeff_lines}, open(prefix + '.json', 'w'))
    print("wrote %s.json: %d reduction rows, %d coefficient rows" % (prefix, len(red_lines), len(coeff_lines)))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
