"""
Fourier coefficients of paramodular eigenforms for smf_qexp_coeffs and smf_qexp_reduction,
from Brandon Williams' eigenform files (Eigenforms_Weight{k}/level_{N}_weight_{k}_label_{i}.txt).

A file gives the Hecke field Q(b) of the eigenform and its Fourier expansion
    sum a(n, r, m) q^n r^r s^m + O(q, s)^prec,
where a(n, r, m) is the coefficient at the index [[N n, r/2], [r/2, m]], i.e. at the binary
quadratic form m x^2 + r x y + N n y^2.  For each file this script

  1. computes T(p) for the good p in (2, 3, 5) on the expansion (the action of
     Paramodular_Qexp/HeckeOperator), checks that the expansion is an eigenvector and reads off
     lambda_p in Q(b);
  2. finds the newform orbit of the db with the same dimension and an embedding of Q(b) into its
     Hecke field taking these lambda_p to the stored ones (exactly one orbit must match);
  3. writes every coefficient in the integral basis of the Hecke ring H of that orbit, normalized
     to unit content (the coefficients generate H; the convention of Poor-Yuen and
     Poor-Shurman-Yuen): scaled by a rational number so that the coordinates are coprime
     integers, and, if the ideal they generate is still proper, divided by a generator of it when
     it is principal and the quotients stay in H.  The scaling is then unique up to a unit of H;
     the first nonzero coordinate at the index of least discriminant is made positive.  The
     index of the ideal the coefficients generate is reported (1: unit content);
  4. reduces the indices with fourier_tables.legendre_reduce to (qf_legendre, qf_tag, det), and
     checks that det^k a(T) is constant on each class (a(u^T T u) = det(u)^k a(T)).

For each level it writes the classes (smf_qexp_reduction: the representative is the index of
least n + m in the class, then the least (m, r, n); classes sorted by discriminant) and for each
orbit the coefficient at the representative (smf_qexp_coeffs).  Nothing is written to the db.

    sage -python -m smf_lmfdb.db_tables.qexp_from_eigenforms OUTPREFIX FILE [FILE ...]
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict

from sage.all import (QQ, ZZ, PolynomialRing, NumberField, LaurentPolynomialRing, Matrix, vector,
                      sage_eval, lcm, gcd, prime_range, MatrixSpace, BinaryQF, pari)
from lmfdb import db

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fourier_tables import P1N_lex_minimize, legendre_reduce  # noqa: E402

CHECK_TAGS = os.environ.get('CHECK_TAGS') == '1'

INT_MAX = 2**31 - 1


# ---------------------------------------------------------------- the eigenform files

def parse(fname):
    N, k = map(int, re.search(r'level_(\d+)_weight_(\d+)_label_\d+\.txt$', fname).groups())
    lines = open(fname).read().split("\n")
    head = lines[0].strip()
    if head == "Hecke field: QQ":
        K = QQ
    else:
        K = NumberField(PolynomialRing(QQ, 'x')(head.replace("Hecke field: b = root of ", "")), 'b')
    ser = max(lines, key=len)
    m = re.search(r'\s*\+\s*O\(q, s\)\^(\d+)\s*$', ser)
    assert m, "no precision in %s" % fname
    prec = int(m.group(1))
    R = LaurentPolynomialRing(K, 'q,r,s')
    loc = dict(zip(['q', 'r', 's'], R.gens()))
    if K is not QQ:
        loc['b'] = K.gen()
    f = sage_eval(ser[:m.start()], locals=loc)
    # (A, B, C) = (N n, r, m): the index [[A, B/2], [B/2, C]]
    d = {(N * int(e[0]), int(e[1]), int(e[2])): K(c) for e, c in f.dict().items()}
    return N, k, K, d, prec


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


def hecke(p, d, N, k, prec):
    """{index: (T(p) f)(index)} for the indices where it is determined by the known coefficients."""
    fac = {'one': 1, 'p2k3': p**(2 * k - 3), 'pk2': p**(k - 2)}
    out = {}
    for m in d:
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
    Tp = hecke(p, d, N, k, prec)
    lams = {Tp[m] / d[m] for m in Tp if d[m] != 0}
    assert len(lams) == 1, "T(%d): not an eigenvector (%d ratios)" % (p, len(lams))
    return lams.pop(), len(Tp)


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


def process_file(fname, stats):
    N, k, K, d, prec = parse(fname)
    lam, nchk = {}, {}
    for p in (2, 3, 5):
        if N % p:
            lam[p], nchk[p] = eigenvalue(p, d, N, k, prec)
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
    divided = False
    if idx != 1:
        # the content ideal: if its O_K-ideal is principal and the quotients stay in H, divide
        # by a generator (unique up to a unit)
        els = [Kd(list(vector(QQ, v) * B)) for v in zspan(list(co.values()))]
        J = Kd.ideal(els)
        if J.is_principal():
            al = J.gens_reduced()[0]
            new = {m: list(vector(QQ, list(Kd(list(vector(QQ, v) * B)) / al)) * Bi) for m, v in co.items()}
            if all(x.denominator() == 1 for v in new.values() for x in v):
                co = {m: [ZZ(x) for x in v] for m, v in new.items()}
                if next(x for x in co[first] if x) < 0:
                    co = {m: [-x for x in v] for m, v in co.items()}
                idx, divided = content_index(list(co.values()), Kd, B), True
    info = {'file': os.path.basename(fname), 'label': f['label'], 'level': N, 'weight': k, 'dim': B.nrows(),
            'prec': prec, 'nonzero': len(co), 'lambda_checked_at': {str(p): n for p, n in nchk.items()},
            'content_index': int(idx), 'divided_by_generator': divided, 'max_abs': int(max(abs(x) for v in co.values() for x in v))}
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
