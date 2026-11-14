"""
Correct the Klingen-Eisenstein series (aut_rep_type 'Q', level 1) and the
Yoshida lifts (aut_rep_type 'Y', level 2) in smf_newforms, smf_newspaces and
smf_hecke_nf.

At each embedding, the Klingen-Eisenstein series of a classical newform f of
level 1 and weight j+k has the eigenvalues KE_* of
Hecke_Eigenvalues_Klingen_Eisenstein in terms of a_p(f), a_{p^2}(f), and spinor
L-function L(f, s) L(f, s-k+2); the Yoshida lift of a pair (f, g) of level-2
newforms of weights j+2k-2 and j+2 has the eigenvalues Y_* of
Hecke_Eigenvalues_Yoshida in terms of a_p(f), a_p(g), and L-function
L(f, s) L(g, s-k+2).  (Both sets of formulas agree with these L-functions: the
Euler factor is 1 - l t + p (l_1 + (p^2+1) l_2) t^2 - ..., l_i = lambda_p_square_i.)
The stored values were computed

  * Q, dimension > 1: by applying the affine KE_*_square formulas once to the
    orbit's Tr(a_{p^2}), so the constant term is counted once instead of dim
    times; in smf_hecke_nf, by applying them to each coordinate, so the constant
    is added to every coordinate instead of to the coordinate of 1;
  * Y, dimension > 1: by applying the Y_* formulas to Tr(a_p(f)) and Tr(a_p(g)),
    which drops the factor dim(g) (resp. dim(f)) and treats Tr(a)^2 as Tr(a^2);
  * the newspace sums the same way, once per space;
  * `traces` (all Q and Y rows) and smf_hecke_nf.an (Q of dimension > 1, all Y)
    from the old Euler factor, which reads lambda_p_square as T(p,2), and for an
    orbit as a product of traces; for Y with values (0, -4, ...) at even n.

This script recomputes, from the classical data (mf_newforms.traces for forms of
dimension 1, mf_hecke_nf.ap otherwise, with the a_n built by the Hecke
recursion and checked against mf_hecke_nf.an),

  smf_newforms:  trace_lambda_p, trace_lambda_p_square[_0,_1,_2], trace_display
                 (rows of dimension > 1), traces (all rows)
  smf_newspaces: eis_Q_lambda_p*, cusp_Y_lambda_p* (sums over the orbits, in
                 the spaces that have Q or Y rows)
  smf_hecke_nf:  lambda_p_square[_0,_1,_2] (Q), an (Q and Y)

with `traces` and `an` NULL at every n not prime to the level, as for the other
types.  The Yoshida rows of weight k = 2 keep their stored square columns: by a
deliberate convention of the generating code they hold p times lambda_p_square_0
and _1 (p^2 times _2 when j = 0); only their traces and an change.

Before writing anything it checks that each row's classical forms have the
expected level, weight, trivial character and dimensions, that the classical
a_n built from a_p agree with the stored mf_hecke_nf.an, that each stored
column is exactly what the old computation gives (so that nothing else is
being overwritten), that the stored smf_hecke_nf eigenvalues that do not change
are the formulas at the embedding matched by lambda_p, and that the new
Dirichlet coefficients satisfy the Euler factor at the good p <= 13.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_qy_traces [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_qy_traces [--datafile PREFIX] --apply
"""
import argparse
from collections import defaultdict, Counter

from sage.all import Matrix, NumberField, PolynomialRing, QQ, ZZ, divisors, factor, prime_range, vector
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import as_int, write_update_file

HT = ['lambda_p', 'lambda_p_square', 'lambda_p_square_0', 'lambda_p_square_1', 'lambda_p_square_2']
SQUARES = HT[1:]
NF_COLS = ['id', 'label', 'level', 'weight', 'dim', 'related_objects', 'traces', 'trace_display'] + ['trace_' + h for h in HT]
HNF_COLS = ['id', 'label', 'an', 'maxp', 'maxp_square', 'field_poly', 'hecke_ring_power_basis',
            'hecke_ring_numerators', 'hecke_ring_denominators'] + HT
RING_COLS = ['field_poly', 'hecke_ring_power_basis', 'hecke_ring_numerators', 'hecke_ring_denominators']
SPACE_PREFIX = {'Q': 'eis_Q_', 'Y': 'cusp_Y_'}
NMAX = 200   # length of smf_newforms.traces


def KE(k, j, p, a, a2):
    p = QQ(p)
    return {'lambda_p': (1 + p**(k-2)) * a,
            'lambda_p_square': (1 + p**(k-2) + p**(2*k-4)) * a2 + p**(2*k+j-4) * (p-1),
            'lambda_p_square_0': (1 + (p-1) * p**(k-3) + p**(2*k-4)) * a2 + p**(k+j-2) * (p**(k-2) * (p-1) - p**(2*k-4) - 1),
            'lambda_p_square_1': p**(k-3) * a2 + p**(k+j-2) * (1 - p**(k-4) + p**(2*k-4)),
            'lambda_p_square_2': p**(2*k+j-6) + 0 * a}


def KE_constant(k, j, p, h):
    """The term of KE_h not involving a_{p^2}."""
    return KE(k, j, p, QQ(0), QQ(0))[h]


def YO(k, j, p, a, b):
    p = QQ(p)
    return {'lambda_p': a + p**(k-2) * b,
            'lambda_p_square': a**2 + p**(k-2) * a * b + p**(2*k-4) * b**2 - p**(2*k+j-4) * (1 + 2*p),
            'lambda_p_square_0': a**2 + p**(2*k-4) * b**2 + a * b * p**(k-3) * (p-1) - 2 * p**(2*k+j-4) * (1+p),
            'lambda_p_square_1': a * b * p**(k-3) + (p**2-1) * p**(2*k+j-6),
            'lambda_p_square_2': p**(2*k+j-6) + 0 * a}


def k2_scale(j, p, h):
    """Factor by which the generating code scales the Yoshida square columns at k = 2."""
    if h in ('lambda_p_square_0', 'lambda_p_square_1'):
        return p
    if h == 'lambda_p_square_2' and j == 0:
        return p**2
    return 1


def ring_basis(h):
    K = NumberField(PolynomialRing(QQ, 'x')([int(c) for c in h['field_poly']]), 'nu')
    nu = K.gen()
    if h['hecke_ring_power_basis']:
        return K, [nu**i for i in range(K.degree())]
    return K, [sum(int(c) * nu**i for i, c in enumerate(num)) / int(den)
               for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])]


def coords(x, K, B):
    """Coordinates of x in the basis B of K, as integers."""
    M = Matrix(QQ, [b.list() for b in B])
    c = vector(QQ, K(x).list()) * M.inverse()
    assert all(y in ZZ for y in c), "non-integral coordinates %s" % c
    return [int(y) for y in c]


class Classical:
    """A classical newform orbit: its field and the a_n, 1 <= n <= NMAX, at one embedding."""
    def __init__(self, label):
        f = db.mf_newforms.lookup(label, ['label', 'level', 'weight', 'dim', 'char_order', 'traces', 'hecke_orbit_code'])
        assert f is not None, label
        self.label, self.level, self.weight, self.dim = label, f['level'], f['weight'], f['dim']
        assert f['char_order'] == 1, "%s: character is not trivial" % label
        if self.dim == 1:
            self.K, self.an = QQ, [QQ(x) for x in f['traces'][:NMAX]]
            return
        h = db.mf_hecke_nf.lucky({'hecke_orbit_code': f['hecke_orbit_code']}, ['an', 'ap', 'maxp'] + RING_COLS)
        assert h is not None, "%s: no mf_hecke_nf row" % label
        self.K, B = ring_basis(h)
        assert self.K.degree() == self.dim, label
        elt = lambda c: sum(QQ(x) * b for x, b in zip(c, B))
        ap = dict(zip(prime_range(h['maxp'] + 1), [elt(c) for c in h['ap']]))
        an = [None] * (NMAX + 1)
        an[1] = self.K(1)
        for n in range(2, NMAX + 1):
            fa = factor(n)
            if len(fa) > 1:
                an[n] = 1
                for p, e in fa:
                    an[n] *= an[p**e]
                continue
            p, e = fa[0]
            if e == 1:
                an[n] = ap[p]
            else:
                chi_p = 0 if self.level % p == 0 else 1
                an[n] = ap[p] * an[p**(e-1)] - chi_p * p**(self.weight - 1) * an[p**(e-2)]
        self.an = an[1:]
        stored = [elt(c) for c in h['an']]
        assert self.an[:len(stored)] == stored, "%s: a_n from the Hecke recursion differ from mf_hecke_nf.an" % label
        assert [ZZ(x.trace()) for x in self.an] == [ZZ(x) for x in f['traces'][:NMAX]], "%s: traces" % label

    def a(self, n):
        return self.an[n - 1]


def tr(x, L):
    return QQ(x) if L is QQ else QQ(L(x).trace())


def lift(r, cl):
    """(L, A, B, info): the field of the lift and the a_n of the two classical forms in it."""
    k, j = r['weight']
    labels = ['.'.join(x.split('/')[-4:]) for x in r['related_objects']]
    t = r['aut_rep_type']
    if t == 'Q':
        assert len(labels) == 1, r['label']
        f = cl(labels[0])
        assert (f.level, f.weight) == (1, j + k), "%s: classical form %s" % (r['label'], f.label)
        assert f.dim == r['dim'], r['label']
        return f.K, f.an, f.an, (f, f)
    assert len(labels) == 2, r['label']
    f, g = cl(labels[0]), cl(labels[1])
    assert (f.level, f.weight, g.level, g.weight) == (2, j + 2*k - 2, 2, j + 2), \
        "%s: classical forms %s, %s" % (r['label'], f.label, g.label)
    assert f.dim * g.dim == r['dim'] and min(f.dim, g.dim) == 1, r['label']
    if g.dim == 1:
        return f.K, f.an, [f.K(x) for x in g.an] if f.K is not QQ else g.an, (f, g)
    return g.K, [g.K(x) for x in f.an], g.an, (f, g)


def eigenvalues(r, L, A, B, p):
    k, j = r['weight']
    if r['aut_rep_type'] == 'Q':
        return KE(k, j, p, A[p-1], A[p*p-1] if p*p <= NMAX else None) if p*p <= NMAX else {'lambda_p': (1 + QQ(p)**(k-2)) * A[p-1]}
    return YO(k, j, p, A[p-1], B[p-1])


def dirichlet(r, A, B):
    """a_n of L(f,s) L(g,s-k+2) at one embedding, None at n not prime to the level."""
    k = r['weight'][0]
    N = r['level']
    return [None if ZZ(n).gcd(N) != 1 else sum(A[d-1] * B[n//d - 1] * QQ(n//d)**(k-2) for d in divisors(n))
            for n in range(1, NMAX + 1)]


def check_euler(r, ev, c):
    """c_{p^2} = l_p^2 - p (l_1 + (p^2+1) l_2) and c_p = l_p, at the good p <= 13."""
    for p in prime_range(14):
        if r['level'] % p == 0:
            continue
        e = ev[p]
        assert c[p-1] == e['lambda_p'], r['label']
        assert c[p*p-1] == e['lambda_p']**2 - p * (e['lambda_p_square_1'] + (p**2 + 1) * e['lambda_p_square_2']), \
            "%s: the Dirichlet coefficient a_%d does not satisfy the Euler factor" % (r['label'], p*p)


def old_trace_columns(r, L, info, ev, h, ps):
    """What the generating code stored in smf_newforms.trace_h at the primes ps."""
    k, j = r['weight']
    t, d = r['aut_rep_type'], r['dim']
    if t == 'Q':
        new = [tr(ev[p][h], L) for p in ps]
        return new if h == 'lambda_p' else [x - (d - 1) * KE_constant(k, j, p, h) for x, p in zip(new, ps)]
    f, g = info
    old = [YO(k, j, p, QQ(f.a(p).trace() if f.dim > 1 else f.a(p)), QQ(g.a(p).trace() if g.dim > 1 else g.a(p)))[h] for p in ps]
    if k == 2:
        old = [x * k2_scale(j, p, h) for x, p in zip(old, ps)]
    return old


def main(apply=False, datafile=None):
    rows = list(db.smf_newforms.search({'aut_rep_type': {'$in': ['Q', 'Y']}}, NF_COLS + ['aut_rep_type']))
    hnf = {h['label']: h for h in db.smf_hecke_nf.search({'label': {'$in': [r['label'] for r in rows]}}, HNF_COLS)}
    print("%d Q and Y newforms, %d with smf_hecke_nf rows" % (len(rows), len(hnf)))
    cache = {}
    def cl(label):
        if label not in cache:
            cache[label] = Classical(label)
        return cache[label]

    nf_rows, hnf_rows, stats = [], [], Counter()
    new_vals = {}
    for r in rows:
        t, (k, j), d, N = r['aut_rep_type'], r['weight'], r['dim'], r['level']
        L, A, B, info = lift(r, cl)
        ev = {p: eigenvalues(r, L, A, B, p) for p in prime_range(NMAX)}
        c = dirichlet(r, A, B)
        check_euler(r, ev, c)
        out = {'id': r['id']}
        # trace columns
        for h in HT:
            stored = [as_int(x) for x in r['trace_' + h]]
            ps = prime_range(NMAX)[:len(stored)]
            assert old_trace_columns(r, L, info, ev, h, ps) == stored, \
                "%s: stored trace_%s is not what the old computation gives" % (r['label'], h)
            if t == 'Y' and k == 2 and h in SQUARES:
                new = stored    # the k = 2 convention, kept
            else:
                new = [int(tr(ev[p][h], L)) for p in ps]
            out['trace_' + h] = new
        out['trace_display'] = out['trace_lambda_p'][:len(r['trace_display'])]
        assert r['trace_display'] == [as_int(x) for x in r['trace_lambda_p'][:len(r['trace_display'])]], r['label']
        # Dirichlet coefficients
        stored = [as_int(x) for x in r['traces']]
        assert len(stored) == NMAX and stored[0] == 1, r['label']
        for i, p in enumerate(prime_range(NMAX)):
            if N % p:
                assert stored[p-1] == as_int(r['trace_lambda_p'][i]), "%s: stored traces at %d" % (r['label'], p)
        out['traces'] = [None if x is None else int(tr(x, L)) for x in c]
        new_vals[r['label']] = out
        if any(out[col] != [as_int(x) for x in r[col]] for col in out if col != 'id'):
            nf_rows.append(out)
            stats['%s dim%s newforms change' % (t, '1' if d == 1 else '>1')] += 1

        # smf_hecke_nf
        h = hnf.get(r['label'])
        if h is None:
            stats['%s dim%s without smf_hecke_nf' % (t, '1' if d == 1 else '>1')] += 1
            continue
        K2, B2 = ring_basis(h)
        assert K2.degree() == d, r['label']
        elt = lambda cc: sum(QQ(x) * b for x, b in zip(cc, B2))
        lam = dict(zip(prime_range(h['maxp'] + 1), h['lambda_p']))
        # the embedding of L into K2 matching the stored lambda_p (which are right)
        embs = [lambda x: K2(QQ(x))] if L is QQ else L.embeddings(K2)
        embs = [phi for phi in embs if all(phi(ev[p]['lambda_p']) == elt(lam[p]) for p in lam)]
        assert embs, "%s: no embedding matches the stored lambda_p" % r['label']
        phi = embs[0]
        hrow = {'id': h['id']}
        for hname in SQUARES:
            ps = prime_range(h['maxp_square'] + 1)
            stored_c = [[int(x) for x in cc] for cc in h[hname]]
            if t == 'Q' and d > 1:
                f = info[0]
                Bf = ring_basis(db.mf_hecke_nf.lucky({'label': f.label}, RING_COLS))[1]
                a2c = [coords(f.a(p*p), f.K, Bf) for p in ps]
                # the affine KE formula applied to each coordinate of a_{p^2}
                old = [[int(KE(k, j, p, QQ(x), QQ(x))[hname]) for x in a2c[i]] for i, p in enumerate(ps)]
                assert old == stored_c, "%s: stored smf_hecke_nf %s is not the old computation" % (r['label'], hname)
                hrow[hname] = [coords(phi(ev[p][hname]), K2, B2) for p in ps]
            else:
                scale = (lambda p: k2_scale(j, p, hname)) if (t == 'Y' and k == 2) else (lambda p: 1)
                assert all(elt(cc) == phi(ev[p][hname] * scale(p)) for p, cc in zip(ps, stored_c)), \
                    "%s: stored smf_hecke_nf %s is not the formula" % (r['label'], hname)
                hrow[hname] = stored_c
        stored_an = h['an']
        hrow['an'] = ['NULL' if x is None else coords(phi(x), K2, B2) for x in c[:len(stored_an)]]
        assert [None if x == 'NULL' else ZZ(K2(elt(x)).trace()) for x in hrow['an']] == out['traces'][:len(stored_an)], r['label']
        if any(hrow[col] != h[col] for col in hrow if col != 'id'):
            hnf_rows.append(hrow)
            stats['%s dim%s smf_hecke_nf change' % (t, '1' if d == 1 else '>1')] += 1
    print(dict(stats))

    # newspace sums
    by_space = defaultdict(list)
    for r in rows:
        by_space[(r['aut_rep_type'], '.'.join(r['label'].split('.')[:-1]))].append(r)
    sp_rows = []
    for t in 'QY':
        pre = SPACE_PREFIX[t]
        labels = [sp for (tt, sp) in by_space if tt == t]
        dimcols = ['eis_Q_dim', 'new_eis_Q_dim'] if t == 'Q' else ['cusp_Y_dim', 'new_cusp_Y_dim']
        spaces = {s['label']: s for s in db.smf_newspaces.search({'label': {'$in': labels}}, ['label'] + dimcols + [pre + h for h in HT])}
        missing = sorted(set(labels) - set(spaces))
        if missing:
            print("%d %s spaces are not in smf_newspaces, skipped: %s" % (len(missing), t, missing))
        labels = [sp for sp in labels if sp in spaces]
        # spaces with no rows of this type are not touched (the level-2 spaces of
        # family P have eis_Q_dim > 0 but no Q orbits and empty eis_Q columns)
        others = db.smf_newspaces.count({dimcols[0]: {'$gt': 0}, 'label': {'$nin': labels}})
        print("%d spaces with %s > 0 and no %s rows, not touched" % (others, dimcols[0], t))
        for sp in sorted(labels):
            s, rs = spaces[sp], by_space[(t, sp)]
            D = sum(r['dim'] for r in rs)
            if all(s.get(pre + h) is None for h in HT):
                stats['%s spaces without stored sums' % t] += 1
                continue
            # (the cusp_Y_dim of the level-2 spaces is not the dimension of the Yoshida lifts,
            # but none of those spaces stores cusp_Y sums)
            assert all(s.get(c) == D for c in dimcols), "%s: %s is not the sum of the orbit dimensions" % (sp, dimcols)
            k, j = rs[0]['weight']
            out, changed = {'label': sp}, False
            for h in HT:
                col = pre + h
                if s.get(col) is None:
                    out[col] = None
                    continue
                stored = [as_int(x) for x in s[col]]
                total = [sum(new_vals[r['label']]['trace_' + h][i] for r in rs) for i in range(len(stored))]
                ps = prime_range(NMAX)[:len(stored)]
                if t == 'Q':
                    old = total if h == 'lambda_p' else [x - (D - 1) * KE_constant(k, j, p, h) for x, p in zip(total, ps)]
                else:
                    old = []
                    for p in ps:
                        Ta = sum(QQ(cache[lab].a(p).trace() if cache[lab].dim > 1 else cache[lab].a(p))
                                 for r in rs for lab in ['.'.join(r['related_objects'][0].split('/')[-4:])])
                        Tb = sum(QQ(cache[lab].a(p).trace() if cache[lab].dim > 1 else cache[lab].a(p))
                                 for r in rs for lab in ['.'.join(r['related_objects'][1].split('/')[-4:])])
                        v = YO(k, j, p, Ta, Tb)[h]
                        old.append(v * (k2_scale(j, p, h) if k == 2 else 1))
                assert old == stored, "%s: stored %s is not what the old computation gives" % (sp, col)
                out[col] = total
                changed = changed or total != stored
            if changed:
                sp_rows.append(out)
                stats['%s newspaces change' % t] += 1
    print(dict(stats))
    print("%d newforms, %d newspaces, %d smf_hecke_nf rows change" % (len(nf_rows), len(sp_rows), len(hnf_rows)))
    for x in nf_rows[:2]:
        print("   ", x['id'], "trace_lambda_p[:4] =", x['trace_lambda_p'][:4], "traces[:9] =", x['traces'][:9])

    base = datafile or "smf_qy_fix"
    todo = []
    if nf_rows:
        fname = base + "_newforms.txt"
        write_update_file(fname, db.smf_newforms, 'id', nf_rows, ['trace_' + h for h in HT] + ['trace_display', 'traces'])
        todo.append((db.smf_newforms, fname, 'id'))
    if sp_rows:
        fname = base + "_newspaces.txt"
        cols = [SPACE_PREFIX[t] + h for t in 'QY' for h in HT]
        for x in sp_rows:
            for col in cols:
                x.setdefault(col, None)   # \N: update_from_file keeps the current value
        write_update_file(fname, db.smf_newspaces, 'label', sp_rows, cols)
        todo.append((db.smf_newspaces, fname, 'label'))
    if hnf_rows:
        fname = base + "_hecke_nf.txt"
        write_update_file(fname, db.smf_hecke_nf, 'id', hnf_rows, SQUARES + ['an'])
        todo.append((db.smf_hecke_nf, fname, 'id'))
    if not todo:
        print("nothing to correct")
        return
    if not apply:
        print("dry run; pass --apply to write")
        return
    # Each update_from_file builds a replacement table; undo with <table>.reload_revert().
    for table, fname, key in todo:
        table.update_from_file(fname, label_col=key)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update files (default: ./smf_qy_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
