"""
Correct the traces of the Saito-Kurokawa lifts (aut_rep_type 'P') in
smf_newforms and their sums in smf_newspaces.

A lift F of a classical newform orbit f of weight w = 2k-2 has, at each
embedding, Hecke eigenvalues given by the SK_* formulas of
Hecke_Eigenvalues_Saito_Kurokawa in terms of a_p(f), and spinor L-function
L(F, s) = zeta(s-k+1) zeta(s-k+2) L(f, s).  The stored values were computed

  * by applying the SK_* formulas once to Tr(a_p) instead of summing them over
    the embeddings (wrong for every orbit of dimension > 1: the constant terms
    are counted once and the square formulas are not linear),
  * with `p^2-1` in SK_lambda_p_square_0, which in a .py module is p XOR 1,
  * and `traces` from the Euler factor of the old Get_All_Dirichlet_Coeffs_Up_To,
    which reads trace_lambda_p_square as a trace of T(p,2) (it is one of T(p^2))
    and multiplies traces of the orbit together.

This script recomputes, from the classical data (mf_newforms.traces,
mf_hecke_nf.ap, mf_hecke_charpolys),

  smf_newforms:  trace_lambda_p, trace_lambda_p_square[_0,_1,_2], trace_display,
                 traces
  smf_newspaces: cusp_P_lambda_p, cusp_P_lambda_p_square[_0,_1,_2] (sums over the
                 P orbits of the space, for the columns that are stored)
  smf_hecke_nf:  lambda_p, lambda_p_square[_0,_1,_2], an (the eigenvalues and
                 Dirichlet coefficients, in the Hecke ring basis of the classical
                 form; the stored ones were converted with the wrong inverse basis
                 of the old get_nf_basis(two=True), and have the right traces but
                 are the wrong elements)

Conventions, matching the G rows:  `traces` and `an` are NULL at every n with a
prime factor dividing the level, and an entry of a square column is NULL when
Tr(a_p^2) is not available from the classical data (orbits without mf_hecke_nf
coordinates, p > 97).  The values of the lambda columns at p | N are the SK
formulas, as before, now summed over the embeddings.

Not touched (reported): rows lifting a form f of root number +1, which has no
lift (the old code took Fricke sign -1 for every weight, which is right only for
k odd, and no condition at level 1: the level-1 rows of odd k, and
2.K.31.4.0.a.b), and the level-2 rows of family P.

Before writing anything it checks that every row's classical form has the
expected level, weight, dimension and trivial character, that the sources for Tr(a_p^2) agree where they overlap, that the stored
dimension-1 rows equal the formulas (with p XOR 1 for _square_0), that the
stored trace_lambda_p of the other rows equals the formula applied to Tr(a_p),
and that the stored new_cusp_P_dim of each affected space is the sum of its
orbits' dimensions.  Labels do not change: the P orbits of a space are ordered by their
classical labels, not by their traces.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_sk_traces [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_sk_traces [--datafile PREFIX] --apply
"""
import argparse
from collections import defaultdict, Counter

import json

from sage.all import Matrix, NumberField, PolynomialRing, QQ, ZZ, divisors, prime_range, vector
from lmfdb import db
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import (
    SK_lambda_p, SK_lambda_p_square, SK_lambda_p_square_0, SK_lambda_p_square_1, SK_lambda_p_square_2)

PRIMES = list(prime_range(200))
SK = {'lambda_p': SK_lambda_p,
      'lambda_p_square': SK_lambda_p_square,
      'lambda_p_square_0': SK_lambda_p_square_0,
      'lambda_p_square_1': SK_lambda_p_square_1,
      'lambda_p_square_2': SK_lambda_p_square_2}
HT = list(SK)
NF_COLS = ['id', 'label', 'family', 'space_label', 'level', 'weight', 'dim', 'hecke_orbit',
           'related_objects', 'trace_display', 'traces'] + ['trace_' + h for h in HT]
SP_COLS = ['label', 'new_cusp_P_dim'] + ['cusp_P_' + h for h in HT]


def square_0_as_stored(k, p, a):
    # SK_lambda_p_square_0 as it ran when the data was made: p^2-1 is p ^ (2-1) = p XOR 1
    return a**2 + p**(k - 3) * ((p ^ 1)) * (a + p**(k - 1))


def coefficients(F, k, p):
    """c0, c1, c2 with F(k, p, a) = c0 + c1 a + c2 a^2 (all SK_* are of degree <= 2 in a)."""
    c0, s, m = F(k, p, 0), F(k, p, 1), F(k, p, -1)
    c1, c2 = (s - m) // 2, (s + m - 2 * c0) // 2
    assert all(F(k, p, a) == c0 + c1 * a + c2 * a**2 for a in [2, -3, 7]), (F, k, p)
    return c0, c1, c2


def classical_label(r):
    return '.'.join(r['related_objects'][0].split('holomorphic/')[1].split('/'))


def nf_basis(h):
    K = NumberField(PolynomialRing(QQ, 'x')([int(c) for c in h['field_poly']]), 'nu')
    nu = K.gen()
    d = K.degree()
    if h['hecke_ring_power_basis']:
        return [nu**i for i in range(d)]
    return [sum(int(c) * nu**i for i, c in enumerate(num)) / int(den)
            for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])]


def trace_of_squares(f, hnf, charpolys):
    """{p: Tr(a_p^2)} for the primes p < 200 for which the classical data determines it."""
    N, w, d, tr = f['level'], int(f['weight']), f['dim'], f['traces']
    t2 = {}
    sources = Counter()

    def put(p, v, src):
        if p in t2:
            assert t2[p] == v, "%s: Tr(a_%d^2) is %s from one source, %s from %s" % (f['label'], p, t2[p], v, src)
        else:
            t2[p] = v
            sources[src] += 1
    for p in PRIMES:
        if p * p <= len(tr):  # a_{p^2} = a_p^2 - p^(w-1) for p not dividing N (trivial character), a_p^2 otherwise
            put(p, tr[p * p - 1] + (d * p**(w - 1) if N % p else 0), 'traces')
    h = hnf.get(f['hecke_orbit_code'])
    if h is not None:
        B = nf_basis(h)
        for p, c in zip(prime_range(h['maxp'] + 1), h['ap']):
            if p < 200:
                a = sum(QQ(x) * b for x, b in zip(c, B))
                put(p, ZZ((a * a).trace()), 'mf_hecke_nf')
    for p, fac in charpolys.get(f['hecke_orbit_code'], {}).items():
        s = 0
        for coeffs, mult in fac:  # coefficients from low to high degree, monic
            e = len(coeffs) - 1
            e1 = -coeffs[e - 1]
            e2 = coeffs[e - 2] if e >= 2 else 0
            s += mult * (e1**2 - 2 * e2)
        put(p, ZZ(s), 'mf_hecke_charpolys')
    return t2


def expected_row(r, f, t2):
    k, N, d, tr = r['weight'][0], r['level'], r['dim'], f['traces']
    out = {}
    for h, F in SK.items():
        vals = []
        for p in PRIMES:
            c0, c1, c2 = coefficients(F, k, p)
            if c2 and p not in t2:
                vals.append(None)
            else:
                vals.append(d * c0 + c1 * tr[p - 1] + (c2 * t2[p] if c2 else 0))
        out['trace_' + h] = vals
    out['trace_display'] = out['trace_lambda_p'][:4]
    # L(F, s) = zeta(s-k+1) zeta(s-k+2) L(f, s), at n prime to N
    out['traces'] = [None if ZZ(n).gcd(N) != 1 else
                     sum(d1**(k - 1) * sum(d2**(k - 2) * tr[n // d1 // d2 - 1] for d2 in divisors(n // d1))
                         for d1 in divisors(n))
                     for n in range(1, len(r['traces']) + 1)]
    return out


def as_int(x):
    return None if x is None or x == 'NULL' else int(x)


def pg(x, typ=None):
    # A NULL field is \N in the COPY file; the literal NULL is only valid inside an array.
    if x is None:
        return "\\N"
    if typ == 'jsonb':
        return json.dumps(x)
    if isinstance(x, (list, tuple)):
        return "{" + ",".join("NULL" if y is None else pg(y) for y in x) + "}"
    return str(x)


def write_update_file(fname, table, key, rows, cols):
    types = [table.col_type[c] for c in [key] + cols]
    with open(fname, 'w') as F:
        F.write("|".join([key] + cols) + "\n" + "|".join(types) + "\n\n")
        for r in rows:
            F.write("|".join(pg(r[c], t) for c, t in zip([key] + cols, types)) + "\n")
    print("wrote %s (%d rows)" % (fname, len(rows)))


def main(apply=False, datafile=None):
    rows = list(db.smf_newforms.search({'aut_rep_type': 'P'}, NF_COLS))
    labels = sorted({classical_label(r) for r in rows})
    cl = {f['label']: f for f in db.mf_newforms.search(
        {'label': {'$in': labels}},
        ['label', 'level', 'weight', 'dim', 'char_orbit_index', 'fricke_eigenval', 'hecke_orbit_code', 'traces'])}
    assert len(cl) == len(labels), "classical forms not found: %s" % sorted(set(labels) - set(cl))[:5]
    # A lift needs f to have root number (-1)^(w/2) * fricke = -1, w = 2k-2.  The old code
    # took fricke = -1 for every k (right only for k odd), and no condition at level 1.
    invalid = [r for r in rows if not (r['family'] == 'P' and r['level'] == 2) and
               (-1)**(r['weight'][0] - 1) * cl[classical_label(r)]['fricke_eigenval'] != -1]
    skipped = invalid + [r for r in rows if r['family'] == 'P' and r['level'] == 2]
    skip_spaces = {r['space_label'] for r in skipped}
    rows = [r for r in rows if r['space_label'] not in skip_spaces]
    print("%d P rows to check; not touched: %d rows in %d spaces" % (len(rows), len(skipped), len(skip_spaces)))
    print("    lifts of a form with root number +1 (no such lift exists): %d: %s"
          % (len(invalid), ", ".join(sorted(r['label'] for r in invalid if r['level'] > 1)
                                    + ["%d at level 1" % sum(r['level'] == 1 for r in invalid)])))
    print("    level 2, family P (not checked): %d" % (len(skipped) - len(invalid)))
    for r in rows:
        f = cl[classical_label(r)]
        k = r['weight'][0]
        assert r['weight'][1] == 0 and (f['level'], int(f['weight']), f['dim'], f['char_orbit_index']) == \
            (r['level'], 2 * k - 2, r['dim'], 1), "%s does not match its classical form %s" % (r['label'], f['label'])
        assert len(f['traces']) >= 200, f['label']
    codes = [f['hecke_orbit_code'] for f in cl.values()]
    hnf = {h['hecke_orbit_code']: h for h in db.mf_hecke_nf.search(
        {'hecke_orbit_code': {'$in': codes}},
        ['hecke_orbit_code', 'ap', 'maxp', 'field_poly', 'hecke_ring_power_basis',
         'hecke_ring_numerators', 'hecke_ring_denominators'])}
    no_hnf = [f['hecke_orbit_code'] for f in cl.values() if f['hecke_orbit_code'] not in hnf and f['dim'] > 1]
    charpolys = defaultdict(dict)
    if no_hnf:
        for x in db.mf_hecke_charpolys.search({'hecke_orbit_code': {'$in': no_hnf}},
                                              ['hecke_orbit_code', 'p', 'charpoly_factorization']):
            charpolys[x['hecke_orbit_code']][x['p']] = x['charpoly_factorization']
    print("%d classical orbits, %d with mf_hecke_nf, charpolys for %d of the other %d of dimension > 1"
          % (len(cl), len(hnf), len(charpolys), len(no_hnf)))

    # new values, and checks that the stored values are what the diagnosed bugs produce
    t2_cache = {}
    nf_rows, unknown = [], Counter()
    for r in rows:
        f = cl[classical_label(r)]
        if f['label'] not in t2_cache:
            t2_cache[f['label']] = trace_of_squares(f, hnf, charpolys)
        t2 = t2_cache[f['label']]
        new = expected_row(r, f, t2)
        k, d, tr = r['weight'][0], r['dim'], f['traces']
        stored = {c: [as_int(x) for x in r[c]] for c in new}
        assert stored['trace_lambda_p'] == [SK_lambda_p(k, p, tr[p - 1]) for p in PRIMES], \
            "%s: stored trace_lambda_p is not SK_lambda_p(Tr a_p)" % r['label']
        if d == 1:
            for h in ['lambda_p_square', 'lambda_p_square_1', 'lambda_p_square_2']:
                assert stored['trace_' + h] == new['trace_' + h], "%s: stored trace_%s is wrong" % (r['label'], h)
            # level 1 was made by SK_Hecke_p_square_0, which has p**2; the others by SK_lambda_p_square_0
            sq0 = new['trace_lambda_p_square_0'] if r['level'] == 1 else \
                [square_0_as_stored(k, p, tr[p - 1]) for p in PRIMES]
            assert stored['trace_lambda_p_square_0'] == sq0, \
                "%s: stored trace_lambda_p_square_0 is not what the old code computes" % r['label']
        for c, v in new.items():
            if any(x is None for x in v) and c != 'traces':
                unknown[c] += 1
        if any(stored[c] != new[c] for c in new):
            row = dict(r)
            row.update(new)
            nf_rows.append(row)
    print("%d smf_newforms rows change; rows with NULL entries for lack of Tr(a_p^2): %s"
          % (len(nf_rows), dict(unknown)))

    # labels: the P orbits of a space follow the order of their classical labels
    by_space = defaultdict(list)
    for r in rows:
        by_space[r['space_label']].append(r)
    for sp, rs in by_space.items():
        by_orbit = [classical_label(r) for r in sorted(rs, key=lambda r: r['hecke_orbit'])]
        assert by_orbit == sorted(by_orbit, key=lambda l: (int(l.split('.')[0]), int(l.split('.')[1]), l.split('.')[2],
                                                             len(l.split('.')[3]), l.split('.')[3])), \
            "%s: P orbits are not in the order of their classical labels" % sp

    # newspace sums
    new_by_label = {x['label']: x for x in nf_rows}
    sp_rows = []
    spaces = {s['label']: s for s in db.smf_newspaces.search({'label': {'$in': list(by_space)}}, SP_COLS)}
    assert len(spaces) == len(by_space), sorted(set(by_space) - set(spaces))[:5]
    bad_dim = [(sp, spaces[sp].get('new_cusp_P_dim'), sum(r['dim'] for r in rs)) for sp, rs in by_space.items()
               if spaces[sp].get('new_cusp_P_dim') != sum(r['dim'] for r in rs)]
    assert not bad_dim, "new_cusp_P_dim is not the sum over the P orbits (space, stored, sum): %s" % sorted(bad_dim)
    for sp, rs in sorted(by_space.items()):
        s = spaces[sp]
        out, changed = {'label': sp}, False
        for h in HT:
            col = 'cusp_P_' + h
            if s.get(col) is None:
                continue
            stored = [as_int(x) for x in s[col]]
            vals = [new_by_label.get(r['label'], r)['trace_' + h] for r in rs]
            vals = [[as_int(x) for x in v] for v in vals]
            total = [None if any(v[i] is None for v in vals) else sum(v[i] for v in vals)
                     for i in range(len(stored))]
            out[col] = total
            changed = changed or (total != stored)
        if changed:
            for h in HT:  # columns not stored stay NULL; update_from_file keeps the current value
                out.setdefault('cusp_P_' + h, s.get('cusp_P_' + h))
            sp_rows.append(out)
    print("%d of %d newspaces change" % (len(sp_rows), len(by_space)))

    hnf_rows = hecke_nf_rows(rows, cl)
    for x in nf_rows[:3]:
        print("   ", x['label'], "trace_lambda_p[:3] =", x['trace_lambda_p'][:3], "traces[:6] =", x['traces'][:6])

    if not (nf_rows or sp_rows or hnf_rows):
        print("nothing to correct")
        return
    base = datafile or "smf_sk_fix"
    nf_file, sp_file, hnf_file = base + "_newforms.txt", base + "_newspaces.txt", base + "_hecke_nf.txt"
    todo = []
    if nf_rows:
        write_update_file(nf_file, db.smf_newforms, 'id', nf_rows,
                          ['trace_' + h for h in HT] + ['trace_display', 'traces'])
        todo.append((db.smf_newforms, nf_file, 'id'))
    if sp_rows:
        write_update_file(sp_file, db.smf_newspaces, 'label', sp_rows, ['cusp_P_' + h for h in HT])
        todo.append((db.smf_newspaces, sp_file, 'label'))
    if hnf_rows:
        write_update_file(hnf_file, db.smf_hecke_nf, 'id', hnf_rows, HT + ['an'])
        todo.append((db.smf_hecke_nf, hnf_file, 'id'))
    if not apply:
        print("dry run; pass --apply to write")
        return
    # Each update_from_file builds a replacement table; undo with <table>.reload_revert().
    for table, fname, key in todo:
        table.update_from_file(fname, label_col=key)


HNF_RING = ['field_poly', 'hecke_ring_power_basis', 'hecke_ring_numerators', 'hecke_ring_denominators']


def hecke_nf_rows(rows, cl):
    """
    Corrected smf_hecke_nf rows of the P orbits in rows: the SK formulas applied
    to the classical eigenvalues at each embedding, in the Hecke ring basis that
    the row stores (that of the classical form).
    """
    by_label = {r['label']: r for r in rows}
    hs = list(db.smf_hecke_nf.search({'label': {'$in': list(by_label)}},
                                     ['id', 'label', 'maxp', 'maxp_square', 'an'] + HNF_RING + HT))
    cls = {m['label']: m for m in db.mf_hecke_nf.search(
        {'label': {'$in': list({classical_label(by_label[h['label']]) for h in hs})}},
        ['label', 'an', 'ap', 'maxp'] + HNF_RING)}
    out, trace_kept = [], Counter()
    for h in hs:
        r = by_label[h['label']]
        m = cls[classical_label(r)]
        assert all(h[c] == m[c] for c in HNF_RING), "%s: Hecke ring differs from that of %s" % (h['label'], m['label'])
        k, N = r['weight'][0], r['level']
        B = nf_basis(m)
        K = B[0].parent()
        Minv = Matrix(QQ, [list(K(b)) for b in B]).inverse()

        def coords(x):
            c = vector(QQ, list(K(x))) * Minv
            assert all(t.is_integer() for t in c), (h['label'], x)
            return [int(t) for t in c]

        def elt(c):
            return sum(QQ(t) * b for t, b in zip(c, B))
        ap = dict(zip(prime_range(m['maxp'] + 1), [elt(c) for c in m['ap']]))
        an = [elt(c) for c in m['an']]  # an[n-1] = a_n
        new = {'id': h['id']}
        for ht, F in SK.items():
            ps = prime_range((h['maxp'] if ht == 'lambda_p' else h['maxp_square']) + 1)
            assert len(ps) == len(h[ht]), (h['label'], ht)
            new[ht] = [coords(F(k, p, ap[p])) for p in ps]
        assert len(h['an']) <= len(an), h['label']
        new['an'] = ['NULL' if ZZ(n).gcd(N) != 1 else
                     coords(sum(d1**(k - 1) * sum(d2**(k - 2) * an[n // d1 // d2 - 1] for d2 in divisors(n // d1))
                                for d1 in divisors(n)))
                     for n in range(1, len(h['an']) + 1)]
        old_tr = [elt(c).trace() for c in h['lambda_p']]
        trace_kept[old_tr == [elt(c).trace() for c in new['lambda_p']]] += 1
        if any(new[c] != h[c] for c in HT + ['an']):
            out.append(new)
    print("%d smf_hecke_nf rows of P orbits, %d change; stored lambda_p has the right traces in %d"
          % (len(hs), len(out), trace_kept[True]))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update files (default: ./smf_sk_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
