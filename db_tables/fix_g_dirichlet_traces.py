"""
Correct the `traces` column (Dirichlet coefficients of the spinor L-function)
of the paramodular newforms of type G of dimension > 1.

The stored values were computed from the orbit's traces of T(p) and T(p,2) with
the Euler factor and a_{mn} = a_m a_n, as if the orbit were a single form; for
an orbit of dimension d > 1 that is right only at n = p (and wrong at n = 1,
where the trace is d).  The trace over the orbit of a_n is the trace of the
per-embedding coefficient, which smf_hecke_nf stores in `an` (in the Hecke
ring basis) for the orbits that have eigenvalue coordinates.  This script sets

  traces[n-1] = Tr(an[n-1])                     for orbits with an smf_hecke_nf row,
  traces[n-1] = d at n = 1, unchanged at n = p, NULL otherwise   for the others,

keeping NULL wherever n has a prime factor dividing the level.

Before writing anything it checks, for every orbit, that the stored `traces`
are what the old computation gives from the stored trace columns, and, for the
orbits with an smf_hecke_nf row, that the per-embedding an are consistent
with the eigenvalues: a_p = lambda_p, a_{p^2} from the Euler factor at the
good p <= 7, and Tr(a_p) = trace_lambda_p.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_g_dirichlet_traces [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_g_dirichlet_traces [--datafile PREFIX] --apply
"""
import argparse
from collections import Counter

from sage.all import NumberField, PolynomialRing, QQ, ZZ, is_prime, prime_range
from lmfdb import db
from smf_lmfdb.db_tables.fix_degree2_traces import dirichlet_traces, as_int, write_update_file

NF_COLS = ['id', 'label', 'level', 'weight', 'dim', 'traces', 'trace_lambda_p', 'trace_lambda_p_square']
HNF_COLS = ['label', 'an', 'maxp', 'maxp_square', 'lambda_p', 'lambda_p_square', 'field_poly',
            'hecke_ring_power_basis', 'hecke_ring_numerators', 'hecke_ring_denominators']


def hecke_ring_basis(h):
    K = NumberField(PolynomialRing(QQ, 'x')([int(c) for c in h['field_poly']]), 'nu')
    nu = K.gen()
    if h['hecke_ring_power_basis']:
        return [nu**i for i in range(K.degree())]
    return [sum(int(c) * nu**i for i, c in enumerate(num)) / int(den)
            for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])]


def traces_from_an(f, h):
    """Tr(a_n) from the smf_hecke_nf row h, after checking it against the eigenvalues."""
    B = hecke_ring_basis(h)
    assert len(B) == f['dim'], f['label']

    def elt(c):
        return None if c in ('NULL', None) else sum(QQ(x) * b for x, b in zip(c, B))
    an = [elt(c) for c in h['an']]
    lam = dict(zip(prime_range(h['maxp'] + 1), [elt(c) for c in h['lambda_p']]))
    lam2 = dict(zip(prime_range(h['maxp_square'] + 1), [elt(c) for c in h['lambda_p_square']]))
    k, j = f['weight']
    mu = j + 2 * k - 3
    N = f['level']
    for p, a in lam.items():
        if p <= len(an):
            assert an[p - 1] == a, "%s: a_%d is not lambda_%d" % (f['label'], p, p)
    for p in [2, 3, 5, 7]:
        if N % p and p * p <= len(an) and lam.get(p) is not None and lam2.get(p) is not None:
            a, b = lam[p], lam2[p]
            assert an[p * p - 1] == a * a - p**(mu - 2) * (b + p**2 + 1), \
                "%s: a_%d is not given by the Euler factor" % (f['label'], p * p)
    tr = [None if x is None else ZZ(x.trace()) for x in an]
    tl = [as_int(x) for x in f['trace_lambda_p']]
    for i, p in enumerate(prime_range(len(an) + 1)):
        if i < len(tl):
            assert tr[p - 1] == tl[i], "%s: Tr(a_%d) is not trace_lambda_p" % (f['label'], p)
    return tr


def main(apply=False, datafile=None):
    forms = list(db.smf_newforms.search({'aut_rep_type': 'G', 'family': 'K', 'dim': {'$gt': 1}}, NF_COLS))
    hnf = {h['label']: h for h in db.smf_hecke_nf.search({'label': {'$in': [f['label'] for f in forms]}}, HNF_COLS)}
    print("%d G newforms of dimension > 1, %d with smf_hecke_nf rows" % (len(forms), len(hnf)))
    rows, stats = [], Counter()
    for f in forms:
        old = [as_int(x) for x in f['traces']]
        ap = [as_int(x) for x in f['trace_lambda_p']]
        ap2 = [as_int(x) for x in f['trace_lambda_p_square']]
        assert dirichlet_traces(f, ap, ap2, len(old)) == old, \
            "%s: stored traces are not what the old computation gives" % f['label']
        h = hnf.get(f['label'])
        if h is not None and h.get('an'):
            tr = traces_from_an(f, h)
            assert len(tr) >= len(old), f['label']
            new = tr[:len(old)]
            stats['from smf_hecke_nf'] += 1
        else:
            new = [f['dim'] if n == 1 else (old[n - 1] if is_prime(n) else None) for n in range(1, len(old) + 1)]
            stats['no coordinates: composite n set to NULL'] += 1
        # NULL stays NULL (n not prime to the level)
        assert all(new[i] is None for i in range(len(old)) if old[i] is None), f['label']
        if new != old:
            row = dict(f)
            row['traces'] = new
            rows.append(row)
    print(dict(stats))
    print("%d newforms change" % len(rows))
    for r in rows[:3]:
        print("   ", r['label'], "traces[:8] =", r['traces'][:8])
    if not rows:
        print("nothing to correct")
        return
    fname = (datafile or "smf_g_traces_fix") + "_newforms.txt"
    write_update_file(fname, db.smf_newforms, 'id', rows, ['traces'])
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newforms.reload_revert().
    db.smf_newforms.update_from_file(fname, label_col='id')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update file (default: ./smf_g_traces_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
