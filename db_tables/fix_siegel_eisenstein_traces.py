"""
Correct `traces` of the Siegel Eisenstein series (aut_rep_type 'F') in smf_newforms.

The spin L-function of E_k (level 1, weight (k, 0)) is
    zeta(s) zeta(s-k+2) zeta(s-k+1) zeta(s-2k+3),
so a_n = sum over n = abcd of b^(k-2) c^(k-1) d^(2k-3), and every a_n is positive.
The stored traces of the 27 level-1 rows (E_4, ..., E_20 in K, P and S) and of the 9
level-2 rows (Gamma2_lifts.siegel_eisenstein copies the level-1 traces) were computed
by the Get_All_Dirichlet_Coeffs_Up_To of before 1295f51, which reads lambda_p_square
as T(p,2): it applied T(p^2) = l_p^2 - p^(mu-1) - (T(p,2) + p^2 + 1) p^(mu-2) to the
stored T(p^2), so a_4 of E_4 is -10407 instead of 1565.  Their lambda_p and
lambda_p_square are right; they have no smf_hecke_nf rows, and their spaces store no
`traces` sums.

Rows whose traces are already right are skipped.  Before writing anything it checks, for
each other row:
  - dimension 1, j = 0;
  - trace_lambda_p and trace_lambda_p_square are those of E_k;
  - the stored traces are exactly what the old computation gives (so that nothing
    else is being overwritten);
  - the new traces agree with the current Get_All_Dirichlet_Coeffs_Up_To on the
    stored lambda_p, lambda_p_square, wherever that is defined.

The new values are written to PREFIX.json ({label: [old, new]}).

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_siegel_eisenstein_traces [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_siegel_eisenstein_traces [--datafile PREFIX] --apply
"""
import argparse
import json

from sage.all import ZZ, PowerSeriesRing, factor, prime_range
from lmfdb import db

COLS = ['id', 'label', 'level', 'weight', 'dim', 'traces', 'trace_lambda_p', 'trace_lambda_p_square']


def euler_inverse(k, p, r):
    """The coefficient of t^r in 1 / ((1-t)(1-p^(k-2)t)(1-p^(k-1)t)(1-p^(2k-3)t))."""
    t = PowerSeriesRing(ZZ, 't', default_prec=r + 1).gen()
    return (1 / ((1 - t) * (1 - p**(k-2) * t) * (1 - p**(k-1) * t) * (1 - p**(2*k-3) * t)))[r]


def spin(k, n):
    a = ZZ(1)
    for p, r in factor(n):
        a *= euler_inverse(k, p, r)
    return a


def from_euler(k, lam, lam2, n, old):
    """a_n from the Euler factor 1 - l t + (l^2 - T - p^(mu-1)) t^2 - p^mu l t^3 + p^(2mu) t^4, with T
    = T(p^2) (current code) or T computed from lambda_p_square read as T(p,2) (old code); None
    when an eigenvalue is missing."""
    mu = 2 * k - 3
    a = ZZ(1)
    for p, r in factor(n):
        if p not in lam or (r > 1 and p not in lam2):
            return None
        l = lam[p]
        if r == 1:
            a *= l
            continue
        T = lam2[p]
        if old:
            T = l**2 - p**(mu-1) - (T + p**2 + 1) * p**(mu-2)
        t = PowerSeriesRing(ZZ, 't', default_prec=r + 1).gen()
        a *= (1 / (1 - l * t + (l**2 - T - p**(mu-1)) * t**2 - p**mu * l * t**3 + p**(2*mu) * t**4))[r]
    return a


def main(apply=False, datafile=None):
    base = datafile or "smf_siegel_eisenstein_traces"
    rows = list(db.smf_newforms.search({'aut_rep_type': 'F'}, COLS))
    out = {}
    for r in rows:
        k, j = r['weight']
        assert j == 0 and r['dim'] == 1, r['label']
        ps = prime_range(1000)
        lam = dict(zip(ps, map(ZZ, r['trace_lambda_p'])))
        lam2 = {p: ZZ(x) for p, x in zip(ps, r['trace_lambda_p_square']) if x is not None}
        assert all(x == spin(k, p) for p, x in lam.items()), "%s: lambda_p is not that of E_%d" % (r['label'], k)
        # T(p^2) of E_k: a_{p^2} - p^(mu-1)
        assert all(x == spin(k, p * p) - p**(2*k-4) for p, x in lam2.items()), \
            "%s: lambda_p_square is not that of E_%d" % (r['label'], k)
        old = r['traces']
        n = len(old)
        if [int(x) for x in old] == [int(spin(k, m)) for m in range(1, n + 1)]:
            print("%s: already correct" % r['label'])
            continue
        expect_old = [from_euler(k, lam, lam2, m, True) for m in range(1, n + 1)]
        assert all(a == b for a, b in zip(old, expect_old)), \
            "%s: the stored traces are not the old computation at n = %s" % (
                r['label'], [m for m, (a, b) in enumerate(zip(old, expect_old), 1) if a != b][:5])
        new = [spin(k, m) for m in range(1, n + 1)]
        cur = [from_euler(k, lam, lam2, m, False) for m in range(1, n + 1)]
        assert all(c is None or c == x for c, x in zip(cur, new)), r['label']
        changed = [m for m in range(1, n + 1) if old[m-1] != new[m-1]]
        out[r['label']] = {'id': r['id'], 'old': [None if x is None else int(x) for x in old],
                           'new': [int(x) for x in new]}
        print("%s: %d traces, %d change (first n: %s), %d were NULL" % (
            r['label'], n, len(changed), changed[:6], sum(1 for x in old if x is None)))
    print("%d rows" % len(out))
    with open(base + ".json", 'w') as F:
        json.dump(out, F)
    print("wrote %s.json" % base)
    if not apply:
        print("dry run; pass --apply to write")
        return
    for label, v in out.items():
        db.smf_newforms.update({'label': label}, {'traces': v['new']}, restat=False)
    db.smf_newforms.stats.refresh_stats(total=False)
    print("updated")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--datafile", help="prefix for the output file (default: ./smf_siegel_eisenstein_traces)")
    parser.add_argument("--apply", action="store_true", help="write to the database")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
