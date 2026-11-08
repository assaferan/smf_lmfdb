"""
Correct the traces of the weight (3,0) paramodular newforms of type G whose
stored trace_lambda_p / trace_lambda_p_square are not the traces of their
Hecke eigenvalues, and relabel the spaces whose orbit order changes.

omf5_data computed the traces of 1147 degree-2 orbits (250 of type G) with the
Hecke-ring basis B_0, B_1 misread as 1, B_1 (see assaferan/omf5_data#2).
The eigenvalue coordinates, and so smf_hecke_nf, are correct; this script
recomputes the traces from smf_hecke_nf and rewrites:

  smf_newforms:  trace_lambda_p, trace_lambda_p_square, trace_display, traces
                 and, where the order of the G orbits by [dim] + trace_lambda_p
                 changes, label, hecke_orbit, hecke_orbit_code
  smf_hecke_nf:  label, hecke_orbit_code of the relabelled orbits
  smf_newspaces: cusp_G_lambda_p, cusp_G_lambda_p_square of the affected spaces

Before writing anything it checks that recomputing every unaffected form
reproduces its stored values, that the current labels follow the order of the
stored traces, and that the newspace sums equal the sums over their forms.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_degree2_traces [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_degree2_traces [--datafile PREFIX] --apply
"""
import argparse
from collections import defaultdict

from sage.all import NumberField, PolynomialRing, QQ, factor, prime_range
from lmfdb import db
from smf_lmfdb.db_tables.common_populate import MAX_P, base_26
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code

NF_COLS = ['id', 'label', 'space_label', 'level', 'weight', 'dim', 'hecke_orbit',
           'hecke_orbit_code', 'char_orbit_index', 'trace_lambda_p',
           'trace_lambda_p_square', 'trace_display', 'traces', 'atkin_lehner_eigenvals']
HNF_COLS = ['id', 'label', 'field_poly', 'hecke_ring_power_basis', 'hecke_ring_numerators',
            'hecke_ring_denominators', 'lambda_p', 'lambda_p_square']
SP_COLS = ['label', 'cusp_G_lambda_p', 'cusp_G_lambda_p_square']


def as_int(x):
    return None if x is None or x == 'NULL' else int(x)


def traces_from_coordinates(h):
    """Traces of the eigenvalues stored in an smf_hecke_nf row, or None if it has none."""
    if not h['lambda_p']:
        return None
    K = NumberField(PolynomialRing(QQ, 'x')([int(c) for c in h['field_poly']]), 'nu')
    nu = K.gen()
    d = K.degree()
    if h['hecke_ring_power_basis']:
        basis = [nu**i for i in range(d)]
    else:
        basis = [sum(int(c) * nu**i for i, c in enumerate(num)) / int(den)
                 for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])]
    trB = [b.trace() for b in basis]

    def tr(coords):
        out = []
        for c in coords:
            if c == 'NULL' or c is None:
                out.append(None)
            else:
                t = sum(QQ(x) * b for x, b in zip(c, trB))
                assert t.is_integer(), (h['label'], t)
                out.append(int(t))
        return out
    return tr(h['lambda_p']), tr(h['lambda_p_square'])


def dirichlet_traces(rec, ap, ap2, length):
    """
    The `traces` column: a_n of the spinor L-series, from the traces ap of T(p)
    and ap2 of T(p,2) at good primes, via the Euler factor
        1 - a t + p^(mu-2) (b + p^2 + 1) t^2 - p^mu a t^3 + p^(2 mu) t^4,  mu = j + 2k - 3,
    and NULL whenever n has a prime factor dividing the level (a_p is unknown there).

    This is what the stored values satisfy.  (Get_All_Dirichlet_Coeffs_Up_To no
    longer reproduces them: it uses `^`, which is XOR in a .py module.)
    """
    k, j = rec['weight']
    mu = j + 2 * k - 3
    N = rec['level']
    primes = list(prime_range(MAX_P + 1))
    def a_prime_power(p, r):
        i = primes.index(p)
        if N % p == 0 or i >= len(ap) or ap[i] is None:
            return None
        if r >= 2 and (i >= len(ap2) or ap2[i] is None):
            return None
        a = ap[i]
        b = ap2[i] if r >= 2 else 0
        Q = [1, -a, p**(mu - 2) * (b + p**2 + 1), -p**mu * a, p**(2 * mu)]
        c = [1]  # coefficients of 1/Q
        for m in range(1, r + 1):
            c.append(-sum(Q[s] * c[m - s] for s in range(1, min(m, 4) + 1)))
        return c[r]
    out = []
    for n in range(1, length + 1):
        v = 1
        for p, r in factor(n):
            x = a_prime_power(p, r) if p <= MAX_P else None
            if x is None:
                v = None
                break
            v *= x
        out.append(v)
    return out


def sort_key(dim, ap):
    return [dim] + [float('-inf') if x is None else x for x in ap]


def pg(x):
    if isinstance(x, (list, tuple)):
        return "{" + ",".join(pg(y) for y in x) + "}"
    return "NULL" if x is None else str(x)


def write_update_file(fname, table, key, rows, cols):
    types = [table.col_type[c] for c in [key] + cols]
    with open(fname, 'w') as F:
        F.write("|".join([key] + cols) + "\n" + "|".join(types) + "\n\n")
        for r in rows:
            F.write("|".join(pg(r[c]) for c in [key] + cols) + "\n")
    print("wrote %s (%d rows)" % (fname, len(rows)))


def main(apply=False, datafile=None):
    query = {'family': 'K', 'weight': [3, 0], 'aut_rep_type': 'G'}
    forms = list(db.smf_newforms.search(query, NF_COLS))
    hnf = {h['label']: h for h in db.smf_hecke_nf.search({'label': {'$in': [f['label'] for f in forms]}}, HNF_COLS)}
    print("%d G newforms of weight (3,0), %d with smf_hecke_nf rows" % (len(forms), len(hnf)))

    # 1. new traces, and a check that the `traces` formula reproduces every stored
    #    row from its stored (possibly wrong) trace_lambda_p / trace_lambda_p_square
    changed = {}
    reproduced = mismatched = 0
    for f in forms:
        old_ap = [as_int(x) for x in f['trace_lambda_p']]
        old_ap2 = [as_int(x) for x in f['trace_lambda_p_square']]
        if dirichlet_traces(f, old_ap, old_ap2, len(f['traces'])) == [as_int(x) for x in f['traces']]:
            reproduced += 1
        else:
            mismatched += 1
        new = traces_from_coordinates(hnf[f['label']]) if f['label'] in hnf else None
        if new is None:
            continue
        ap, ap2 = new
        ap2 = ap2[:len(old_ap2)]
        if (ap, ap2) != (old_ap, old_ap2):
            changed[f['label']] = (ap, ap2)
    print("%d newforms have wrong traces" % len(changed))
    print("stored `traces` reproduced by the formula: %d, not reproduced: %d" % (reproduced, mismatched))
    assert mismatched == 0, "the `traces` formula does not reproduce the stored values"

    # 2. new orbit order in each affected space
    by_space = defaultdict(list)
    for f in forms:
        by_space[f['space_label']].append(f)
    nf_rows, hnf_rows, relabelled = [], [], {}
    for sp, fs in by_space.items():
        if not any(f['label'] in changed for f in fs):
            continue
        old_ap = lambda f: [as_int(x) for x in f['trace_lambda_p']]
        new_ap = lambda f: changed[f['label']][0] if f['label'] in changed else old_ap(f)
        current = sorted(fs, key=lambda f: f['hecke_orbit'])
        assert [f['hecke_orbit'] for f in current] == list(range(1, len(fs) + 1)), sp
        assert current == sorted(fs, key=lambda f: sort_key(f['dim'], old_ap(f))), \
            "%s: current labels do not follow the stored traces" % sp
        # sorted() is stable, so orbits whose keys tie keep their current order, as
        # they keep their file order in parse_omf5
        new_order = sorted(current, key=lambda f: sort_key(f['dim'], new_ap(f)))
        keys = [sort_key(f['dim'], new_ap(f)) for f in new_order]
        for a, b, fa, fb in zip(keys, keys[1:], new_order, new_order[1:]):
            if a == b:
                print("    tie in %s: %s, %s keep their order" % (sp, fa['label'], fb['label']))
        for X, f in enumerate(new_order, start=1):
            row = dict(f)
            if f['label'] in changed:
                ap, ap2 = changed[f['label']]
                row.update(trace_lambda_p=ap, trace_lambda_p_square=ap2, trace_display=ap[:4],
                           traces=dirichlet_traces(f, ap, ap2, len(f['traces'])))
            if X != f['hecke_orbit']:
                k, j = f['weight']
                row.update(hecke_orbit=X, label=sp + '.' + base_26(X),
                           hecke_orbit_code=make_orbit_code(2, 'K', f['level'], k, j, f['char_orbit_index'], X))
                relabelled[f['label']] = row['label']
                hnf_rows.append({'id': hnf[f['label']]['id'], 'label': row['label'],
                                 'hecke_orbit_code': row['hecke_orbit_code']})
            if f['label'] in changed or X != f['hecke_orbit']:
                nf_rows.append(row)
    print("%d spaces affected, %d of them relabelled (%d orbits get a new label)"
          % (len({r['space_label'] for r in nf_rows}), len({l.rsplit('.', 1)[0] for l in relabelled}), len(relabelled)))
    for old, new in list(relabelled.items())[:6]:
        print("   ", old, "->", new)

    # 3. newspace sums
    sp_rows = []
    delta = defaultdict(lambda: None)
    for sp in {r['space_label'] for r in nf_rows}:
        s = db.smf_newspaces.lucky({'label': sp}, SP_COLS)
        fs = by_space[sp]
        for col, fcol, i in [('cusp_G_lambda_p', 'trace_lambda_p', 0), ('cusp_G_lambda_p_square', 'trace_lambda_p_square', 1)]:
            stored = [as_int(x) for x in s[col]]
            old_sum = [None if all(as_int(f[fcol][n]) is None for f in fs) else sum(as_int(f[fcol][n]) or 0 for f in fs)
                       for n in range(len(stored))]
            assert old_sum == stored, "%s: %s is not the sum over its G forms" % (sp, col)
            new_vals = [changed[f['label']][i] if f['label'] in changed else [as_int(x) for x in f[fcol]] for f in fs]
            s[col] = [None if all(v[n] is None for v in new_vals) else sum(v[n] or 0 for v in new_vals)
                      for n in range(len(stored))]
        sp_rows.append(s)
    print("%d newspaces to update" % len(sp_rows))

    base = datafile or "smf_degree2_fix"
    nf_file, hnf_file, sp_file = base + "_newforms.txt", base + "_hecke_nf.txt", base + "_newspaces.txt"
    write_update_file(nf_file, db.smf_newforms, 'id', nf_rows,
                      ['label', 'hecke_orbit', 'hecke_orbit_code', 'trace_lambda_p',
                       'trace_lambda_p_square', 'trace_display', 'traces'])
    write_update_file(hnf_file, db.smf_hecke_nf, 'id', hnf_rows, ['label', 'hecke_orbit_code'])
    write_update_file(sp_file, db.smf_newspaces, 'label', sp_rows, ['cusp_G_lambda_p', 'cusp_G_lambda_p_square'])
    if not apply:
        print("dry run; pass --apply to write")
        return
    # Each update_from_file builds a replacement table; undo with <table>.reload_revert().
    db.smf_newforms.update_from_file(nf_file, label_col='id')
    if hnf_rows:
        db.smf_hecke_nf.update_from_file(hnf_file, label_col='id')
    db.smf_newspaces.update_from_file(sp_file, label_col='label')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update files (default: ./smf_degree2_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
