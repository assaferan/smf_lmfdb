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

Before writing anything it checks that the formula for `traces` reproduces
every stored row, that only forms of degree 2 change, that the current labels
follow the order of the stored traces, that smf_newforms and smf_hecke_nf agree
on labels and orbit codes, and that the sums of the unaffected newspaces equal
the sums over their forms.

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
           'trace_lambda_p_square', 'trace_display', 'traces']
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
    hnf_list = list(db.smf_hecke_nf.search({'label': {'$in': [f['label'] for f in forms]}}, HNF_COLS + ['hecke_orbit_code']))
    hnf = {h['label']: h for h in hnf_list}
    print("%d G newforms of weight (3,0), %d with smf_hecke_nf rows" % (len(forms), len(hnf)))
    # The two tables are joined by label.  After a partial apply (smf_newforms
    # relabelled, smf_hecke_nf not) a label would point at another orbit's
    # coordinates, so refuse to run unless labels and orbit codes agree.
    codes = {f['label']: f['hecke_orbit_code'] for f in forms}
    bad = [h['label'] for h in hnf_list if codes[h['label']] != h['hecke_orbit_code']]
    assert not bad, "smf_newforms and smf_hecke_nf disagree on label/hecke_orbit_code: %s" % bad[:5]

    # 1. new traces from the coordinates; check that the `traces` formula
    #    reproduces every stored row from its stored trace_lambda_p / _square
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
        assert len(ap) == len(old_ap) and len(ap2) >= len(old_ap2), \
            "%s: smf_hecke_nf has %d/%d eigenvalues, smf_newforms %d/%d traces" % (
                f['label'], len(ap), len(ap2), len(old_ap), len(old_ap2))
        ap2 = ap2[:len(old_ap2)]
        if (ap, ap2) != (old_ap, old_ap2):
            changed[f['label']] = (ap, ap2)
    print("%d newforms have wrong traces" % len(changed))
    print("stored `traces` reproduced by the formula: %d, not reproduced: %d" % (reproduced, mismatched))
    assert mismatched == 0, "the `traces` formula does not reproduce the stored values"
    # The known bug affects only degree 2; anything else changing means the
    # data or the basis convention has changed, and the script should not guess.
    dims = {f['label']: f['dim'] for f in forms}
    other = [l for l in changed if dims[l] != 2]
    assert not other, "traces of forms of degree != 2 would change: %s" % other[:5]

    # 2. new orbit order in each affected space
    by_space = defaultdict(list)
    for f in forms:
        by_space[f['space_label']].append(f)
    nf_rows, hnf_rows, relabelled = [], [], {}
    old_ap = lambda f: [as_int(x) for x in f['trace_lambda_p']]
    new_ap = lambda f: changed[f['label']][0] if f['label'] in changed else old_ap(f)
    for sp, fs in by_space.items():
        if not any(f['label'] in changed for f in fs):
            continue
        current = sorted(fs, key=lambda f: f['hecke_orbit'])
        assert [f['hecke_orbit'] for f in current] == list(range(1, len(fs) + 1)), sp
        # sorted() is stable, so orbits whose keys tie keep their current order, as
        # they keep their file order in parse_omf5
        assert current == sorted(current, key=lambda f: sort_key(f['dim'], old_ap(f))), \
            "%s: current labels do not follow the stored traces" % sp
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
                if f['label'] in hnf:  # forms of large degree have no smf_hecke_nf row
                    hnf_rows.append({'id': hnf[f['label']]['id'], 'label': row['label'],
                                     'hecke_orbit_code': row['hecke_orbit_code']})
            if f['label'] in changed or X != f['hecke_orbit']:
                nf_rows.append(row)
    print("%d spaces affected, %d of them relabelled (%d orbits get a new label)"
          % (len({r['space_label'] for r in nf_rows}), len({l.rsplit('.', 1)[0] for l in relabelled}), len(relabelled)))
    for old, new in list(relabelled.items())[:6]:
        print("   ", old, "->", new)

    # 3. newspace sums: compare every space with the (corrected) traces of its G
    #    forms, not only the spaces of changed forms, so that a rerun also repairs
    #    a space left behind by a partial apply
    def sums(fs, i, fcol, length):
        vals = [changed[f['label']][i] if f['label'] in changed else [as_int(x) for x in f[fcol]] for f in fs]
        return [None if all(v[n] is None for v in vals) else sum(v[n] or 0 for v in vals) for n in range(length)]
    sp_rows = []
    for s in db.smf_newspaces.search({'label': {'$in': list(by_space)}}, SP_COLS):
        fs = by_space[s['label']]
        new = {}
        for col, fcol, i in [('cusp_G_lambda_p', 'trace_lambda_p', 0), ('cusp_G_lambda_p_square', 'trace_lambda_p_square', 1)]:
            stored = [as_int(x) for x in s[col]]
            if not any(f['label'] in changed for f in fs):
                assert sums(fs, i, fcol, len(stored)) == stored, "%s: %s is not the sum over its G forms" % (s['label'], col)
            new[col] = sums(fs, i, fcol, len(stored))
            if new[col] != stored:
                s[col] = new[col]
                if s not in sp_rows:
                    sp_rows.append(s)
    print("%d newspaces to update" % len(sp_rows))

    if not (nf_rows or hnf_rows or sp_rows):
        print("nothing to correct")
        return
    base = datafile or "smf_degree2_fix"
    nf_file, hnf_file, sp_file = base + "_newforms.txt", base + "_hecke_nf.txt", base + "_newspaces.txt"
    todo = []
    if nf_rows:
        write_update_file(nf_file, db.smf_newforms, 'id', nf_rows,
                          ['label', 'hecke_orbit', 'hecke_orbit_code', 'trace_lambda_p',
                           'trace_lambda_p_square', 'trace_display', 'traces'])
        todo.append((db.smf_newforms, nf_file, 'id'))
    if hnf_rows:
        write_update_file(hnf_file, db.smf_hecke_nf, 'id', hnf_rows, ['label', 'hecke_orbit_code'])
        todo.append((db.smf_hecke_nf, hnf_file, 'id'))
    if sp_rows:
        write_update_file(sp_file, db.smf_newspaces, 'label', sp_rows, ['cusp_G_lambda_p', 'cusp_G_lambda_p_square'])
        todo.append((db.smf_newspaces, sp_file, 'label'))
    if not apply:
        print("dry run; pass --apply to write")
        return
    # Each update_from_file builds a replacement table; undo with <table>.reload_revert().
    # The updates are not atomic: if one fails, rerun the script (the checks above
    # refuse an inconsistent label state; revert smf_newforms first in that case).
    for table, fname, key in todo:
        table.update_from_file(fname, label_col=key)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update files (default: ./smf_degree2_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
