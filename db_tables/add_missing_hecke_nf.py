"""
Add the smf_hecke_nf rows of the orbits of dimension <= 20 that have none.

The populate script (common_populate.write_data) skipped every entry without
hecke_ring_numerators.  mf_hecke_nf leaves the numerators NULL when the Hecke ring
basis is the power basis, so no row was written for

  * the Saito-Kurokawa lifts in family K of the classical newforms with a power
    basis: every lift of dimension 1 and some of dimension 2;
  * the orbits of levels 1 and 2 made without a basis: the Siegel Eisenstein
    series, the Klingen Eisenstein series and Saito-Kurokawa lifts of dimension 1,
    and the form of weight (20, 0).

The rows of the lifts in family K are made by SK_all_evs (Hecke_Eigenvalues/
paramodular/Hecke_Eigenvalues_Saito_Kurokawa.py) from the mf_hecke_nf row of the
classical form, as the other lifts.  The other orbits have dimension 1: their
eigenvalues are the traces stored in smf_newforms.  Every new row gives its
basis explicitly, as the stored rows of family K do (for the power basis, the
identity), so it starts with 1.

Before writing anything it checks, for every new row, that the traces of lambda_p,
lambda_p_square, lambda_p_square_{0,1,2} and an are the traces stored in
smf_newforms wherever both are known, that it has the columns of the stored rows,
and that no row of its label or orbit code exists.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.add_missing_hecke_nf [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.add_missing_hecke_nf [--datafile PREFIX] --apply
Undo: delete the rows of the labels in PREFIX_hecke_nf.json.
"""
import argparse
import json
from collections import Counter

from sage.all import QQ, PolynomialRing, NumberField, prime_range
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable
from smf_lmfdb.db_tables.add_level2_lifts import nulls
from smf_lmfdb.db_tables.common_populate import MAX_RING_DEGREE
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import SK_all_evs

MAXP = 97
HT = {'lambda_p': 'trace_lambda_p', 'lambda_p_square': 'trace_lambda_p_square',
      'lambda_p_square_0': 'trace_lambda_p_square_0', 'lambda_p_square_1': 'trace_lambda_p_square_1',
      'lambda_p_square_2': 'trace_lambda_p_square_2'}
NF_COLS = ['label', 'hecke_orbit_code', 'level', 'weight', 'degree', 'family', 'char_orbit_index', 'dim',
           'aut_rep_type', 'related_objects', 'traces'] + list(HT.values())


def identity_basis(e):
    d = len(e['field_poly']) - 1
    ident = [[int(i == j) for j in range(d)] for i in range(d)]
    e['hecke_ring_power_basis'] = False
    e['hecke_ring_numerators'] = e['hecke_ring_inverse_numerators'] = ident
    e['hecke_ring_denominators'] = e['hecke_ring_inverse_denominators'] = [1] * d
    return e


def classical_label(f):
    src = [x for x in (f.get('related_objects') or []) if x.startswith('ModularForm/GL2/Q/holomorphic/')]
    assert len(src) == 1, "%s: related_objects %s" % (f['label'], f.get('related_objects'))
    return '.'.join(src[0].strip('/').split('/')[-4:])


def from_sk(f):
    """The row of a Saito-Kurokawa lift in family K, by SK_all_evs from its classical form."""
    k = f['weight'][0]
    evs = SK_all_evs(k, {'label': classical_label(f)})
    assert len(evs) == 1, "%s: %d rows from SK_all_evs" % (f['label'], len(evs))
    e = evs[0]
    if e.get('hecke_ring_power_basis'):
        e = identity_basis(e)
    return e


def from_traces(f):
    """The row of an orbit of dimension 1, from its traces ('NULL' where a trace is not stored)."""
    assert f['dim'] == 1, f['label']
    ps = prime_range(MAXP + 1)
    e = {'field_poly': [0, 1], 'hecke_ring_rank': 1, 'hecke_ring_cyclotomic_generator': 0}
    e = identity_basis(e)
    for ht, tc in HT.items():
        t = f.get(tc)
        if t is None:
            continue
        e[ht] = ['NULL' if x is None else [x] for x in t[:len(ps)]]
    e['maxp'] = ps[len(e['lambda_p']) - 1]
    e['maxp_square'] = ps[len(e['lambda_p_square']) - 1] if 'lambda_p_square' in e else None
    if f.get('traces') is not None:      # the form of weight (20, 0) has only trace_lambda_p
        e['an'] = ['NULL' if x is None else [x] for x in f['traces'][:MAXP + 1]]
    return e


def traces_of(e):
    """{column: list of traces} of a row (None where the row has 'NULL')."""
    K = NumberField(PolynomialRing(QQ, 'y')(e['field_poly']), 'a') if len(e['field_poly']) > 2 else QQ
    a = K.gen() if K is not QQ else None
    if K is QQ:
        B = [QQ(n[0]) / d for n, d in zip(e['hecke_ring_numerators'], e['hecke_ring_denominators'])]
    else:
        B = [sum(QQ(c) * a**i for i, c in enumerate(n)) / d for n, d in zip(e['hecke_ring_numerators'], e['hecke_ring_denominators'])]

    def tr(c):
        if c in (None, 'NULL'):
            return None
        x = sum(QQ(ci) * b for ci, b in zip(c, B))
        return QQ(x) if K is QQ else x.trace()
    return {col: [tr(c) for c in e[col]] for col in list(HT) + ['an'] if e.get(col) is not None}


def check(e, f):
    """The disagreements of a new row with the traces of its orbit."""
    bad = []
    t = traces_of(e)
    for ht, tc in list(HT.items()) + [('an', 'traces')]:
        if ht not in t or f.get(tc) is None:
            continue
        for i, (x, y) in enumerate(zip(t[ht], f[tc])):
            if x is not None and y is not None and x != y:
                bad.append((ht, i, x, y))
    first = (e['hecke_ring_numerators'][0], e['hecke_ring_denominators'][0])
    if first != ([1] + [0] * (len(e['field_poly']) - 2), 1):
        bad.append(('basis does not start with 1', first))
    return bad


def main(apply=False, datafile=None):
    base = datafile or "smf_missing_hecke_nf"
    cols = [c for c in db.smf_hecke_nf.search_cols if c != 'id']
    have = set(db.smf_hecke_nf.distinct('hecke_orbit_code'))
    have_labels = set(db.smf_hecke_nf.distinct('label'))
    forms = [f for f in db.smf_newforms.search({'dim': {'$lte': MAX_RING_DEGREE}}, NF_COLS) if f['hecke_orbit_code'] not in have]
    ref = list(db.smf_hecke_nf.search({'family': 'K', 'level': {'$gt': 2}}, cols, limit=300))
    always = [c for c in cols if all(r.get(c) is not None for r in ref)]
    print("%d orbits of dimension <= %d without an smf_hecke_nf row" % (len(forms), MAX_RING_DEGREE))
    rows, stats, problems = [], Counter(), []
    for f in forms:
        if (f['family'], f['aut_rep_type']) == ('K', 'P') and f['level'] > 1:
            e, how = from_sk(f), 'SK_all_evs'
        else:
            e, how = from_traces(f), 'traces'
        for c in ['level', 'weight', 'degree', 'family', 'char_orbit_index', 'label', 'hecke_orbit_code']:
            e[c] = f[c]
        assert f['label'] not in have_labels, "%s: a row of this label exists" % f['label']
        bad = check(e, f)
        row = {c: e.get(c) if c == 'an' else nulls(e.get(c)) for c in cols}
        # a column is required unless the newform has no traces to make it from
        absent = {c for c, tc in list(HT.items()) + [('an', 'traces'), ('maxp_square', 'trace_lambda_p_square')] if f.get(tc) is None}
        missing = [c for c in always if row.get(c) is None and c not in absent]
        if missing:
            bad.append(('missing columns', missing))
        if bad:
            problems.append((f['label'], bad[:3]))
            continue
        stats[(f['family'], f['aut_rep_type'], f['dim'], how)] += 1
        rows.append(row)
    for k_, v in sorted(stats.items()):
        print("   %4d  family %s, type %s, dimension %d, from %s" % (v, *k_))
    print("%d rows to add, %d orbits with problems" % (len(rows), len(problems)))
    for p in problems[:10]:
        print("   PROBLEM", p)
    with open(base + "_hecke_nf.json", 'w') as F:
        json.dump(jsonable(rows), F)
    print("wrote %s_hecke_nf.json" % base)
    if problems:
        print("not writing: fix the problems first")
        return
    if not apply:
        print("dry run; pass --apply to write")
        return
    db.smf_hecke_nf.insert_many(jsonable(rows))    # psycopg cannot adapt Sage numbers
    db.smf_hecke_nf.stats.refresh_stats(total=False)
    print("inserted %d rows" % len(rows))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--datafile", help="prefix for the output file (default: ./smf_missing_hecke_nf)")
    parser.add_argument("--apply", action="store_true", help="write to the database")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
