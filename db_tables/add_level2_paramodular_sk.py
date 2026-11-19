"""
Add the Saito-Kurokawa lifts of paramodular level 2 (family K, weight (k,0)).

No form of family K and level 2 is stored, although the spaces 2.K.2.k.0.a record
the dimension of their SK lifts (new_cusp_P_dim): ten spaces, k = 8, 11, 12, 14,
..., 20, of total dimension 11.  The lifts are those of the classical newforms of
level 2, weight 2k-2, trivial character and Fricke sign (-1)^k (root number -1),
as for every other paramodular level (Hecke_Eigenvalues/paramodular/
Hecke_Eigenvalues_Saito_Kurokawa.py, SK_query).

For each space with new_cusp_P_dim > 0 this builds the smf_newforms and smf_hecke_nf
rows of its lifts with that generator, numbered by (dim, trace_lambda_p) as in the
other spaces, and checks, before writing anything:
  - that no form of the space is stored;
  - that the dimensions of the lifts add up to new_cusp_P_dim;
  - that every lift has an smf_hecke_nf row (an mf_hecke_nf row of its classical form);
  - that the traces of the lifts add up to the stored cusp_P sums of the space
    (cusp_P_lambda_p, cusp_P_lambda_p_square), where those are not NULL;
  - that the new rows have the columns of the stored SK rows of family K.
A space whose new cusp space the lifts fill (new_cusp_dim = new_cusp_P_dim; there
new_cusp_G_dim is 0) gets num_forms, the number of its lifts, as the complete spaces
of the other levels; the others keep num_forms NULL.

The new rows are written to PREFIX_newforms.json and PREFIX_hecke_nf.json.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.add_level2_paramodular_sk [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.add_level2_paramodular_sk [--datafile PREFIX] --apply
"""
import argparse
import json

from sage.all import prime_range
from lmfdb import db
from smf_lmfdb.db_tables.common_populate import (
    base_26, common_entry_values, entry_add_common_columns, make_space_label)
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable
from smf_lmfdb.db_tables.add_level2_lifts import make_orbit_code, nulls, sort_key
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import (
    SK_query, SK_all_forms, SK_all_evs)

LEVEL = 2
FAMILY = 'K'
SUMS = {'cusp_P_lambda_p': 'trace_lambda_p', 'cusp_P_lambda_p_square': 'trace_lambda_p_square'}


def full_row(k, fields, hecke_orbit, cols):
    """A complete smf_newforms row, as smf_newforms_populate.entry_add_columns makes it."""
    e = common_entry_values(k, 0, LEVEL, FAMILY)
    e.update(fields)
    e = entry_add_common_columns(e, {'id': None})
    e['space_label'] = make_space_label(e)
    e['hecke_orbit'] = hecke_orbit
    e['label'] = e['space_label'] + '.' + base_26(hecke_orbit)
    e['hecke_orbit_code'] = make_orbit_code(2, FAMILY, LEVEL, k, 0, 1, hecke_orbit)
    e['trace_hash'] = e['analytic_rank'] = e['qexp_display'] = None
    e['analytic_rank_proved'] = False
    e['embedded_related_objects'] = []
    e['fricke_eigenval'] = 1
    for _, s in e['atkin_lehner_eigenvals']:
        e['fricke_eigenval'] *= s
    e['trace_display'] = e['trace_lambda_p'][:4]
    e['multiplicity'] = 1
    return {c: nulls(e.get(c)) for c in cols}


def full_ev(k, ev, row, cols):
    # 'an' (jsonb) keeps the string 'NULL' at the n not prime to the level, as in every stored row
    e = common_entry_values(k, 0, LEVEL, FAMILY)
    e.update(ev)
    e['label'], e['hecke_orbit_code'] = row['label'], row['hecke_orbit_code']
    if e.get('hecke_ring_power_basis'):
        # every stored row of family K gives its basis explicitly, also when it is the power
        # basis; the coordinates are the same
        d = len(e['field_poly']) - 1
        ident = [[int(i == j) for j in range(d)] for i in range(d)]
        e['hecke_ring_power_basis'] = False
        e['hecke_ring_numerators'] = e['hecke_ring_inverse_numerators'] = ident
        e['hecke_ring_denominators'] = e['hecke_ring_inverse_denominators'] = [1] * d
    return {c: e.get(c) if c == 'an' else nulls(e.get(c)) for c in cols}


def main(apply=False, datafile=None):
    base = datafile or "smf_level2_paramodular_sk"
    nf_cols = [c for c in db.smf_newforms.search_cols if c != 'id']
    hnf_cols = [c for c in db.smf_hecke_nf.search_cols if c != 'id']
    # the columns every stored SK row of family K has
    ref = list(db.smf_newforms.search({'family': FAMILY, 'aut_rep_type': 'P'}, nf_cols, limit=200))
    always = [c for c in nf_cols if all(r.get(c) is not None for r in ref)]
    ref_h = list(db.smf_hecke_nf.search({'label': {'$in': [r['label'] for r in ref]}}, hnf_cols))
    always_h = [c for c in hnf_cols if ref_h and all(h.get(c) is not None for h in ref_h)]
    new_rows, new_evs, num_forms = [], [], {}
    spaces = list(db.smf_newspaces.search({'level': LEVEL, 'family': FAMILY, 'new_cusp_P_dim': {'$gt': 0}},
                                          ['label', 'weight', 'new_cusp_P_dim', 'new_cusp_G_dim', 'new_cusp_dim', 'num_forms']
                                          + list(SUMS)))
    for s in spaces:
        k, j = s['weight']
        assert j == 0, s['label']
        assert db.smf_newforms.count({'space_label': s['label']}) == 0, "%s: forms are stored" % s['label']
        q = SK_query(k, LEVEL)
        forms = SK_all_forms(k, q, True)
        evs = SK_all_evs(k, q)
        assert sum(f['dim'] for f in forms) == s['new_cusp_P_dim'], \
            "%s: lifts of dimension %d, new_cusp_P_dim %d" % (s['label'], sum(f['dim'] for f in forms), s['new_cusp_P_dim'])
        assert len(evs) == len(forms), "%s: %d forms, %d smf_hecke_nf rows" % (s['label'], len(forms), len(evs))
        # SK_all_forms and SK_all_evs run through the same mf_newforms search
        pairs = sorted(zip(forms, evs), key=lambda fe: sort_key(fe[0]['dim'], fe[0]['trace_lambda_p']))
        rows = []
        for i, (f, ev) in enumerate(pairs):
            row = full_row(k, f, i + 1, nf_cols)
            missing = [c for c in always if row.get(c) is None]
            assert not missing, "%s: no value for %s" % (row['label'], missing)
            h = full_ev(k, ev, row, hnf_cols)
            missing = [c for c in always_h if h.get(c) is None]
            assert not missing, "%s: no smf_hecke_nf value for %s" % (row['label'], missing)
            rows.append(row)
            new_rows.append(row)
            new_evs.append(h)
        for sc, fc in SUMS.items():
            stored = s.get(sc)
            if stored is None:
                continue
            n = len(stored)
            got = [sum(r[fc][i] for r in rows) if all(r[fc][i] is not None for r in rows) else None for i in range(n)]
            agree = all(a is None or b is None or a == b for a, b in zip(stored, got))
            assert agree, "%s: %s differs from the sum over the lifts" % (s['label'], sc)
        complete = s.get('new_cusp_dim') == sum(r['dim'] for r in rows)
        if complete:
            assert s.get('num_forms') is None and not s.get('new_cusp_G_dim'), s['label']
            num_forms[s['label']] = len(rows)
        print("%s: %s, AL %s; new_cusp_G_dim %s, %s"
              % (s['label'], ', '.join("%s (dim %d, lifts %s)" % (r['label'], r['dim'], r['related_objects'][0].split('/')[-4:])
                                       for r in rows),
                 ' '.join(r['atkin_lehner_string'] for r in rows), s.get('new_cusp_G_dim'),
                 "complete: num_forms %d" % len(rows) if complete else "incomplete: num_forms stays NULL"))
    print("%d spaces, %d new smf_newforms rows, %d new smf_hecke_nf rows; num_forms set in %d spaces"
          % (len(spaces), len(new_rows), len(new_evs), len(num_forms)))
    with open(base + "_newforms.json", 'w') as F:
        json.dump(jsonable(new_rows), F)
    with open(base + "_hecke_nf.json", 'w') as F:
        json.dump(jsonable(new_evs), F)
    print("wrote %s_newforms.json and %s_hecke_nf.json" % (base, base))
    if not apply:
        print("dry run; pass --apply to write")
        return
    # delete by label to undo
    db.smf_newforms.insert_many(jsonable(new_rows))
    db.smf_hecke_nf.insert_many(jsonable(new_evs))    # psycopg cannot adapt Sage numbers
    for label, n in num_forms.items():
        db.smf_newspaces.update({'label': label}, {'num_forms': n}, restat=False)
    for t in (db.smf_newforms, db.smf_hecke_nf, db.smf_newspaces):
        t.stats.refresh_stats(total=False)
    print("inserted; num_forms set")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--datafile", help="prefix for the output files (default: ./smf_level2_paramodular_sk)")
    parser.add_argument("--apply", action="store_true", help="write to the database")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
