"""
Add the paramodular newforms of a weight (k, j) other than (3, 0) from omf5_data
(hecke_evs_4_0, hecke_evs_3_2): the forms of type G of omf5_data and, for j = 0, the
Saito-Kurokawa lifts, as smf_newforms_populate and smf_hecke_nf_populate make them for
(3, 0), into the spaces 2.K.N.k.j.a.

A level is done only if the whole space is filled (no partly populated spaces):
  - the new G dimension of omf5_data is the space's new_cusp_G_dim;
  - the Saito-Kurokawa lifts have dimension new_cusp_P_dim, and their Atkin-Lehner
    signs the dimensions of Saito_Kurokawa_new_lift_dim;
  - the forms have dimension new_cusp_dim (new_cusp_Y_dim is 0);
and then the space gets what smf_newspaces_populate gives the spaces of (3, 0): num_forms,
ALdims_G, ALdims_P, ALdims and the traces of the Hecke eigenvalues of each type
(Hecke_Eigenvalues_Traces_paramodular).  Otherwise the level is skipped and listed.

The orbits are numbered as for (3, 0): the G orbits by (dim, trace_lambda_p), then the
Saito-Kurokawa lifts.  Every row gets multiplicity 1, as the stored paramodular rows.
A level that already has stored forms (31 in weight (4, 0), 19 in (3, 2)) is left as it
is, since its rows have data omf5_data does not (values at the bad primes, T(p^2) for
p <= 29, q-expansions); its stored forms must be the regenerated ones (labels, types,
dimensions and the traces of T(p) at the good primes).  The values at the primes
dividing the level are not in omf5_data: they stay NULL, for add_bad_primes.py --weight.

The rows are written to PREFIX_newforms.json, PREFIX_hecke_nf.json and PREFIX_updates.json
(the changes to stored rows and to the spaces), and the num_forms to PREFIX_num_forms.json.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.add_paramodular_weights --weight K,J [--datafile PREFIX] [N ...]            # dry run
    sage -python -m smf_lmfdb.db_tables.add_paramodular_weights --weight K,J [--datafile PREFIX] [N ...] --apply
Without levels, every level of the weight's folder of omf5_data.

After --apply, run fill_trace_tables.py, which adds the new forms to smf_hecke_traces and
smf_hecke_newspace_traces (the trace search pages).
"""
import argparse
import json
import os
import re

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import common_entry_values, fill_nulls, MAX_RING_DEGREE
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable
from smf_lmfdb.db_tables.add_level2_lifts import nulls
from smf_lmfdb.db_tables import smf_newforms_populate, smf_hecke_nf_populate
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import (
    Hecke_Eigenforms_paramodular, Hecke_Eigenvalues_paramodular_G, Hecke_Eigenvalues_Traces_paramodular,
    al_str_to_num, weight_folder)
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import SK_all_evs
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import Saito_Kurokawa_new_lift_dim, cmf_ALdims_index
from sage.all import divisors, is_squarefree


def levels_of(k, j):
    folder = os.path.join(weight_folder(k, j), 'data')
    pat = re.compile(r'hecke_ev_%d_%d_(\d+)\.dat$' % (k, j))
    return sorted(int(m.group(1)) for m in map(pat.match, os.listdir(folder)) if m)


def ring_basis(e):
    """As common_populate.write_data: an explicit basis for the power basis; None above MAX_RING_DEGREE."""
    if 'hecke_ring_numerators' in e:
        return e
    d = len(e['field_poly']) - 1
    if e.get('hecke_ring_power_basis') or d == 1:
        ident = [[int(r == c) for c in range(d)] for r in range(d)]
        e['hecke_ring_power_basis'] = False
        e['hecke_ring_numerators'] = e['hecke_ring_inverse_numerators'] = ident
        e['hecke_ring_denominators'] = e['hecke_ring_inverse_denominators'] = [1] * d
        return e
    if d > MAX_RING_DEGREE:
        return None
    raise ValueError("%s: an smf_hecke_nf entry without a Hecke ring basis" % e.get('label'))


def level_rows(k, j, N):
    """(smf_newforms rows, smf_hecke_nf rows) of the level, numbered as the populate code does."""
    nf_types = db.smf_newforms.col_type
    hnf_types = db.smf_hecke_nf.col_type
    entry = common_entry_values(k, j, N, 'K')
    forms = Hecke_Eigenforms_paramodular(k, j, N)
    rows = []
    for i, f in enumerate(forms, 1):
        e = entry.copy()
        e.update(f)
        e = smf_newforms_populate.entry_add_columns(e, {'id': None, 'num_forms': i})
        e['multiplicity'] = 1
        e = fill_nulls(e, db.smf_newforms)
        rows.append({c: nulls(e.get(c)) for c in nf_types if c != 'id'})
    # the entries in the order of the forms: the G orbits (Hecke_Eigenvalues_paramodular_G keeps
    # the order of parse_omf5), then each Saito-Kurokawa lift from its classical form (SK_all_evs
    # skips a form without an mf_hecke_nf row, so Hecke_Eigenvalues_PM_SK_all_evs would shift them)
    evs = Hecke_Eigenvalues_paramodular_G(k, j, N)
    nG = sum(1 for f in forms if f['aut_rep_type'] == 'G')
    assert len(evs) == nG and all(f['aut_rep_type'] == 'G' for f in forms[:nG]), N
    for f in forms[nG:]:
        assert f['aut_rep_type'] == 'P' and len(f['related_objects']) == 1, N
        cmf = '.'.join(f['related_objects'][0].split('/')[-4:])
        ev = SK_all_evs(k, {'label': cmf})
        assert len(ev) <= 1, cmf
        evs.append(ev[0] if ev else None)
    hrows = []
    for i, ev in enumerate(evs, 1):
        if ev is None:
            continue
        e = entry.copy()
        e.update(ev)
        e = ring_basis(e)
        if e is None:
            continue
        e = smf_hecke_nf_populate.entry_add_columns(e, {'id': None, 'num_forms': i})
        e = fill_nulls(e, db.smf_hecke_nf)
        r = {c: (e.get(c) if c == 'an' else nulls(e.get(c))) for c in hnf_types if c != 'id'}
        # the newform of the same number
        f = rows[i - 1]
        assert r['label'] == f['label'] and r['field_poly'] == f['field_poly'], (N, i)
        hrows.append(r)
    return rows, hrows


def same(a, b):
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if a is None or b is None:
        return a is None and b is None
    try:
        return float(a) == float(b) and int(a) == int(b)
    except (TypeError, ValueError, OverflowError):
        return str(a) == str(b)


def plan_level(k, j, N):
    label = '2.K.%d.%d.%d.a' % (N, k, j)
    sp = db.smf_newspaces.lucky({'label': label}, ['id'] + list(db.smf_newspaces.search_cols))
    if sp is None:
        return None, 'no space'
    rows, hrows = level_rows(k, j, N)
    dims = {t: sum(r['dim'] for r in rows if r['aut_rep_type'] == t) for t in 'GP'}
    if dims['G'] != sp.get('new_cusp_G_dim'):
        return None, 'G dimension %d, new_cusp_G_dim %s' % (dims['G'], sp.get('new_cusp_G_dim'))
    if dims['P'] != (sp.get('new_cusp_P_dim') or 0):
        return None, 'Saito-Kurokawa dimension %d, new_cusp_P_dim %s' % (dims['P'], sp.get('new_cusp_P_dim'))
    if sp.get('new_cusp_Y_dim'):
        return None, 'new_cusp_Y_dim %s' % sp.get('new_cusp_Y_dim')
    if dims['G'] + dims['P'] != sp.get('new_cusp_dim'):
        return None, 'dimension %d, new_cusp_dim %s' % (dims['G'] + dims['P'], sp.get('new_cusp_dim'))
    # the space columns, as smf_newspaces_populate
    traces, dim_G_new, al_dims_G = Hecke_Eigenvalues_Traces_paramodular(k, j, N)
    assert dim_G_new == dims['G']
    al_dims_P = [0] * len(al_dims_G)
    for d in divisors(N):
        if j == 0 and is_squarefree(d):
            al_dims_P[cmf_ALdims_index(d, N)] = int(Saito_Kurokawa_new_lift_dim(k, N, al=d))
    al_rows_P = [0] * len(al_dims_G)
    for r in rows:
        if r['aut_rep_type'] == 'P':
            al_rows_P[cmf_ALdims_index(al_str_to_num(r['atkin_lehner_string'], N), N)] += r['dim']
    if al_rows_P != al_dims_P:
        return None, 'Saito-Kurokawa Atkin-Lehner dimensions %s, Saito_Kurokawa_new_lift_dim %s' % (al_rows_P, al_dims_P)
    space_set = {c: nulls(v) for c, v in traces.items()}
    space_set.update({'num_forms': len(rows), 'ALdims_G': [int(x) for x in al_dims_G], 'ALdims_P': al_dims_P,
                      'ALdims': [int(x) + y for x, y in zip(al_dims_G, al_dims_P)]})
    space_set = {c: v for c, v in space_set.items() if not same(v, sp.get(c))}
    stored = {r['label']: r for r in db.smf_newforms.search({'space_label': label},
                                                            ['id', 'label', 'aut_rep_type', 'dim', 'trace_lambda_p'])}
    if stored:
        by_label = {r['label']: r for r in rows}
        for lab, s in stored.items():
            r = by_label.get(lab)
            assert r is not None, "%s: stored, but not regenerated" % lab
            assert (s['aut_rep_type'], s['dim']) == (r['aut_rep_type'], r['dim']), lab
            assert all(t is None or same(t, u) for t, u in zip(r['trace_lambda_p'], s['trace_lambda_p'])), \
                "%s: the stored traces of T(p) are not those of omf5_data" % lab
        assert len(stored) == len(rows), "%s: %d stored forms, %d regenerated" % (label, len(stored), len(rows))
        return None, 'already stored, checked'
    new, new_h, updates = rows, hrows, []
    if space_set:
        updates.append({'table': 'smf_newspaces', 'id': sp['id'], 'label': label, 'set': space_set})
    return {'space': label, 'new': new, 'new_h': new_h, 'updates': updates, 'num_forms': len(rows),
            'stored_num_forms': sp.get('num_forms'), 'dims': dims}, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--weight', required=True, help='k,j')
    ap.add_argument('--datafile', default='smf_paramodular_weights')
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('levels', nargs='*', type=int)
    a = ap.parse_args()
    k, j = map(int, a.weight.split(','))
    assert (k, j) != (3, 0) and weight_folder(k, j), "no omf5_data folder for the weight %s" % a.weight
    levels = a.levels or levels_of(k, j)
    plans, skipped = [], {}
    for N in levels:
        plan, why = plan_level(k, j, N)
        if plan is None:
            skipped[N] = why
            if why != 'no space':
                print("level %d: %s" % (N, why if why.startswith('already') else 'SKIPPED, ' + why), flush=True)
            continue
        plans.append(plan)
        print("%s: %d new forms (%d hecke_nf), dims %s, space columns %s, num_forms %s -> %d"
              % (plan['space'], len(plan['new']), len(plan['new_h']), plan['dims'],
                 sorted(c for u in plan['updates'] for c in u['set']), plan['stored_num_forms'], plan['num_forms']),
              flush=True)
    new = [r for p in plans for r in p['new']]
    new_h = [r for p in plans for r in p['new_h']]
    updates = [u for p in plans for u in p['updates']]
    num_forms = {p['space']: p['num_forms'] for p in plans}
    print("weight (%d,%d): %d levels, %d new smf_newforms rows, %d new smf_hecke_nf rows, %d updates, num_forms in %d spaces; "
          "%d levels without a space, %d skipped"
          % (k, j, len(plans), len(new), len(new_h), len(updates), len(num_forms),
             sum(1 for w in skipped.values() if w == 'no space'),
             sum(1 for w in skipped.values() if w not in ('no space', 'already stored, checked'))))
    for name, data in (('newforms', new), ('hecke_nf', new_h), ('updates', updates), ('num_forms', num_forms)):
        with open('%s_%s.json' % (a.datafile, name), 'w') as F:
            json.dump(jsonable(data), F)
    if not a.apply:
        print("dry run; pass --apply to write")
        return
    for u in updates:
        db[u['table']].update({'id': u['id']}, jsonable(u['set']), restat=False)
    # the new rows: delete them by label to undo
    if new:
        db.smf_newforms.insert_many(jsonable(new), restat=False)
    if new_h:
        db.smf_hecke_nf.insert_many(jsonable(new_h), restat=False)
    for t in (db.smf_newforms, db.smf_hecke_nf, db.smf_newspaces):
        t.stats.refresh_stats(total=False)
    print("written")


if __name__ == '__main__':
    main()
