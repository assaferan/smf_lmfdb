"""
Add the missing lifts to the spaces 2.P.2.k.j.a (Gamma[2], k >= 3), and give
every lift its multiplicity (smf_newforms.multiplicity).

Each space stored only the lifts of multiplicity 1 (the Yoshida lifts of level-2
pairs of opposite signs) and the Saito-Kurokawa lifts of the level-2 forms with
w_2 = +1, although its dimension columns count all of them.  The eigensystems,
their multiplicities and the conventions are in Hecke_Eigenvalues/Gamma2_lifts.py.
A form of multiplicity m is stored once, with `multiplicity` m.

In each space the rows are numbered in the order of the level-1 spaces: by type
(F, Q, P, Y), then by (dim, trace_lambda_p).  Adding rows can therefore relabel
the stored ones; their smf_hecke_nf rows follow.  Spaces of weight k = 2, whose
dimensions are conjectural, are left alone.

Before writing anything it checks, in each space,
- that the sum of dim * multiplicity over each type is the stored new_*_dim;
- that every stored row is one of the eigensystems, and that regenerating it and
  its smf_hecke_nf row gives what is stored (the columns that differ are
  reported, and the space is left alone unless --refresh lists them);
- that the stored labels, hecke_orbit_codes and smf_hecke_nf rows are consistent.
Spaces that fail a check are reported and left alone.

The new rows are written to PREFIX_newforms.json and PREFIX_hecke_nf.json, the
relabelling to PREFIX_newforms_update.txt and PREFIX_hecke_nf_update.txt.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.add_level2_lifts [--datafile PREFIX] [--spaces LABEL,...]          # dry run
    sage -python -m smf_lmfdb.db_tables.add_level2_lifts [--datafile PREFIX] --apply
"""
import argparse
import sys
import time
import json
import traceback
from collections import Counter, defaultdict

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import (
    base_26, common_entry_values, entry_add_common_columns, make_space_label)
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable
from smf_lmfdb.Hecke_Eigenvalues.Gamma2_lifts import (
    LEVEL, TYPE_ORDER, check_sum, eigensystems, lift_data, sources)

DIM_COLS = {'F': 'new_eis_F_dim', 'Q': 'new_eis_Q_dim', 'P': 'new_cusp_P_dim', 'Y': 'new_cusp_Y_dim'}
RELABEL = ['label', 'hecke_orbit', 'hecke_orbit_code']
# columns not compared between a stored row and its regeneration
NOT_COMPARED = {'id', 'label', 'hecke_orbit', 'hecke_orbit_code', 'multiplicity'}


def sort_key(dim, ap):
    # as in fix_degree2_traces (not imported: it imports smf_newforms_populate)
    return [dim] + [float('-inf') if x is None else x for x in ap]


def make_orbit_code(g, F, N, k, j, i, X):
    # as in smf_newforms_populate (not imported: that module loads the .sage code)
    return g + (ord(F) << 8) + (N << 16) + (k << 32) + (j << 40) + ((i-1) << 48) + ((X-1) << 56)


def nulls(x):
    """The populate code's 'NULL' as None, also inside lists."""
    if isinstance(x, list):
        return [nulls(y) for y in x]
    return None if x == 'NULL' else x


def full_row(k, j, fields, hecke_orbit, mult, cols):
    """A complete smf_newforms row, as smf_newforms_populate.entry_add_columns makes it."""
    e = common_entry_values(k, j, LEVEL, 'P')
    e.update(fields)
    e = entry_add_common_columns(e, {'id': None})
    e['space_label'] = make_space_label(e)
    e['hecke_orbit'] = hecke_orbit
    e['label'] = e['space_label'] + '.' + base_26(hecke_orbit)
    e['hecke_orbit_code'] = make_orbit_code(2, 'P', LEVEL, k, j, 1, hecke_orbit)
    e['trace_hash'] = e['analytic_rank'] = e['qexp_display'] = None
    e['analytic_rank_proved'] = False
    e['embedded_related_objects'] = []
    if not e.get('atkin_lehner_eigenvals') and 'atkin_lehner_eigenvals' not in fields:
        e['atkin_lehner_eigenvals'] = [[2, 1]]
        e['atkin_lehner_string'] = '+'
    e['fricke_eigenval'] = 1
    for _, s in e['atkin_lehner_eigenvals']:
        e['fricke_eigenval'] *= s
    e['trace_display'] = e['trace_lambda_p'][:4]
    e['multiplicity'] = mult
    return {c: nulls(e.get(c)) for c in cols}


def full_ev(k, j, ev, row, cols):
    # 'an' (jsonb) keeps the string 'NULL' at the n not prime to the level, as in every stored row
    e = common_entry_values(k, j, LEVEL, 'P')
    e.update(ev)
    e['label'], e['hecke_orbit_code'] = row['label'], row['hecke_orbit_code']
    return {c: e.get(c) if c == 'an' else nulls(e.get(c)) for c in cols}


def differing(a, b, cols):
    return sorted(c for c in cols if c not in NOT_COMPARED and jsonable(a.get(c)) != jsonable(b.get(c)))


def space_plan(s, stored, hnf, nf_cols, hnf_cols, refresh):
    """(new rows, new hecke_nf rows, relabelled rows, relabelled hecke_nf rows, report) for one space."""
    k, j = s['weight']
    E = eigensystems(k, j)
    sums = check_sum(E)
    for t, c in DIM_COLS.items():
        assert sums.get(t, 0) == (s.get(c) or 0), "%s: %s = %s, the eigensystems give %s" % (s['label'], c, s.get(c), sums.get(t, 0))
    gen = []
    for typ, mult, data in E:
        for fields, ev in lift_data(k, j, data):
            gen.append({'typ': typ, 'mult': mult, 'fields': fields, 'ev': ev, 'src': sources(data)})
    gen.sort(key=lambda g: [TYPE_ORDER.index(g['typ'])] + sort_key(g['fields']['dim'], g['fields']['trace_lambda_p']))
    keys = [[TYPE_ORDER.index(g['typ'])] + sort_key(g['fields']['dim'], g['fields']['trace_lambda_p']) for g in gen]
    assert all(a != b for a, b in zip(keys, keys[1:])), "%s: two eigensystems have the same sort key" % s['label']
    # the stored rows, by type and the forms they lift
    by_src = {}
    for r in stored:
        key = (r['aut_rep_type'], tuple(r['related_objects'] or []))
        assert key not in by_src, "%s: two stored rows lift %s" % (s['label'], key)
        by_src[key] = r
    new_rows, new_evs, upd, upd_ev, diffs = [], [], [], [], Counter()
    matched = set()
    for i, g in enumerate(gen, 1):
        row = full_row(k, j, g['fields'], i, g['mult'], nf_cols)
        key = (g['typ'], tuple(row['related_objects']))
        old = by_src.get(key)
        # the orbits of a split Yoshida lift share their key; none of them is stored
        assert old is None or sum(1 for x in gen if x['typ'] == g['typ'] and x['src'] == g['src']) == 1, key
        if old is None:
            new_rows.append(row)
            if g['ev'] is not None:
                new_evs.append(full_ev(k, j, g['ev'], row, hnf_cols))
            continue
        matched.add(key)
        d = differing(old, row, nf_cols)
        for c in d:
            diffs[(g['typ'], c)] += 1
        bad = [c for c in d if c not in refresh]
        assert not bad, "%s: the stored %s differs from its regeneration in %s" % (s['label'], old['label'], bad)
        u = {'id': old['id'], 'multiplicity': g['mult']}
        u.update({c: row[c] for c in RELABEL + d})
        upd.append(u)
        h = hnf.get(old['hecke_orbit_code'])
        assert h is None or g['ev'] is not None, "%s has an smf_hecke_nf row but no regeneration" % old['label']
        if h is None and g['ev'] is not None:
            new_evs.append(full_ev(k, j, g['ev'], row, hnf_cols))    # the stored row had none
        if h is not None:
            assert h['label'] == old['label'], (h['label'], old['label'])
            new_ev = full_ev(k, j, g['ev'], row, hnf_cols)
            dh = differing(h, new_ev, hnf_cols)
            for c in dh:
                diffs[(g['typ'], 'hecke_nf.' + c)] += 1
            bad = [c for c in dh if 'hecke_nf.' + c not in refresh]
            assert not bad, "%s: its smf_hecke_nf row differs from the regeneration in %s" % (old['label'], bad)
            ue = {'id': h['id'], 'label': row['label'], 'hecke_orbit_code': row['hecke_orbit_code']}
            ue.update({c: new_ev[c] for c in dh})
            upd_ev.append(ue)
    missing = set(by_src) - matched
    assert not missing, "%s: stored rows that are no lift eigensystem: %s" % (s['label'], [by_src[m]['label'] for m in missing])
    return new_rows, new_evs, upd, upd_ev, diffs


def main(apply=False, datafile=None, spaces=None, refresh=()):
    base = datafile or "smf_level2_lifts"
    nf_cols = sorted(c for c in db.smf_newforms.search_cols if c != 'id') + ['id']
    has_mult = 'multiplicity' in db.smf_newforms.search_cols
    if not has_mult:
        print("smf_newforms has no column multiplicity yet: it is in the files, and --apply needs it")
        nf_cols = [c for c in nf_cols if c != 'id'] + ['multiplicity', 'id']
    hnf_cols = sorted(c for c in db.smf_hecke_nf.search_cols if c != 'id')
    query = {'level': LEVEL, 'family': 'P'}
    if spaces:
        query['label'] = {'$in': spaces}
    stored = defaultdict(list)
    for r in db.smf_newforms.search({'level': LEVEL, 'family': 'P'}, [c for c in nf_cols if c in db.smf_newforms.search_cols or c == 'id']):
        stored[r['space_label']].append(r)
    codes = [r['hecke_orbit_code'] for L in stored.values() for r in L]
    hnf = {h['hecke_orbit_code']: h for h in db.smf_hecke_nf.search({'hecke_orbit_code': {'$in': codes}}, ['id'] + hnf_cols)}
    new_rows, new_evs, upd, upd_ev = [], [], [], []
    diffs, skipped, report = Counter(), [], Counter()
    for s in db.smf_newspaces.search(query, ['label', 'weight'] + list(DIM_COLS.values())):
        k, j = s['weight']
        if k < 3:
            report['k = 2 (conjectural), left alone'] += 1
            continue
        t0 = time.time()
        try:
            plan = space_plan(s, stored[s['label']], hnf, nf_cols, hnf_cols, set(refresh))
            sys.stderr.write("%s planned in %.0fs\n" % (s['label'], time.time() - t0))
            sys.stderr.flush()
        except (AssertionError, NotImplementedError, ValueError) as err:
            skipped.append((s['label'], repr(err)))
            if not isinstance(err, (AssertionError, NotImplementedError)):
                traceback.print_exc()
            continue
        nr, ne, u, ue, d = plan
        new_rows += nr
        new_evs += ne
        upd += u
        upd_ev += ue
        diffs += d
        report['spaces planned'] += 1
        for r in nr:
            report['new %s rows' % r['aut_rep_type']] += 1
        report['relabelled rows'] += sum(1 for x in u if x['label'] != next(r['label'] for r in stored[s['label']] if r['id'] == x['id']))
    for x, n in sorted(report.items()):
        print("%6d %s" % (n, x))
    print("%d new smf_newforms rows, %d new smf_hecke_nf rows; %d stored rows updated (%d smf_hecke_nf)"
          % (len(new_rows), len(new_evs), len(upd), len(upd_ev)))
    if diffs:
        print("stored rows that differ from their regeneration, by (type, column):", dict(diffs))
    if skipped:
        print("%d spaces left alone:" % len(skipped))
        for x in skipped:
            print("   ", x)
    with open(base + "_newforms.json", 'w') as F:
        json.dump(jsonable(new_rows), F)
    with open(base + "_hecke_nf.json", 'w') as F:
        json.dump(jsonable(new_evs), F)
    print("wrote %s_newforms.json and %s_hecke_nf.json" % (base, base))
    upd_cols = sorted({c for u in upd for c in u} - {'id'})
    for u in upd:
        for c in upd_cols:
            u.setdefault(c, None)    # \N: keep the stored value
    if has_mult and upd:
        write_update_file(base + "_newforms_update.txt", db.smf_newforms, 'id', upd, upd_cols)
    elif upd:
        with open(base + "_newforms_update.json", 'w') as F:
            json.dump(jsonable(upd), F)
        print("wrote %s_newforms_update.json (the update file needs the column multiplicity)" % base)
    ev_cols = sorted({c for u in upd_ev for c in u} - {'id'})
    for u in upd_ev:
        for c in ev_cols:
            u.setdefault(c, None)
    if upd_ev:
        write_update_file(base + "_hecke_nf_update.txt", db.smf_hecke_nf, 'id', upd_ev, ev_cols)
    if not apply:
        print("dry run; pass --apply to write")
        return
    assert has_mult, "add the column first: db.smf_newforms.add_column('multiplicity', 'smallint')"
    assert not skipped or spaces, "some spaces fail a check; pass --spaces to apply the others"
    # update_from_file builds replacement tables; undo with reload_revert().  Relabel first,
    # so that no new row shares a label with a stored one.
    if upd:
        db.smf_newforms.update_from_file(base + "_newforms_update.txt", label_col='id')
    if upd_ev:
        db.smf_hecke_nf.update_from_file(base + "_hecke_nf_update.txt", label_col='id')
    # The new rows: delete them by label to undo.
    if new_rows:
        db.smf_newforms.insert_many([{c: v for c, v in r.items() if c != 'id'} for r in new_rows])
    if new_evs:
        db.smf_hecke_nf.insert_many(new_evs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write to the database")
    parser.add_argument("--datafile", help="prefix for the output files (default: ./smf_level2_lifts)")
    parser.add_argument("--spaces", help="comma-separated space labels (default: all of level 2, family P)")
    parser.add_argument("--refresh", help="comma-separated columns of the stored rows to replace by their regeneration "
                        "(hecke_nf.COL for a column of smf_hecke_nf)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile,
         spaces=args.spaces.split(',') if args.spaces else None,
         refresh=args.refresh.split(',') if args.refresh else ())
