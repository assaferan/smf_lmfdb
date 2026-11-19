"""
Add the spaces 2.P.2.2.j.a (Gamma[2], weight (2, j), j even, j <= JMAX) and their lifts.

No space of family P, level 2 and k = 2 is stored, though three of its forms are
(2.P.2.2.{12,18,20}.a.a, the s[1^6] Yoshida lifts).  The cusp part is
Clery-van der Geer (arXiv:1709.01748): Thm 1.1 gives its Yoshida lifts (with their
S_6 multiplicities; Hecke_Eigenvalues/Gamma2_lifts.py, eigensystems), Prop. 1.5
gives S_{j,2}(Gamma[2]) = 0 for j < 12, and Conj. 1.2 says that the lifts span it.
So cusp_dim = cusp_Y_dim is the dimension of the lifts, proved for j < 12 and
conjectural for j >= 12.  The level-1 part is 0 (S_{j,2}(Gamma_2) = 0 for j <= 52,
Thm 1.3, and M_{j,2}(Gamma_2) = S_{j,2}(Gamma_2)), so every dimension is new.

The non-cusp part:
  j = 0: M_2(Gamma[2]) has dimension 5 (Igusa; the s[2^3] spanned by the theta
         fourth powers) and S_2(Gamma[2]) = 0: eis_F_dim = 5, with no form stored;
  j > 0: --eis-j sets its dimension (default 0, as M_{j,2} = S_{j,2} at level 1).

dims_proved is False where a dimension is conjectural: j >= 12, and j > 0 unless
--eis-proved.

For j > 0 the lifts go through add_level2_lifts.space_plan: new rows for the lifts
not stored, and the stored ones relabelled and refreshed (--refresh) as there.
num_forms is set where the forms fill the space (eis 0), else NULL.

The new space rows are written to PREFIX_newspaces.json, the forms as in
add_level2_lifts.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.add_level2_weight2 [--datafile PREFIX] [--refresh COLS]          # dry run
    sage -python -m smf_lmfdb.db_tables.add_level2_weight2 [--datafile PREFIX] [--refresh COLS] --apply
"""
import argparse
import json
from collections import Counter

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import common_entry_values, entry_add_common_columns, make_space_label
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.db_tables.add_level2_lifts import DIM_COLS, make_orbit_code, nulls, space_plan
from smf_lmfdb.Hecke_Eigenvalues.Gamma2_lifts import check_sum, eigensystems, save_polredabs_cache

K = 2
LEVEL = 2
JMAX = 20
TYPES = ['eis_F', 'eis_Q', 'cusp_P', 'cusp_Y', 'cusp_G']


def space_row(j, eis_j, eis_proved, cols):
    E = eigensystems(K, j)
    sums = check_sum(E)
    assert set(sums) <= {'Y'}, (j, sums)
    d = {t: 0 for t in TYPES}
    d['cusp_Y'] = int(sums.get('Y', 0))
    if j < 12:
        assert d['cusp_Y'] == 0, "j = %d: lifts of dimension %d, but S_{j,2}(Gamma[2]) = 0 (Prop. 1.5)" % (j, d['cusp_Y'])
    if j == 0:
        d['eis_F'] = 5
    else:
        d['eis_Q'] = eis_j
    e = common_entry_values(K, j, LEVEL, 'P')
    e = entry_add_common_columns(e, {'id': None})
    e['label'] = make_space_label(e)
    e['hecke_orbit_code'] = make_orbit_code(2, 'P', LEVEL, K, j, 1, 1)
    for t, v in d.items():
        e[t + '_dim'] = e['new_' + t + '_dim'] = v
        e['old_' + t + '_dim'] = 0
    for part, ts in (('eis', TYPES[:2]), ('cusp', TYPES[2:]), ('total', TYPES)):
        e[part + '_dim'] = e['new_' + part + '_dim'] = sum(d[t] for t in ts)
        e['old_' + part + '_dim'] = 0
    e['dims_proved'] = j == 0 or (j < 12 and eis_proved)
    e['num_forms'] = None
    return E, {c: nulls(e.get(c)) for c in cols}


def main(apply=False, datafile=None, refresh=(), eis_j=0, eis_proved=False, jmax=JMAX):
    base = datafile or "smf_level2_weight2"
    ns_cols = [c for c in db.smf_newspaces.search_cols if c != 'id']
    has_proved = 'dims_proved' in ns_cols
    if not has_proved:
        print("smf_newspaces has no column dims_proved yet: it is in the files, and --apply needs it")
        ns_cols.append('dims_proved')
    nf_cols = sorted(c for c in db.smf_newforms.search_cols if c != 'id') + ['id']
    hnf_cols = sorted(c for c in db.smf_hecke_nf.search_cols if c != 'id')
    stored_spaces = set(db.smf_newspaces.search({'level': LEVEL, 'family': 'P', 'weight': {'$in': [[K, j] for j in range(0, jmax + 1, 2)]}}, 'label'))
    assert not stored_spaces, "spaces already stored: %s" % sorted(stored_spaces)
    stored = list(db.smf_newforms.search({'level': LEVEL, 'family': 'P', 'space_label': {'$like': '2.P.2.2.%'}},
                                         [c for c in nf_cols if c in db.smf_newforms.search_cols or c == 'id']))
    hnf = {h['hecke_orbit_code']: h for h in db.smf_hecke_nf.search(
        {'hecke_orbit_code': {'$in': [r['hecke_orbit_code'] for r in stored]}}, ['id'] + hnf_cols)}
    spaces, new_rows, new_evs, upd, upd_ev, diffs = [], [], [], [], [], Counter()
    for j in range(0, jmax + 1, 2):
        E, s = space_row(j, eis_j, eis_proved, ns_cols)
        mine = [r for r in stored if r['space_label'] == s['label']]
        rows = []
        if j > 0:
            nr, ne, u, ue, d = space_plan(s, mine, hnf, nf_cols, hnf_cols, set(refresh))
            save_polredabs_cache()
            new_rows += nr
            new_evs += ne
            upd += u
            upd_ev += ue
            diffs += d
            rows = nr + u
            if s['new_eis_dim'] == 0:
                s['num_forms'] = len(rows)
        else:
            assert not mine, s['label']
        spaces.append(s)
        print("%s: total %d (eis %d, cusp %d = Y %d), dims_proved %s, num_forms %s; %d new forms, %d stored; %s"
              % (s['label'], s['total_dim'], s['eis_dim'], s['cusp_dim'], s['cusp_Y_dim'], s['dims_proved'], s['num_forms'],
                 len(rows) - len(mine), len(mine),
                 ', '.join("mult %d x %s" % (m, '+'.join(x['label'].split('.')[-1] if isinstance(x, dict) and 'label' in x else '?'
                                                        for x in data[1:])) for _, m, data in E) or 'no lifts'))
    print("%d new spaces, %d new smf_newforms rows, %d new smf_hecke_nf rows; %d stored rows updated (%d smf_hecke_nf)"
          % (len(spaces), len(new_rows), len(new_evs), len(upd), len(upd_ev)))
    if diffs:
        print("stored rows that differ from their regeneration, by (type, column):", dict(diffs))
    for name, data in (('newspaces', spaces), ('newforms', new_rows), ('hecke_nf', new_evs), ('newforms_update', upd), ('hecke_nf_update', upd_ev)):
        with open("%s_%s.json" % (base, name), 'w') as F:
            json.dump(jsonable(data), F)
    print("wrote %s_{newspaces,newforms,hecke_nf,newforms_update,hecke_nf_update}.json" % base)
    if not apply:
        print("dry run; pass --apply to write")
        return
    assert has_proved, "add the column first: db.smf_newspaces.add_column('dims_proved', 'boolean')"
    upd_cols = sorted({c for u in upd for c in u} - {'id'})
    for u in upd:
        for c in upd_cols:
            u.setdefault(c, None)
    ev_cols = sorted({c for u in upd_ev for c in u} - {'id'})
    for u in upd_ev:
        for c in ev_cols:
            u.setdefault(c, None)
    # update_from_file builds replacement tables (undo with reload_revert); relabel first.
    if upd:
        write_update_file(base + "_newforms_update.txt", db.smf_newforms, 'id', upd, upd_cols)
        db.smf_newforms.update_from_file(base + "_newforms_update.txt", label_col='id')
    if upd_ev:
        write_update_file(base + "_hecke_nf_update.txt", db.smf_hecke_nf, 'id', upd_ev, ev_cols)
        db.smf_hecke_nf.update_from_file(base + "_hecke_nf_update.txt", label_col='id')
    # the new rows: delete them by label to undo
    db.smf_newspaces.insert_many(jsonable(spaces))
    if new_rows:
        db.smf_newforms.insert_many(jsonable([{c: v for c, v in r.items() if c != 'id'} for r in new_rows]))
    if new_evs:
        db.smf_hecke_nf.insert_many(jsonable(new_evs))
    for t in (db.smf_newspaces, db.smf_newforms, db.smf_hecke_nf):
        t.stats.refresh_stats(total=False)
    print("inserted")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write to the database")
    parser.add_argument("--datafile", help="prefix for the output files (default: ./smf_level2_weight2)")
    parser.add_argument("--refresh", help="as in add_level2_lifts")
    parser.add_argument("--eis-j", type=int, default=0, help="dimension of the non-cusp part for j > 0 (default 0)")
    parser.add_argument("--eis-proved", action="store_true", help="that dimension is proved")
    parser.add_argument("--jmax", type=int, default=JMAX)
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile, refresh=args.refresh.split(',') if args.refresh else (),
         eis_j=args.eis_j, eis_proved=args.eis_proved, jmax=args.jmax)
