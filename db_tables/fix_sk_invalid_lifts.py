"""
Remove the Saito-Kurokawa rows that are not lifts, store the lift that was
missing in their place, and mark the level 1 and 2 lifts as cusp forms.

A classical newform f of weight 2k-2 has a Saito-Kurokawa lift only if its root
number (-1)^(k-1) * (Fricke sign) is -1.  The old generator took Fricke sign -1
for every k, which is right only for k odd, and imposed nothing at level 1, so

  * at level 1 it "lifted" the forms of weight 2k-2 for k odd (k = 7, ..., 19,
    in each of the families K, S and P): 21 rows in spaces of dimension 0,
    12 of them with smf_hecke_nf rows;
  * at level 31 and k = 4 it lifted 31.6.a.b (Fricke -1, dimension 8) instead
    of 31.6.a.a (Fricke +1, dimension 5), as 2.K.31.4.0.a.b, whose space has
    new_cusp_P_dim 5.

The SK rows of levels 1 and 2 also have is_cuspidal = False.

This script

  * deletes the 21 level-1 rows from smf_newforms and their smf_hecke_nf rows,
    after saving them to a JSON file, and sets num_forms = 0 and the cusp_P
    sums to 0 in their spaces (as in the empty space 2.K.1.5.0.a);
  * replaces 2.K.31.4.0.a.b (same id, label and orbit code) by the lift of
    31.6.a.a, computed by the Saito-Kurokawa generator, in smf_newforms and
    smf_hecke_nf, and sets the cusp_P sums of its space (and, if they include
    the old lift, the traces in smf_hecke_newspace_traces) accordingly;
  * sets is_cuspidal = True on the remaining SK rows of levels 1 and 2.

Before writing anything it checks that each row to delete lifts a form of root
number +1, is the only row of its space and that space has dimension 0; that
2.K.31.4.0.a.b lifts 31.6.a.b and that 31.6.a.a is the only form the generator
lifts there, of the dimension the space records; and that the stored sums are
the old row's values.  Deletions come last, after all updates.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_sk_invalid_lifts [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_sk_invalid_lifts [--datafile PREFIX] --apply
"""
import argparse
import json
from collections import Counter

from sage.all import prime_range
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import as_int, write_update_file
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import (
    SK_query, SK_all_forms, SK_all_evs, HECKE_TYPES)

SPACE_SUMS = ['cusp_P_' + h for h in HECKE_TYPES]
IDENTITY = {'id', 'label', 'hecke_orbit_code', 'space_label', 'level', 'weight', 'family', 'degree',
            'char_orbit_index', 'char_orbit_label', 'hecke_orbit'}
SPACE_31, ROW_31, OLD_31, NEW_31 = '2.K.31.4.0.a', '2.K.31.4.0.a.b', '31.6.a.b', '31.6.a.a'


def classical_label(r):
    return '.'.join(r['related_objects'][0].split('/')[-4:])


def jsonable(x):
    if isinstance(x, (list, tuple)):
        return [jsonable(y) for y in x]
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    if x is None or isinstance(x, (bool, str)):
        return x
    try:
        return int(x)
    except (TypeError, ValueError):
        return str(x)


def main(apply=False, datafile=None):
    base = datafile or "smf_sk_invalid"
    nf_cols = db.smf_newforms.search_cols
    hnf_cols = db.smf_hecke_nf.search_cols

    # 1. the level-1 rows of k odd
    P1 = list(db.smf_newforms.search({'aut_rep_type': 'P', 'level': 1}, ['id', 'label', 'space_label', 'weight', 'related_objects']))
    doomed = [r for r in P1 if r['weight'][0] % 2 == 1]
    cl = {f['label']: f for f in db.mf_newforms.search({'label': {'$in': [classical_label(r) for r in doomed]}},
                                                      ['label', 'level', 'weight', 'fricke_eigenval'])}
    for r in doomed:
        f = cl[classical_label(r)]
        k = r['weight'][0]
        assert (f['level'], int(f['weight']), r['weight'][1]) == (1, 2*k - 2, 0), r['label']
        assert (-1)**(k-1) * f['fricke_eigenval'] == 1, "%s lifts a form of root number -1" % r['label']
    sp_labels = sorted(set(r['space_label'] for r in doomed))
    assert len(sp_labels) == len(doomed), "two rows to delete share a space"
    spaces = {s['label']: s for s in db.smf_newspaces.search({'label': {'$in': sp_labels}},
                                                            ['label', 'num_forms', 'total_dim', 'cusp_P_dim', 'new_cusp_P_dim', 'hecke_orbit_code'] + SPACE_SUMS)}
    sp_rows = []
    for sp in sp_labels:
        s = spaces[sp]
        assert (s['total_dim'], s['cusp_P_dim'], s['new_cusp_P_dim'], s['num_forms']) == (0, 0, 0, 1), sp
        assert db.smf_newforms.count({'space_label': sp}) == 1, sp
        assert all(as_int(t['trace_an']) == 0 for t in db.smf_hecke_newspace_traces.search({'hecke_orbit_code': s['hecke_orbit_code']}, ['trace_an'])), sp
        out = {'label': sp, 'num_forms': 0}
        for c in SPACE_SUMS:
            out[c] = None if s.get(c) is None else [0] * len(s[c])
        sp_rows.append(out)
    doomed_labels = [r['label'] for r in doomed]
    doomed_full = list(db.smf_newforms.search({'label': {'$in': doomed_labels}}, nf_cols + ['id']))
    doomed_hnf = list(db.smf_hecke_nf.search({'label': {'$in': doomed_labels}}, hnf_cols + ['id']))
    print("%d level-1 rows of odd k to delete (k = %s), %d smf_hecke_nf rows; %d spaces set to empty"
          % (len(doomed), sorted(set(r['weight'][0] for r in doomed)), len(doomed_hnf), len(sp_rows)))

    # 2. 2.K.31.4.0.a.b
    old = db.smf_newforms.lookup(ROW_31, nf_cols + ['id'])
    assert classical_label(old) == OLD_31, old['related_objects']
    f_old = db.mf_newforms.lookup(OLD_31, ['fricke_eigenval'])
    assert (-1)**(4-1) * f_old['fricke_eigenval'] == 1, "%s lifts a form of root number -1" % ROW_31
    new_forms = SK_all_forms(4, SK_query(4, 31), True)
    assert len(new_forms) == 1 and classical_label(new_forms[0]) == NEW_31, [x['related_objects'] for x in new_forms]
    new = new_forms[0]
    s31 = db.smf_newspaces.lookup(SPACE_31, ['label', 'new_cusp_P_dim', 'cusp_P_dim', 'hecke_orbit_code'] + SPACE_SUMS)
    assert s31['new_cusp_P_dim'] == s31['cusp_P_dim'] == new['dim'], (s31['new_cusp_P_dim'], new['dim'])
    assert [r['label'] for r in db.smf_newforms.search({'space_label': SPACE_31, 'aut_rep_type': 'P'}, ['label'])] == [ROW_31]
    nf31 = {'id': old['id']}
    for c, v in new.items():
        if c in nf_cols and c not in IDENTITY:
            nf31[c] = v
    nf31['trace_display'] = new['trace_lambda_p'][:4]
    changed = sorted(c for c in nf31 if c != 'id' and jsonable(nf31[c]) != jsonable(old.get(c)))
    print("%s: lift of %s (dim %d) replaced by the lift of %s (dim %d); columns changed: %s"
          % (ROW_31, OLD_31, old['dim'], NEW_31, new['dim'], changed))
    # \N in the update file keeps the stored value, so a new NULL cannot overwrite an old value
    lost = sorted(c for c in nf31 if nf31[c] is None and old.get(c) is not None)
    assert not lost, "%s: new value NULL but stored value present in %s" % (ROW_31, lost)
    kept = sorted(c for c in nf_cols if c not in nf31 and c not in IDENTITY)
    print("    columns kept as stored:", kept)

    old_h = db.smf_hecke_nf.lookup(ROW_31, hnf_cols + ['id'])
    new_evs = SK_all_evs(4, SK_query(4, 31))
    assert len(new_evs) == 1, len(new_evs)
    hnf31 = {'id': old_h['id']}
    for c, v in new_evs[0].items():
        if c in hnf_cols and c not in IDENTITY:
            hnf31[c] = v
    print("    smf_hecke_nf columns changed:", sorted(c for c in hnf31 if c != 'id' and jsonable(hnf31[c]) != jsonable(old_h.get(c))))

    # its space: the P sums are those of the one P orbit
    s31_row = {'label': SPACE_31, 'num_forms': None}
    for h in HECKE_TYPES:
        c = 'cusp_P_' + h
        if s31.get(c) is None:
            s31_row[c] = None
            continue
        stored = [as_int(x) for x in s31[c]]
        assert stored == [as_int(x) for x in old['trace_' + h][:len(stored)]], "%s: stored %s is not the old row's" % (SPACE_31, c)
        s31_row[c] = new['trace_' + h][:len(stored)]
    sp_rows.append(s31_row)
    # the newspace traces of T(p): do they include the old lift?
    G = [r for r in db.smf_newforms.search({'space_label': SPACE_31}, ['label', 'aut_rep_type', 'trace_lambda_p']) if r['aut_rep_type'] != 'P']
    tr_rows = []
    nst = {t['n']: t for t in db.smf_hecke_newspace_traces.search({'hecke_orbit_code': s31['hecke_orbit_code']}, ['id', 'n', 'trace_an'])}
    primes = prime_range(200)
    def at(rows_, n):
        return sum(as_int(r['trace_lambda_p'][primes.index(n)]) for r in rows_)
    with_old = all(as_int(t['trace_an']) == at(G, n) + as_int(old['trace_lambda_p'][primes.index(n)]) for n, t in nst.items() if n in primes)
    without = all(as_int(t['trace_an']) == at(G, n) for n, t in nst.items() if n in primes)
    if not nst:
        print("    smf_hecke_newspace_traces has no rows for %s" % SPACE_31)
    else:
        print("    smf_hecke_newspace_traces of %s (%d rows): %s" % (SPACE_31, len(nst),
              "include the old lift, corrected" if with_old else ("do not include the P orbit, not touched" if without else "neither G nor G + old P: NOT TOUCHED, check by hand")))
    if nst and with_old:
        assert all(n in primes for n in nst), "composite n in smf_hecke_newspace_traces"
        for n, t in nst.items():
            tr_rows.append({'id': t['id'], 'trace_an': at(G, n) + as_int(new['trace_lambda_p'][primes.index(n)])})

    # 3. is_cuspidal
    cusp_rows = [{'id': r['id'], 'is_cuspidal': True}
                 for r in db.smf_newforms.search({'aut_rep_type': 'P', 'level': {'$lte': 2}, 'is_cuspidal': False}, ['id', 'label'])
                 if r['label'] not in doomed_labels]
    print("%d SK rows of level 1 or 2 get is_cuspidal = True" % len(cusp_rows))

    # files
    backup = base + "_deleted.json"
    with open(backup, 'w') as F:
        json.dump({'smf_newforms': jsonable(doomed_full), 'smf_hecke_nf': jsonable(doomed_hnf),
                   'replaced_smf_newforms': jsonable(old), 'replaced_smf_hecke_nf': jsonable(old_h)}, F)
    print("wrote %s (the rows deleted or replaced)" % backup)
    todo = []
    nf_upd = [dict(nf31)] + cusp_rows
    nf_upd_cols = sorted(set(c for x in nf_upd for c in x if c != 'id'))
    for x in nf_upd:
        for c in nf_upd_cols:
            x.setdefault(c, None)   # \N keeps the stored value
    fname = base + "_newforms.txt"
    write_update_file(fname, db.smf_newforms, 'id', nf_upd, nf_upd_cols)
    todo.append((db.smf_newforms, fname, 'id'))
    fname = base + "_newspaces.txt"
    write_update_file(fname, db.smf_newspaces, 'label', sp_rows, ['num_forms'] + SPACE_SUMS)
    todo.append((db.smf_newspaces, fname, 'label'))
    fname = base + "_hecke_nf.txt"
    write_update_file(fname, db.smf_hecke_nf, 'id', [hnf31], sorted(c for c in hnf31 if c != 'id'))
    todo.append((db.smf_hecke_nf, fname, 'id'))
    if tr_rows:
        fname = base + "_newspace_traces.txt"
        write_update_file(fname, db.smf_hecke_newspace_traces, 'id', tr_rows, ['trace_an'])
        todo.append((db.smf_hecke_newspace_traces, fname, 'id'))
    if not apply:
        print("dry run; pass --apply to write")
        return
    # Each update_from_file builds a replacement table; undo with <table>.reload_revert().
    for table, fname, key in todo:
        table.update_from_file(fname, label_col=key)
    # Deletions last; the deleted rows are in the JSON backup.
    db.smf_hecke_nf.delete({'label': {'$in': doomed_labels}})
    db.smf_newforms.delete({'label': {'$in': doomed_labels}})
    print("deleted %d smf_newforms and %d smf_hecke_nf rows" % (len(doomed_labels), len(doomed_hnf)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update and backup files (default: ./smf_sk_invalid)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
