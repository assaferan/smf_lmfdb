"""
Remove the 67 extra copies of the Siegel Eisenstein series E_4 in each of the
level-1 spaces 2.K.1.4.0.a, 2.S.1.4.0.a, 2.P.1.4.0.a, and set num_forms = 0 in
the level-1 spaces of odd k and odd j that claim one form.

* create_entries (smf_newforms_populate) made the level-1 and level-2 forms of
  weight (k, j) once for every triple (k, j, N) outside the paramodular range
  (N square, or N >= 1000), not once per weight.  Only (4, 0) has such triples
  with forms to copy: the 68 squarefree N in [1000, 1111] (N_bound(6) = 1111),
  so each of the three level-1 spaces of weight (4, 0), of dimension 1, has 68
  identical rows of type F, labelled .a to .cp.  They have no smf_hecke_nf,
  smf_hecke_traces or smf_hecke_newspace_traces rows, and the space sums were
  computed separately (they are right).  The row of hecke_orbit 1 (.a) is kept.
* KE_num_forms counted the level-1 newforms of weight j + k without asking j to
  be even; for j odd there are no Klingen-Eisenstein series (M_{k,j}(Sp(4,Z)) = 0),
  so the 69 spaces of k, j odd whose weight j + k has newforms say num_forms = 1
  (in each of the three families).  Their dimensions are all 0 and they have no rows.

Before writing anything it checks that each deleted row is the kept row up to
id, label, hecke_orbit and hecke_orbit_code, that nothing else refers to it, and
that each space set to num_forms = 0 has dimension 0, no rows, and the old count
of KE_num_forms.  The deleted rows are saved to PREFIX_deleted.json first, and the
deletion runs last.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_level1_dups [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_level1_dups [--datafile PREFIX] --apply
"""
import argparse
import json

from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable

E4_SPACES = ['2.%s.1.4.0.a' % F for F in 'KSP']
DIFFERING = {'id', 'label', 'hecke_orbit', 'hecke_orbit_code'}
DIM_COLS = ['total_dim', 'cusp_dim', 'eis_dim', 'eis_F_dim', 'eis_Q_dim', 'cusp_P_dim', 'cusp_Y_dim', 'cusp_G_dim']


def main(apply=False, datafile=None):
    base = datafile or "smf_level1_dups"
    cols = db.smf_newforms.search_cols

    # 1. the copies of E_4
    doomed = []
    for sp in E4_SPACES:
        s = db.smf_newspaces.lookup(sp, ['total_dim', 'eis_F_dim', 'num_forms'])
        assert (s['total_dim'], s['eis_F_dim'], s['num_forms']) == (1, 1, 1), (sp, s)
        rows = list(db.smf_newforms.search({'space_label': sp}, cols))
        keep = [r for r in rows if r['hecke_orbit'] == 1]
        assert len(keep) == 1 and keep[0]['label'] == sp + '.a', sp
        keep = keep[0]
        assert keep['aut_rep_type'] == 'F' and keep['trace_lambda_p'][0] == 45, (sp, keep['trace_lambda_p'][:1])  # 1 + 2^2 + 2^3 + 2^5
        extra = [r for r in rows if r['hecke_orbit'] != 1]
        assert len(extra) == 67, (sp, len(extra))
        for r in extra:
            diff = sorted(c for c in cols if c not in DIFFERING and jsonable(r.get(c)) != jsonable(keep.get(c)))
            assert not diff, "%s differs from %s in %s" % (r['label'], keep['label'], diff)
        doomed += extra
    codes = [r['hecke_orbit_code'] for r in doomed]
    for t in ['smf_hecke_nf', 'smf_hecke_traces', 'smf_hecke_newspace_traces']:
        n = getattr(db, t).count({'hecke_orbit_code': {'$in': codes}})
        assert n == 0, "%s has %d rows for the copies" % (t, n)
    print("%d copies of E_4 to delete (%d per space, keeping %s)" % (len(doomed), len(doomed) // 3, ", ".join(sp + '.a' for sp in E4_SPACES)))

    # 2. num_forms of the empty spaces of odd k and odd j
    ke_count = {}
    for f in db.mf_newforms.search({'level': 1, 'char_order': 1}, ['weight']):
        ke_count[f['weight']] = ke_count.get(f['weight'], 0) + 1
    rows_per_space = {}
    for f in db.smf_newforms.search({'level': 1}, ['space_label']):
        rows_per_space[f['space_label']] = rows_per_space.get(f['space_label'], 0) + 1
    sp_out = []
    for s in db.smf_newspaces.search({'level': 1}, ['label', 'weight', 'num_forms'] + DIM_COLS):
        k, j = s['weight']
        if k % 2 == 0 or j % 2 == 0 or not s.get('num_forms'):
            continue
        old_ke = ke_count.get(j + k, 0) if k >= 4 else 0
        assert s['num_forms'] == old_ke, "%s: num_forms %s, KE_num_forms gave %d" % (s['label'], s['num_forms'], old_ke)
        assert all(not s[c] for c in DIM_COLS), "%s: dimensions %s" % (s['label'], {c: s[c] for c in DIM_COLS})
        assert rows_per_space.get(s['label'], 0) == 0, s['label']
        sp_out.append({'label': s['label'], 'num_forms': 0})
    print("%d level-1 spaces of odd k and odd j: num_forms -> 0" % len(sp_out))

    fname = base + "_newspaces.txt"
    write_update_file(fname, db.smf_newspaces, 'label', sp_out, ['num_forms'])
    backup = base + "_deleted.json"
    with open(backup, 'w') as F:
        json.dump({'smf_newforms': jsonable(doomed)}, F)
    print("wrote %s (the rows to delete)" % backup)
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newspaces.reload_revert().
    db.smf_newspaces.update_from_file(fname, label_col='label')
    # Deletion last; the deleted rows are in the JSON backup.
    db.smf_newforms.delete({'label': {'$in': [r['label'] for r in doomed]}})
    print("deleted %d smf_newforms rows" % len(doomed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update and backup files (default: ./smf_level1_dups)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
