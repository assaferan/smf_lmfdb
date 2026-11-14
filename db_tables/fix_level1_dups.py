"""
Remove the 67 extra copies of the Siegel Eisenstein series E_4 in each of the
level-1 spaces 2.K.1.4.0.a, 2.S.1.4.0.a, 2.P.1.4.0.a, and empty the level-1
spaces of odd k and odd j that claim Klingen-Eisenstein series.

* create_entries (smf_newforms_populate) made the level-1 and level-2 forms of
  weight (k, j) once for every triple (k, j, N) outside the paramodular range
  (N square, or N >= 1000), not once per weight.  Only (4, 0) has such triples
  with forms to copy: the 68 squarefree N in [1000, 1111] (N_bound(6) = 1111),
  so each of the three level-1 spaces of weight (4, 0), of dimension 1, has 68
  identical rows of type F, labelled .a to .cp.  They have no smf_hecke_nf,
  smf_hecke_traces or smf_hecke_newspace_traces rows, and the space sums were
  computed separately (they are right).  The row of hecke_orbit 1 (.a) is kept.
* KE_num_forms and KE_space_traces took the level-1 newforms of weight j + k
  without asking j to be even; for j odd there are no Klingen-Eisenstein series
  (M_{k,j}(Sp(4,Z)) = 0), so the 69 spaces of k, j odd whose weight j + k has
  newforms (in each of the three families) say num_forms = 1 and store the
  Klingen-Eisenstein traces of those newforms in eis_Q_lambda_p*.  Their
  dimensions are all 0 and they have no rows.  num_forms becomes 0 and the
  eis_Q_lambda_p* columns 0.

Before writing anything it checks that each deleted row is the kept row up to
id, label, hecke_orbit and hecke_orbit_code, that nothing else refers to it, and
that each space it empties has dimension 0, no rows, and exactly the old
KE_num_forms count and KE_space_traces sums.  The deleted rows are saved to
PREFIX_deleted.json first, and the deletion runs last.  Each part does nothing
once it has been applied, so the script can be run again.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_level1_dups [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_level1_dups [--datafile PREFIX] --apply
"""
import argparse
import json
from collections import Counter

from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.db_tables.fix_sk_invalid_lifts import jsonable
from smf_lmfdb.Hecke_Eigenvalues.Klingen_Yoshida import HECKE_TYPES, KE_lifts, KE_query

E4_SPACES = ['2.%s.1.4.0.a' % F for F in 'KSP']
DIFFERING = {'id', 'label', 'hecke_orbit', 'hecke_orbit_code'}
DIM_COLS = ['total_dim', 'cusp_dim', 'eis_dim', 'eis_F_dim', 'eis_Q_dim', 'cusp_P_dim', 'cusp_Y_dim', 'cusp_G_dim']
KE_COLS = ['eis_Q_' + ht for ht in HECKE_TYPES]


def e4_copies(cols):
    """The rows to delete: all but .a in the three spaces of E_4."""
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
        assert len(extra) in (0, 67), (sp, len(extra))
        for r in extra:
            diff = sorted(c for c in cols if c not in DIFFERING and jsonable(r.get(c)) != jsonable(keep.get(c)))
            assert not diff, "%s differs from %s in %s" % (r['label'], keep['label'], diff)
        doomed += extra
    codes = [r['hecke_orbit_code'] for r in doomed]
    for t in ['smf_hecke_nf', 'smf_hecke_traces', 'smf_hecke_newspace_traces']:
        n = getattr(db, t).count({'hecke_orbit_code': {'$in': codes}}) if codes else 0
        assert n == 0, "%s has %d rows for the copies" % (t, n)
    return doomed


def old_ke_sums(k, j, lengths):
    """What KE_space_traces stored before it asked j to be even: the sums over the newforms of weight j + k."""
    sums = {c: [0] * n for c, n in lengths.items()}
    for _, lift in KE_lifts(k, j):
        tr = lift.traces(num_coeffs=0)
        for c in sums:
            v = tr['trace_' + c[len('eis_Q_'):]]
            assert len(v) >= len(sums[c]), (k, j, c)
            sums[c] = [a + b for a, b in zip(sums[c], v)]
    return sums


def odd_spaces():
    """The update rows of the level-1 spaces of odd k and odd j that still claim Klingen-Eisenstein series."""
    rows_per_space = Counter(f['space_label'] for f in db.smf_newforms.search({'level': 1}, ['space_label']))
    out = []
    for s in db.smf_newspaces.search({'level': 1}, ['label', 'weight', 'num_forms'] + DIM_COLS + KE_COLS):
        k, j = s['weight']
        if k % 2 == 0 or j % 2 == 0:
            continue
        nonzero = [c for c in KE_COLS if s.get(c) and any(x not in (0, None) for x in s[c])]
        if not s.get('num_forms') and not nonzero:
            continue
        assert all(not s[c] for c in DIM_COLS), "%s: dimensions %s" % (s['label'], {c: s[c] for c in DIM_COLS})
        assert rows_per_space.get(s['label'], 0) == 0, s['label']
        old_count = db.mf_newforms.count(KE_query(k, j)) if k >= 4 else 0
        assert s.get('num_forms') in (0, old_count), "%s: num_forms %s, the old KE_num_forms %d" % (s['label'], s.get('num_forms'), old_count)
        row = {'label': s['label'], 'num_forms': 0 if s.get('num_forms') else None}
        if nonzero:
            old = old_ke_sums(k, j, {c: len(s[c]) for c in KE_COLS})
            for c in KE_COLS:
                assert s[c] == old[c], "%s: %s is not the old KE_space_traces sum" % (s['label'], c)
                row[c] = [0] * len(s[c])
        out.append(row)
    return out


def main(apply=False, datafile=None):
    base = datafile or "smf_level1_dups"
    doomed = e4_copies(db.smf_newforms.search_cols)
    print("%d copies of E_4 to delete (keeping %s)" % (len(doomed), ", ".join(sp + '.a' for sp in E4_SPACES)))
    sp_out = odd_spaces()
    print("%d level-1 spaces of odd k and odd j to empty: num_forms in %d, eis_Q_lambda_p* in %d"
          % (len(sp_out), sum(r['num_forms'] is not None for r in sp_out), sum('eis_Q_lambda_p' in r for r in sp_out)))
    for r in sp_out:
        for c in KE_COLS:
            r.setdefault(c, None)    # \N: keep the stored value
    if sp_out:
        fname = base + "_newspaces.txt"
        write_update_file(fname, db.smf_newspaces, 'label', sp_out, ['num_forms'] + KE_COLS)
    if doomed:
        backup = base + "_deleted.json"
        with open(backup, 'w') as F:
            json.dump({'smf_newforms': jsonable(doomed)}, F)
        print("wrote %s (the rows to delete)" % backup)
    if not sp_out and not doomed:
        print("nothing to correct")
        return
    if not apply:
        print("dry run; pass --apply to write")
        return
    if sp_out:
        # update_from_file builds a replacement table; undo with db.smf_newspaces.reload_revert().
        db.smf_newspaces.update_from_file(fname, label_col='label')
    if doomed:
        # Deletion last; the deleted rows are in the JSON backup.
        db.smf_newforms.delete({'label': {'$in': [r['label'] for r in doomed]}})
        print("deleted %d smf_newforms rows" % len(doomed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update and backup files (default: ./smf_level1_dups)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
