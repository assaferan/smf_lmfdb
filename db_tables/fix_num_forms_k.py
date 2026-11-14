"""
Correct num_forms of the paramodular newspaces (family K, level > 2).

num_forms_paramodular already counts the Saito-Kurokawa lifts (one per
classical newform orbit in S_{2k-2}^new(Gamma_0(N)) with Fricke sign (-1)^k,
Hecke_Eigenvalues_PM_SK_new_num_forms), and the populate script then added
num_classical_minus_cusp_dim(2k-2, N), the same count, again.  So every space
with lifts has num_forms = #G + 2 #P (950 spaces), and 2.K.19.3.2.a, where j > 0
and there are no lifts, counts the classical forms of weight 4 anyway.
num_forms is set to the number of stored newform orbits, #G + #P.

Before writing anything it checks that every stored num_forms is #G + #P plus
what the populate script added, num_classical_minus_cusp_dim(2k-2, N), and that
for j = 0 the number of stored P orbits is that count.

Run after fix_sk_invalid_lifts (which replaces the P row of 2.K.31.4.0.a).

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_num_forms_k [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_num_forms_k [--datafile PREFIX] --apply
"""
import argparse
from collections import Counter

from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file


def classical_minus(w, N, counts):
    """num_classical_minus_cusp_dim(w, N) of smf_newspaces_populate, from prefetched counts."""
    return counts.get((w, N, (-1)**(w//2 - 1)), 0)


def main(apply=False, datafile=None):
    spaces = list(db.smf_newspaces.search({'family': 'K', 'level': {'$gt': 2}, 'num_forms': {'$exists': True}},
                                          ['label', 'level', 'weight', 'num_forms']))
    orbits = {}
    for f in db.smf_newforms.search({'family': 'K', 'level': {'$gt': 2}}, ['space_label', 'aut_rep_type']):
        orbits.setdefault(f['space_label'], Counter())[f['aut_rep_type']] += 1
    ws = sorted(set(2*s['weight'][0] - 2 for s in spaces))
    counts = Counter()
    for f in db.mf_newforms.search({'weight': {'$in': ws}, 'char_order': 1,
                                    'level': {'$in': sorted(set(s['level'] for s in spaces))}},
                                   ['weight', 'level', 'fricke_eigenval']):
        counts[(f['weight'], f['level'], f['fricke_eigenval'])] += 1
    print("%d paramodular newspaces of level > 2 with num_forms" % len(spaces))
    rows, stats = [], Counter()
    for s in spaces:
        N, (k, j) = s['level'], s['weight']
        o = orbits.get(s['label'], Counter())
        assert set(o) <= {'G', 'P'}, "%s: %s" % (s['label'], dict(o))
        extra = classical_minus(2*k - 2, N, counts)
        if j == 0:
            assert o['P'] == extra, "%s: %d P orbits, %d classical forms with a lift" % (s['label'], o['P'], extra)
        else:
            assert o['P'] == 0, s['label']
        new = o['G'] + o['P']
        assert s['num_forms'] == new + extra, "%s: num_forms %d, #G + #P = %d, populate added %d" % (s['label'], s['num_forms'], new, extra)
        if new != s['num_forms']:
            rows.append({'label': s['label'], 'num_forms': new})
            stats['change'] += 1
        else:
            stats['unchanged'] += 1
    print(dict(stats))
    for x in rows[:3]:
        print("   ", x)
    if not rows:
        print("nothing to correct")
        return
    fname = (datafile or "smf_num_forms_k_fix") + "_newspaces.txt"
    write_update_file(fname, db.smf_newspaces, 'label', rows, ['num_forms'])
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newspaces.reload_revert().
    db.smf_newspaces.update_from_file(fname, label_col='label')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update file (default: ./smf_num_forms_k_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
