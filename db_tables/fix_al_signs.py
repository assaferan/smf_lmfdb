"""
Correct the Atkin-Lehner data of the paramodular newforms computed from
quinary orthogonal modular forms (family K, aut_rep_type G).

These were loaded from omf5_data, whose atkin_lehner_eigenvals encoded the
spinor-character conductor d rather than the paramodular Atkin-Lehner signs.
The stored atkin_lehner_string is correct for weight (3,0); for the single
forms of weight (3,2), level 19 and weight (4,0), level 31 the string was wrong
too (d = 1 at a prime level p_0 gives sign -1, see DPRT Theorem 10.1).

smf_newforms: atkin_lehner_string, atkin_lehner_eigenvals and fricke_eigenval
are rewritten from the (corrected) string.
smf_newspaces: ALdims_G is recomputed from the corrected strings of the G
newforms in each space, and ALdims = ALdims_G + ALdims_P.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_al_signs [--datafile FILE]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_al_signs [--datafile FILE] --apply
"""
import argparse
from collections import defaultdict

from sage.all import divisors, is_squarefree, prod
from lmfdb import db
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import al_str_to_eigenvals, al_str_to_num

STRING_FIXES = {'2.K.19.3.2.a.a': '-', '2.K.31.4.0.a.a': '-'}


def corrected(rec):
    al_str = STRING_FIXES.get(rec['label'], rec['atkin_lehner_string'])
    evs = al_str_to_eigenvals(al_str, rec['level'])
    return {'atkin_lehner_string': al_str, 'atkin_lehner_eigenvals': evs,
            'fricke_eigenval': int(prod(e for _, e in evs))}


def newform_changes(forms):
    changes = {}
    for rec in forms:
        new = corrected(rec)
        if any(rec.get(k) != v for k, v in new.items()):
            changes[rec['label']] = new
    return changes


def newspace_changes(forms):
    # ALdims_G is indexed by the squarefree divisors of N in increasing order;
    # a form with Atkin-Lehner string s counts towards the divisor al_str_to_num(s, N).
    al_dims = defaultdict(dict)
    for rec in forms:
        N = rec['level']
        d = al_str_to_num(corrected(rec)['atkin_lehner_string'], N)
        al_dims[rec['space_label']][d] = al_dims[rec['space_label']].get(d, 0) + rec['dim']
    changes = {}
    cols = ['label', 'level', 'ALdims_G', 'ALdims_P', 'ALdims']
    for sp in db.smf_newspaces.search({'family': 'K', 'ALdims_G': {'$exists': True}}, cols):
        divs = [d for d in divisors(sp['level']) if is_squarefree(d)]
        dims_G = [al_dims[sp['label']].get(d, 0) for d in divs]
        dims = [g + p for g, p in zip(dims_G, sp['ALdims_P'])]
        if dims_G != sp['ALdims_G'] or dims != sp['ALdims']:
            changes[sp['label']] = {'ALdims_G': dims_G, 'ALdims': dims}
    return changes


def write_update_file(fname, table, changes, cols, fmt):
    types = [table.col_type[c] for c in ['label'] + cols]
    with open(fname, 'w') as F:
        F.write("|".join(['label'] + cols) + "\n" + "|".join(types) + "\n\n")
        for label, new in changes.items():
            F.write("|".join([label] + [fmt(new[c]) for c in cols]) + "\n")
    print("wrote", fname)


def pg_array(x):
    if isinstance(x, (list, tuple)):
        return "{" + ",".join(pg_array(y) for y in x) + "}"
    return str(x)


def main(apply=False, datafile=None):
    cols = ['label', 'level', 'dim', 'space_label', 'atkin_lehner_string',
            'atkin_lehner_eigenvals', 'fricke_eigenval']
    # level 1 has no bad primes, so nothing to correct
    forms = [rec for rec in db.smf_newforms.search({'family': 'K', 'aut_rep_type': 'G'}, cols)
             if rec['level'] != 1]
    nf_changes = newform_changes(forms)
    sp_changes = newspace_changes(forms)
    print("%d newforms need correcting" % len(nf_changes))
    for label in list(nf_changes)[:5]:
        print("  ", label, nf_changes[label])
    print("%d newspaces need correcting" % len(sp_changes))
    for label in list(sp_changes)[:5]:
        print("  ", label, sp_changes[label])

    # One bulk update per table: a per-row update() pays lmfdb's lock check
    # (including a fresh devmirror connection) each time, about 4s per row.
    # update_from_file builds a replacement table, so it can be undone with
    # reload_revert().
    base = datafile or "smf_al_fix"
    todo = []
    if nf_changes:
        fname = base + "_newforms.txt"
        write_update_file(fname, db.smf_newforms, nf_changes,
                          ['atkin_lehner_string', 'atkin_lehner_eigenvals', 'fricke_eigenval'],
                          pg_array)
        todo.append((db.smf_newforms, fname))
    if sp_changes:
        fname = base + "_newspaces.txt"
        write_update_file(fname, db.smf_newspaces, sp_changes, ['ALdims_G', 'ALdims'], pg_array)
        todo.append((db.smf_newspaces, fname))
    if not apply:
        print("dry run; pass --apply to write")
        return
    for table, fname in todo:
        table.update_from_file(fname, label_col='label')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update files (default: ./smf_al_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
