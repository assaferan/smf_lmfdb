"""
Correct smf_hecke_newspace_traces.

smf_hecke_newspace_traces_populate stored, for each space and each prime p up
to its bound, trace_an = aps[i-1] where p is the i-th prime (0-based), i.e. the
value of the previous prime (and for p = 2 the value of the last prime). This
script recomputes every row from the same source column the populate script
uses (cusp_G_lambda_p for weight (3,0) at non-square level, cusp_Y_lambda_p
otherwise), now with aps[i], and updates the rows that differ.

It refuses to run if the set of (space, n) rows in the table differs from the
set the populate script would produce: it only corrects values.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_newspace_traces [--datafile FILE]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_newspace_traces [--datafile FILE] --apply
"""
import argparse

from sage.all import is_square, nth_prime, prime_range
from lmfdb import db
from smf_lmfdb.db_tables.common_populate import MAX_P


def as_int(x):
    return None if x is None or x == 'NULL' else int(x)


def expected_rows(space):
    """{n: trace_an} as get_space_entries should have produced it."""
    k, j = space['weight']
    if (k == 3) and (j == 0) and (not is_square(space['level'])):
        aps = space['cusp_G_lambda_p']
    else:
        aps = space['cusp_Y_lambda_p']
    aps = aps or []
    max_p = min(MAX_P, nth_prime(len(aps))) if aps else 0
    return {int(p): as_int(aps[i]) for i, p in enumerate(prime_range(1, max_p + 1))}


def find_changes(table):
    rows = list(table.search({}, ['id', 'hecke_orbit_code', 'n', 'trace_an']))
    codes = {r['hecke_orbit_code'] for r in rows}
    spaces = {s['hecke_orbit_code']: s for s in db.smf_newspaces.search(
        {'hecke_orbit_code': {'$in': list(codes)}},
        ['label', 'hecke_orbit_code', 'level', 'weight', 'cusp_G_lambda_p', 'cusp_Y_lambda_p'])}
    print("%d rows for %d spaces (%d found in smf_newspaces)" % (len(rows), len(codes), len(spaces)))
    assert set(spaces) == codes, "rows for unknown spaces: %s" % sorted(codes - set(spaces))[:5]

    expected = {c: expected_rows(s) for c, s in spaces.items()}
    have = {(r['hecke_orbit_code'], int(r['n'])) for r in rows}
    want = {(c, n) for c, e in expected.items() for n in e}
    assert have == want, "row sets differ: %d missing, %d extra" % (len(want - have), len(have - want))

    # search results leave out NULL fields, hence .get
    wrong = [r for r in rows if as_int(r.get('trace_an')) != expected[r['hecke_orbit_code']][int(r['n'])]]
    changes = [{'id': r['id'], 'trace_an': expected[r['hecke_orbit_code']][int(r['n'])]} for r in wrong]
    bad_spaces = sorted({spaces[r['hecke_orbit_code']]['label'] for r in wrong})
    print("%d rows to correct (%d of them to NULL), in %d spaces"
          % (len(changes), sum(c['trace_an'] is None for c in changes), len(bad_spaces)))
    for lab in bad_spaces[:5]:
        print("   ", lab)
    return changes


def main(apply=False, datafile=None):
    table = db.smf_hecke_newspace_traces
    changes = find_changes(table)
    if not changes:
        print("nothing to correct")
        return
    values = [c for c in changes if c['trace_an'] is not None]
    nulls = [c['id'] for c in changes if c['trace_an'] is None]
    fname = datafile or "smf_newspace_traces_fix.txt"
    if values:
        with open(fname, 'w') as F:
            F.write("id|trace_an\n" + table.col_type['id'] + "|" + table.col_type['trace_an'] + "\n\n")
            for c in values:
                F.write("%s|%s\n" % (c['id'], c['trace_an']))
        print("wrote", fname)
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table (undo with reload_revert()), but it
    # merges with COALESCE(new, current), so it cannot set a value to NULL; those
    # rows are cleared with a single update() instead.
    if values:
        table.update_from_file(fname, label_col='id')
    if nulls:
        table.update({'id': {'$in': nulls}}, {'trace_an': None})
    print("checking the table after the update:")
    left = find_changes(table)
    assert not left, "%d rows are still wrong after the update" % len(left)
    print("all rows correct")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="update file to write (default: ./smf_newspace_traces_fix.txt)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
