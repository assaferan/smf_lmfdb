"""
Fill smf_hecke_traces and smf_hecke_newspace_traces from the stored Dirichlet coefficients, as
the trace search pages read them: trace_an is the trace Tr(a_n) of the coefficient a_n of the
spin Dirichlet series (the n can be primes, prime powers or all n), as for classical forms
(mf_hecke_traces, mf_hecke_newspace_traces).

  smf_hecke_traces:           for each newform and each n <= len(traces) with a stored value,
                              trace_an = traces[n-1] (smf_newforms.traces is Tr(a_n)).
  smf_hecke_newspace_traces:  for each newspace whose newforms are all stored (num_forms of
                              them) and each n at which all of them have a value, the sum of
                              multiplicity * traces[n-1] over its newforms.

smf_hecke_traces was empty; smf_hecke_newspace_traces held, for 2268 spaces, the traces of
lambda_p at the primes over the forms of general type only (cusp_G_lambda_p, or cusp_Y_lambda_p),
which is neither the trace of a_n nor a trace over the whole newspace.  Both tables are
rewritten (reload, with a backup).

Dry run: writes PREFIX_{traces,newspace_traces}.txt and compares with the stored rows.
--apply: reloads both tables from the files and reads every row back.

    sage -python -m smf_lmfdb.db_tables.fill_trace_tables [--datafile PREFIX] [--apply]
"""
import argparse
from collections import Counter, defaultdict

from lmfdb import db

COLS = ['hecke_orbit_code', 'n', 'trace_an']


def newform_rows():
    rows, stats = [], Counter()
    forms = list(db.smf_newforms.search({}, ['label', 'space_label', 'hecke_orbit_code', 'traces', 'multiplicity']))
    for f in forms:
        t = f.get('traces')
        if not t:
            stats['newforms without traces'] += 1
            continue
        for n, x in enumerate(t, 1):
            if x is None:
                stats['a_n not stored (skipped)'] += 1
                continue
            rows.append({'hecke_orbit_code': f['hecke_orbit_code'], 'n': n, 'trace_an': int(x)})
    stats['newforms'] = len(forms)
    return forms, rows, stats


def space_rows(forms):
    rows, stats = [], Counter()
    by_space = defaultdict(list)
    for f in forms:
        by_space[f['space_label']].append(f)
    spaces = list(db.smf_newspaces.search({'num_forms': {'$gt': 0}}, ['label', 'hecke_orbit_code', 'num_forms']))
    for s in spaces:
        fs = by_space.get(s['label'], [])
        if len(fs) != s['num_forms']:
            stats['spaces with num_forms != stored newforms (skipped)'] += 1
            continue
        if any(not f.get('traces') for f in fs):
            stats['spaces with a newform without traces (skipped)'] += 1
            continue
        nmax = min(len(f['traces']) for f in fs)
        added = 0
        for n in range(1, nmax + 1):
            vals = [f['traces'][n - 1] for f in fs]
            if any(v is None for v in vals):
                continue
            rows.append({'hecke_orbit_code': s['hecke_orbit_code'], 'n': n,
                         'trace_an': int(sum((f.get('multiplicity') or 1) * v for f, v in zip(fs, vals)))})
            added += 1
        stats['spaces with traces'] += 1 if added else 0
    stats['newspaces with forms'] = len(spaces)
    return rows, stats


def write(fname, rows):
    with open(fname, 'w') as F:
        F.write("|".join(COLS) + "\nbigint|integer|numeric\n\n")
        for r in rows:
            F.write("%d|%d|%d\n" % (r['hecke_orbit_code'], r['n'], r['trace_an']))
    print("wrote %s (%d rows)" % (fname, len(rows)))


def compare(table, rows):
    old = {(r['hecke_orbit_code'], r['n']): r['trace_an'] for r in table.search({}, COLS)}
    new = {(r['hecke_orbit_code'], r['n']): r['trace_an'] for r in rows}
    same = sum(1 for k in new if k in old and old[k] == new[k])
    changed = sum(1 for k in new if k in old and old[k] != new[k])
    print("%s: %d stored rows, %d new rows: %d same, %d changed, %d added, %d dropped; codes stored %d, new %d"
          % (table.search_table, len(old), len(new), same, changed, sum(1 for k in new if k not in old),
             sum(1 for k in old if k not in new), len({k[0] for k in old}), len({k[0] for k in new})))


def main(apply=False, datafile='smf_traces'):
    forms, nf_rows, s1 = newform_rows()
    sp_rows, s2 = space_rows(forms)
    print(dict(s1)); print(dict(s2))
    compare(db.smf_hecke_traces, nf_rows)
    compare(db.smf_hecke_newspace_traces, sp_rows)
    f1, f2 = datafile + '_traces.txt', datafile + '_newspace_traces.txt'
    write(f1, nf_rows)
    write(f2, sp_rows)
    if not apply:
        print("dry run; pass --apply to reload the two tables")
        return
    for table, fname, rows in [(db.smf_hecke_traces, f1, nf_rows), (db.smf_hecke_newspace_traces, f2, sp_rows)]:
        table.reload(fname)
        got = {(r['hecke_orbit_code'], r['n']): r['trace_an'] for r in table.search({}, COLS)}
        want = {(r['hecke_orbit_code'], r['n']): r['trace_an'] for r in rows}
        bad = [k for k in want if got.get(k) != want[k]] + [k for k in got if k not in want]
        print("%s: read back %d rows, %d mismatches" % (table.search_table, len(got), len(bad)))
        assert not bad, bad[:10]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datafile", default="smf_traces")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
