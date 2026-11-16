"""
Fill in T(p^2) for the 53 paramodular non-lifts left NULL by fix_g_t_p2.

They have weight (3,0), no eigenvalue coordinates in omf5_data, and a level N
with a single prime p < 10 not dividing N, where the Hecke matrices of omf5_data
cannot tell the orbits apart.  The sum over the orbit of lambda_p^2 at that p now
comes from runs of omf5 (Hecke_Eigenvalues/paramodular/omf5_squares_runs.json).

    smf_newforms:  trace_lambda_p_square at p, and a_{p^2} = T(p^2) + dim p^(mu-1) in traces;
    smf_newspaces: cusp_G_lambda_p_square at p, the sum over the G rows of the space.

Before writing it checks that each orbit of the file is a single stored row and a
single orbit of omf5_data, that the trace of T(p,2) of the file is that of
omf5_data, that the values to fill in are NULL and nothing else changes, that
the rest of trace_lambda_p_square is what fix_g_t_p2 gives, and that the space
sums were NULL at p only because of these orbits.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_g_t_p2_omf5 [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_g_t_p2_omf5 [--datafile PREFIX] --apply
"""
import argparse
from collections import defaultdict

from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import parse_omf5
from smf_lmfdb.Hecke_Eigenvalues.paramodular.omf5_squares import (
    PRIMES, _omf5_runs, has_coordinates, mu, sums_of_squares_from_omf5_runs, T_p_square_traces)

NF_COLS = ['id', 'label', 'level', 'weight', 'dim', 'space_label', 'atkin_lehner_string',
           'trace_lambda_p', 'trace_lambda_p_square', 'traces']


def nullify(v):
    return [None if x in ('NULL', None) else int(x) for x in v]


def main(apply=False, datafile=None):
    runs = _omf5_runs()
    nf_out, labels = [], {}
    for o in runs:
        N, p, d = o['level'], o['p'], o['dim']
        i = PRIMES.index(p)
        assert [q for q in [2, 3, 5, 7] if N % q] == [p], N
        rows = [r for r in db.smf_newforms.search(
            {'family': 'K', 'aut_rep_type': 'G', 'level': N, 'weight': [3, 0], 'dim': d,
             'atkin_lehner_string': o['atkin_lehner_string']}, NF_COLS)
            if nullify(r['trace_lambda_p'][:len(o['trace_lambda_p'])]) == o['trace_lambda_p']]
        assert len(rows) == 1, "level %d: %d stored rows for an orbit of the file" % (N, len(rows))
        r = rows[0]
        fs = [f for f in parse_omf5(3, 0, N) if f['aut_rep_type'] == 'G' and f['dim'] == d
              and f['atkin_lehner_string'] == o['atkin_lehner_string']
              and nullify(f['trace_lambda_p']) == nullify(r['trace_lambda_p'])]
        assert len(fs) == 1, "%s: %d orbits of omf5_data" % (r['label'], len(fs))
        f = fs[0]
        assert not has_coordinates(f), r['label']
        sq = sums_of_squares_from_omf5_runs(f, N)    # also checks the trace of T(p,2)
        new = nullify(T_p_square_traces(f, N, 3, 0, sq))
        old = nullify(r['trace_lambda_p_square'])
        assert old[i] is None and new[i] is not None, r['label']
        assert [x for k, x in enumerate(old) if k != i] == [x for k, x in enumerate(new) if k != i], \
            "%s: trace_lambda_p_square differs away from p = %d" % (r['label'], p)
        traces = list(r['traces'])
        a = new[i] + d * p**(mu(3, 0) - 1)
        assert len(traces) >= p**2 and traces[p**2 - 1] is None, r['label']
        traces[p**2 - 1] = a
        nf_out.append({'id': r['id'], 'trace_lambda_p_square': new, 'traces': traces})
        labels[r['label']] = (r['space_label'], i, new[i])
        print("%s  p = %d  sum lambda_p^2 = %d  T(p^2) trace = %d  a_%d = %d"
              % (r['label'], p, sq[p], new[i], p**2, a))
    assert len(labels) == len(runs)

    # the space sums
    by_space = defaultdict(dict)
    for l, (s, i, t) in labels.items():
        by_space[s][l] = (i, t)
    sp_out = []
    for s in db.smf_newspaces.search({'label': {'$in': list(by_space)}}, ['label', 'cusp_G_lambda_p_square']):
        rs = list(db.smf_newforms.search({'space_label': s['label'], 'aut_rep_type': 'G'}, ['label', 'trace_lambda_p_square']))
        cur = nullify(s['cusp_G_lambda_p_square'])
        new = list(cur)
        for i in {i for i, t in by_space[s['label']].values()}:
            assert cur[i] is None, s['label']
            vals = [by_space[s['label']][r['label']][1] if r['label'] in by_space[s['label']]
                    else nullify(r['trace_lambda_p_square'])[i] for r in rs]
            assert None not in vals, "%s: another G row is NULL at %d" % (s['label'], PRIMES[i])
            new[i] = sum(vals)
        sp_out.append({'label': s['label'], 'cusp_G_lambda_p_square': new})
        print("%s  cusp_G_lambda_p_square %s -> %s" % (s['label'], cur, new))
    assert len(sp_out) == len(by_space)
    print("%d newforms, %d newspaces" % (len(nf_out), len(sp_out)))

    prefix = datafile or "smf_g_t_p2_omf5"
    files = [(prefix + "_newforms.txt", db.smf_newforms, 'id', nf_out, ['trace_lambda_p_square', 'traces']),
             (prefix + "_newspaces.txt", db.smf_newspaces, 'label', sp_out, ['cusp_G_lambda_p_square'])]
    for fname, table, key, out, cs in files:
        write_update_file(fname, table, key, out, cs)
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with reload_revert() on each table.
    for fname, table, key, out, cs in files:
        table.update_from_file(fname, label_col=key)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--datafile')
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
