"""
Store T(p^2), not T(p,2), for the paramodular non-lifts (family K, aut_rep_type G).

omf5_data gives, for p <= 7, the eigenvalues of T(p,2) and their traces (its
README), and the populate scripts stored them unchanged as lambda_p_square and
trace_lambda_p_square: the conversion to T(p^2) in Hecke_Eigenvalues_paramodular
never ran (it fails at nth_prime(0)).  Everything else in these columns, the
Saito-Kurokawa lifts of the same spaces included, is T(p^2), and so is the
label of the site.  (Read as T(p,2) the stored values satisfy the Ramanujan
bound at every good prime; read as T(p^2) they do not.)

    T(p^2) = T(p)^2 - p^(mu-1) - (T(p,2) + p^2 + 1) p^(mu-2),    mu = j + 2k - 3.

smf_newforms: trace_lambda_p_square of the G rows.  The trace of T(p^2) needs the
sum of lambda_p^2 over the orbit; it comes from the eigenvalue coordinates when
the orbit has them, and otherwise from the Hecke matrices of omf5_data (see
omf5_squares), except for 53 orbits at levels with a single good prime below
10, where it is set to NULL.  For those orbits without coordinates `traces` has only a_1 and
a_p; a_{p^2} = T(p^2) + p^(mu-1) is filled in for p <= 7.
smf_hecke_nf: lambda_p_square, converted at each embedding.
smf_newspaces: cusp_G_lambda_p_square, the sum over the G rows.

Before writing anything it checks that the stored values are those of omf5_data
(T(p,2), and the same T(p) coordinates), that the sums of squares from the
coordinates and from the Hecke matrices agree wherever there are both, that the
traces of the new hnf rows are the new trace_lambda_p_square, that the stored
a_{p^2} in `traces` (computed from the coordinates by fix_g_dirichlet_traces) are
T(p^2) + p^(mu-1), and that the stored space sums are the sums of the G rows.

It needs the omf5_data submodule (git submodule update --init
Hecke_Eigenvalues/paramodular/omf5_data, from smf_lmfdb).

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_g_t_p2 [--datafile PREFIX] [--procs N]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_g_t_p2 [--datafile PREFIX] [--procs N] --apply
"""
import argparse
from collections import Counter, defaultdict
from multiprocessing import Pool

from sage.all import ZZ
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.db_tables.nf_elt import get_nf_basis, nf_elts_to_lists, nf_lists_to_elements
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import parse_omf5
from smf_lmfdb.Hecke_Eigenvalues.paramodular.omf5_squares import (
    PRIMES, SQUARE_PRIMES, OrbitNotPlaced, has_coordinates, mu, sums_of_squares_from_coordinates,
    sums_of_squares_from_matrices, T_p_square, T_p_square_traces)


def nullify(v):
    return [None if x == 'NULL' else x for x in v]


def orbit_key(al, dim, trace_lambda_p, trace_T_p_2):
    # T(p,2) too: two G orbits of level 748 have the same traces of T(p) for all p < 100
    def t(v):
        return tuple(None if x in ('NULL', None) else int(x) for x in v)
    return (al, int(dim), t(trace_lambda_p), t(trace_T_p_2))


def level_data(args):
    """The G orbits of omf5_data of weight (k,j) and level N, with their sums of squares."""
    k, j, N = args
    out = {}
    for f in parse_omf5(k, j, N):
        if f['aut_rep_type'] != 'G':
            continue
        sq, how = {}, None
        if has_coordinates(f):
            sq, how = sums_of_squares_from_coordinates(f, N), 'coordinates'
            if (k, j) == (3, 0) and sq:
                sq_small = {p: v for p, v in sq.items() if p in SQUARE_PRIMES}
                try:
                    from_matrices = sums_of_squares_from_matrices(f, N)
                except OrbitNotPlaced:
                    # too few good primes below 10 to tell this orbit's block apart
                    how = 'coordinates (matrices cannot tell the orbit apart)'
                else:
                    assert from_matrices == sq_small, \
                        "level %d: the coordinates and the Hecke matrices give different sums of squares" % N
                    how = 'coordinates, = matrices'
        elif any(t not in ('NULL', None) for t in f['trace_lambda_p_square']):
            try:
                sq, how = sums_of_squares_from_matrices(f, N), 'matrices'
            except OrbitNotPlaced:
                # a single good prime below 10: T(p^2) stays unknown
                how = 'none (one good prime below 10), NULL'
        key = orbit_key(f['atkin_lehner_string'], f['dim'], f['trace_lambda_p'], f['trace_lambda_p_square'])
        assert key not in out, "level %d: two G orbits %s" % (N, key)
        out[key] = {'orbit': f, 'sq': sq, 'how': how}
    return (k, j, N), out


def main(apply=False, datafile=None, procs=8):
    cols = ['id', 'label', 'level', 'weight', 'dim', 'space_label', 'atkin_lehner_string',
            'trace_lambda_p', 'trace_lambda_p_square', 'traces']
    rows = list(db.smf_newforms.search({'family': 'K', 'aut_rep_type': 'G', 'trace_lambda_p_square': {'$exists': True}}, cols))
    hnf_cols = ['id', 'label', 'field_poly', 'hecke_ring_power_basis', 'hecke_ring_numerators', 'hecke_ring_denominators',
                'hecke_ring_inverse_numerators', 'hecke_ring_inverse_denominators', 'lambda_p', 'lambda_p_square']
    hnf = {}
    for h in db.smf_hecke_nf.search({'label': {'$in': [r['label'] for r in rows]}}, hnf_cols):
        assert h['label'] not in hnf, h['label']
        hnf[h['label']] = h
    triples = sorted(set((r['weight'][0], r['weight'][1], r['level']) for r in rows))
    print("%d G rows (%d with hnf rows) in %d spaces" % (len(rows), len(hnf), len(triples)))
    with Pool(procs) as pool:
        omf5 = dict(pool.imap_unordered(level_data, triples))

    nf_out, hnf_out, stats = [], [], Counter()
    new_sq = {}
    for r in rows:
        (k, j), N, d = r['weight'], r['level'], r['dim']
        m = mu(k, j)
        key = orbit_key(r['atkin_lehner_string'], d, r['trace_lambda_p'], r['trace_lambda_p_square'])
        assert key in omf5[(k, j, N)], "%s: no G orbit of omf5_data with its traces of T(p) and T(p,2)" % r['label']
        o = omf5[(k, j, N)][key]
        f, sq = o['orbit'], o['sq']
        stats['sum of squares from ' + o['how']] += 1
        assert nullify(r['trace_lambda_p_square']) == nullify(f['trace_lambda_p_square']), \
            "%s: stored trace_lambda_p_square is not T(p,2) of omf5_data" % r['label']
        new = nullify(T_p_square_traces(f, N, k, j, sq))
        if sq == {} and any(t is not None for t in new):
            raise AssertionError(r['label'])
        new_sq[r['label']] = new
        row = {'id': r['id'], 'trace_lambda_p_square': new}
        # a_{p^2} = T(p^2) + p^(mu-1)
        traces = list(r['traces'])
        for i, t in enumerate(new):
            p = PRIMES[i]
            if t is None or p**2 > len(traces):
                continue
            a = t + d*p**(m-1)
            if traces[p**2 - 1] is None:
                assert not has_coordinates(f), "%s: a_%d missing" % (r['label'], p**2)
                traces[p**2 - 1] = a
                stats['a_{p^2} filled in'] += 1
            else:
                assert traces[p**2 - 1] == a, "%s: stored a_%d = %s, T(p^2) gives %s" % (r['label'], p**2, traces[p**2 - 1], a)
                stats['a_{p^2} checked'] += 1
        if traces != r['traces']:
            row['traces'] = traces
        nf_out.append(row)

        h = hnf.get(r['label'])
        if h is None:
            assert not has_coordinates(f) or d == 1, r['label']
            continue
        assert has_coordinates(f), r['label']
        assert nullify(h['lambda_p_square']) == nullify(f['lambda_p_square'][:len(h['lambda_p_square'])]), \
            "%s: hnf lambda_p_square is not T(p,2) of omf5_data" % r['label']
        assert nullify(h['lambda_p'][:len(f['lambda_p'])]) == nullify(f['lambda_p'][:len(h['lambda_p'])]), r['label']
        basis, inv_basis = get_nf_basis(h)
        ap = nf_lists_to_elements(h['lambda_p'], basis)
        bp2 = nf_lists_to_elements(h['lambda_p_square'], basis)
        elts = ['NULL' if b == 'NULL' else T_p_square(ap[i], b, PRIMES[i], k, j) for i, b in enumerate(bp2)]
        lists = nf_elts_to_lists(elts, inv_basis)
        for x in lists:
            assert x == 'NULL' or all(ZZ(c) == c for c in x), "%s: T(p^2) not integral in the Hecke ring basis" % r['label']
        lists = ['NULL' if x == 'NULL' else [int(c) for c in x] for x in lists]
        assert [None if e == 'NULL' else int(e.trace()) if hasattr(e, 'trace') else int(e)*d for e in elts] == new[:len(elts)], \
            "%s: the hnf traces are not the new trace_lambda_p_square" % r['label']
        hnf_out.append({'id': h['id'], 'lambda_p_square': lists})

    # the space sums
    by_space = defaultdict(list)
    for r in rows:
        by_space[r['space_label']].append(r)
    sp_out = []
    for s in db.smf_newspaces.search({'label': {'$in': list(by_space)}}, ['label', 'cusp_G_lambda_p_square']):
        rs = by_space[s['label']]
        n = len(s['cusp_G_lambda_p_square'])

        def total(vals):
            return [None if any(v[i] is None for v in vals) else sum(v[i] for v in vals) for i in range(n)]
        assert nullify(s['cusp_G_lambda_p_square']) == total([nullify(r['trace_lambda_p_square'])[:n] for r in rs]), \
            "%s: cusp_G_lambda_p_square is not the sum of the G rows" % s['label']
        sp_out.append({'label': s['label'], 'cusp_G_lambda_p_square': total([new_sq[r['label']][:n] for r in rs])})
    assert len(sp_out) == len(by_space)

    print(dict(stats))
    for x in nf_out[:2] + hnf_out[:1] + sp_out[:1]:
        print("   ", x)
    prefix = datafile or "smf_g_t_p2_fix"
    files = [(prefix + "_newforms.txt", db.smf_newforms, 'id', nf_out, ['trace_lambda_p_square', 'traces']),
             (prefix + "_hecke_nf.txt", db.smf_hecke_nf, 'id', hnf_out, ['lambda_p_square']),
             (prefix + "_newspaces.txt", db.smf_newspaces, 'label', sp_out, ['cusp_G_lambda_p_square'])]
    for fname, table, key, out, cs in files:
        for x in out:
            for c in cs:
                x.setdefault(c, None)    # \N: keep the stored value
        write_update_file(fname, table, key, out, cs)
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with reload_revert() on each table.
    for fname, table, key, out, cs in files:
        table.update_from_file(fname, label_col=key)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update files (default: ./smf_g_t_p2_fix)")
    parser.add_argument("--procs", type=int, default=8, help="processes for the Hecke matrices")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile, procs=args.procs)
