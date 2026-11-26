"""
Fill smf_qexp_short, the short fully expanded Fourier expansions shown on the newform pages,
from smf_qexp_coeffs (one coefficient per class) and smf_qexp_reduction (the classes).

For each orbit with stored coefficients, of level N and weight (k, 0), and each index
T = [n, r, N m] (the matrix [[n, r/2], [r/2, N m]], i.e. the term q_1^n q_12^r q_2^(N m)) with
1 <= n, m <= B and 4 N n m > r^2, it writes a(T) = (s_T s_rep)^k a(rep): the class of T and the
signs s come from qexp_from_eigenforms.reduce_signed (s = 0: a(T) = 0 in odd weight), and rep is
the representative stored in smf_qexp_reduction.  B is the largest bound <= BMAX for which every
such T has a stored class.

Columns: hecke_orbit_code, qf = [n, r, N m], trace = n + N m, index = 0 (scalar weight), coeff =
the coordinates of a(T), as in smf_qexp_coeffs.

Dry run: writes PREFIX.txt and prints the counts.  --apply: reloads smf_qexp_short from it and
reads it back.  The column coeff must be bigint[] (as smf_qexp_coeffs.coeff).

    sage -python -m smf_lmfdb.db_tables.fill_qexp_short [--bmax B] [--datafile PREFIX] [--procs P] [--apply]
"""
import argparse
import os
import sys
from math import isqrt
from multiprocessing import Pool

from lmfdb import db

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qexp_from_eigenforms import reduce_signed  # noqa: E402


def box(N, B):
    return [(n, r, m) for n in range(1, B + 1) for m in range(1, B + 1)
            for r in range(-isqrt(4 * N * n * m), isqrt(4 * N * n * m) + 1) if 4 * N * n * m > r * r]


def classes_of_level(args):
    """{(n, r, m): (key, s)} for the box of level N, and {key: s_rep}."""
    N, B, reps = args
    out = {}
    for t in box(N, B):
        qf, tag, s = reduce_signed(t, N)
        out[t] = ((tuple(int(x) for x in qf), (int(tag[0]), int(tag[1]))), int(s))
    srep = {}
    for key, rep in reps.items():
        _, _, s = reduce_signed((rep[0], rep[1], rep[2] // N), N)
        srep[key] = int(s) if s else 1
    return N, out, srep


def main(bmax=7, datafile='smf_qexp_short', procs=8, apply=False):
    codes = sorted({r['hecke_orbit_code'] for r in db.smf_qexp_coeffs.search({}, ['hecke_orbit_code'])})
    forms = {f['hecke_orbit_code']: f for f in db.smf_newforms.search({'hecke_orbit_code': {'$in': codes}},
                                                                       ['hecke_orbit_code', 'label', 'level', 'weight', 'family'])}
    print("%d orbits with coefficients, %d with a newform row" % (len(codes), len(forms)))
    levels = sorted({f['level'] for f in forms.values()})
    reps = {N: {(tuple(r['qf_legendre']), tuple(r['qf_tag'])): r['qf_orbit_rep']
                for r in db.smf_qexp_reduction.search({'level': N, 'family': 'K'}, ['qf_legendre', 'qf_tag', 'qf_orbit_rep'])}
            for N in levels}
    with Pool(procs) as pool:
        cls = {N: (out, srep) for N, out, srep in pool.imap_unordered(classes_of_level, [(N, bmax, reps[N]) for N in levels])}
    rows = []
    for code in codes:
        f = forms.get(code)
        if f is None:
            print("  orbit %d: no newform row, skipped" % code)
            continue
        N, (k, j) = f['level'], f['weight']
        assert f['family'] == 'K' and j == 0, f
        co = {(tuple(c['qf_legendre']), tuple(c['qf_tag'])): c['coeff']
              for c in db.smf_qexp_coeffs.search({'hecke_orbit_code': code}, ['qf_legendre', 'qf_tag', 'coeff'])}
        out, srep = cls[N]
        B = bmax
        while B > 0 and any(out[t][0] not in co for t in box(N, B)):
            B -= 1
        if B < bmax:
            print("  %s: complete only up to n, m <= %d" % (f['label'], B))
        dim = len(next(iter(co.values())))
        for (n, r, m) in box(N, B):
            key, s = out[(n, r, m)]
            if s == 0:
                if k % 2:
                    v = [0] * dim
                else:
                    v = [srep[key] ** k * x for x in co[key]]
            else:
                v = [(s * srep[key]) ** k * x for x in co[key]]
            rows.append((code, (n, r, N * m), n + N * m, 0, v))
    print("%d rows for %d orbits (bound n, m <= %d)" % (len(rows), len({r[0] for r in rows}), bmax))
    fname = datafile + '.txt'
    with open(fname, 'w') as F:
        F.write("hecke_orbit_code|qf|trace|index|coeff\nbigint|integer[]|smallint|smallint|bigint[]\n\n")
        for code, qf, tr, idx, v in rows:
            F.write("%d|{%d,%d,%d}|%d|%d|{%s}\n" % (code, qf[0], qf[1], qf[2], tr, idx, ",".join(str(int(x)) for x in v)))
    print("wrote", fname)
    if not apply:
        print("dry run; pass --apply to reload smf_qexp_short")
        return
    table = db.smf_qexp_short
    assert table.col_type['coeff'] == 'bigint[]', "smf_qexp_short.coeff must be bigint[] first"
    table.reload(fname)
    got = table.count()
    print("smf_qexp_short: %d rows (expected %d)" % (got, len(rows)))
    assert got == len(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bmax", type=int, default=7)
    parser.add_argument("--datafile", default="smf_qexp_short")
    parser.add_argument("--procs", type=int, default=8)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    main(args.bmax, args.datafile, args.procs, args.apply)
