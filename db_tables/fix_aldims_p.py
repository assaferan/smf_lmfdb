"""
Correct ALdims_P of the paramodular newspaces, and store ALdims_G, ALdims_P and
ALdims = ALdims_G + ALdims_P in the order the space pages read.

ALdims_P[i] is the dimension of the Saito-Kurokawa lifts (new and old) whose
Atkin-Lehner signs are -1 exactly at the primes of the i-th squarefree divisor
d of N, Saito_Kurokawa_lift_dim(k, N, al=d).  It is read off
mf_newspaces.ALdims, which is ordered by sign vectors (+ before -, primes
increasing), but the dimension formulas indexed it by the squarefree divisors
in increasing order.  The two orders agree at prime levels and differ once N
has two prime factors: for N = 22 the divisors 1, 2, 11, 22 are the sign
vectors ++, -+, +-, --, at indices 0, 2, 1, 3.  DimensionFormulas now uses
cmf_ALdims_index.  The populate script also filled ALdims_P for j > 0, where
there are no lifts (2.K.19.3.2.a); those are set to 0.

All three columns were stored in the order of the squarefree divisors, but the
Atkin-Lehner table of the space page (web_space.ALdim_table) reads them in the
sign-vector order of mf_newspaces.ALdims, so with two or more primes it shows
dimensions under the wrong signs (2.K.161.3.0.a: the form with w_7 = -1 and
the three with w_23 = -1 are swapped).  They are rewritten in sign-vector order,
which is what the populate script now produces.

Before writing anything it checks that every stored ALdims_P is what the old
indexing gives, that the stored ALdims is ALdims_G + ALdims_P, that the stored
ALdims_G is, in divisor order, the dimension of the stored G newforms of each
Atkin-Lehner string (and after reordering, in sign order), and that in the
spaces without old lifts the new ALdims_P is, class by class, the dimension of
the classical newforms of that class with root number -1 (from mf_newforms, not
from the dimension formulas).

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_aldims_p [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_aldims_p [--datafile PREFIX] --apply
"""
import argparse
from collections import Counter

from sage.all import divisors, gcd, is_odd, is_squarefree, prime_divisors, prod
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import CMF_SPACES, Saito_Kurokawa_lift_dim, cmf_ALdims_index


def old_lift_dim(k, N, al):
    """Saito_Kurokawa_lift_dim(k, N, al) with the old indexing of mf_newspaces.ALdims."""
    w = 2*k - 2
    total = 0
    for M in divisors(N):
        f = CMF_SPACES[(w, M)]
        if f['dim'] == 0:
            continue
        a = gcd(al, M) if al != 0 else 0
        if a == 0:
            total += f['plus_dim'] if is_odd(w//2) else f['dim'] - f['plus_dim']
            continue
        divs = [d for d in divisors(M) if is_squarefree(d)]
        if is_odd(w//2 - len(prime_divisors(a))):
            total += f['ALdims'][divs.index(a)]
    return total


def al_string(d, N):
    return ''.join('-' if d % p == 0 else '+' for p in prime_divisors(N))


def sign_order(v, N):
    """v indexed by the squarefree divisors of N in increasing order, reindexed by sign vectors."""
    divs = [d for d in divisors(N) if is_squarefree(d)]
    out = [None] * len(divs)
    for d, x in zip(divs, v):
        out[cmf_ALdims_index(d, N)] = x
    return out


def main(apply=False, datafile=None):
    spaces = list(db.smf_newspaces.search({'family': 'K', 'ALdims_P': {'$exists': True}},
                                          ['label', 'level', 'weight', 'ALdims_G', 'ALdims_P', 'ALdims', 'old_cusp_P_dim']))
    ws = sorted(set(2*s['weight'][0] - 2 for s in spaces))
    levels = sorted(set(M for s in spaces for M in divisors(s['level'])))
    for f in db.mf_newspaces.search({'weight': {'$in': ws}, 'char_orbit_index': 1, 'level': {'$in': levels}},
                                    ['weight', 'level', 'dim', 'plus_dim', 'ALdims']):
        CMF_SPACES[(f['weight'], f['level'])] = f
    classical = {}
    for f in db.mf_newforms.search({'weight': {'$in': ws}, 'char_order': 1, 'level': {'$in': sorted(set(s['level'] for s in spaces))}},
                                   ['weight', 'level', 'dim', 'fricke_eigenval', 'atkin_lehner_eigenvals']):
        classical.setdefault((f['weight'], f['level']), []).append(f)
    G = {}
    for f in db.smf_newforms.search({'family': 'K', 'aut_rep_type': 'G', 'space_label': {'$in': [s['label'] for s in spaces]}},
                                    ['space_label', 'atkin_lehner_string', 'dim']):
        c = G.setdefault(f['space_label'], Counter())
        c[f['atkin_lehner_string']] += f['dim']
    missing = [(w, M) for w in ws for M in levels if (w, M) not in CMF_SPACES]
    assert not missing, "no mf_newspaces row for %s" % missing[:5]
    print("%d paramodular newspaces with ALdims_P" % len(spaces))
    rows, stats = [], Counter()
    for s in spaces:
        N, (k, j) = s['level'], s['weight']
        divs = [d for d in divisors(N) if is_squarefree(d)]
        stored = [int(x) for x in s['ALdims_P']]
        assert len(stored) == len(divs), s['label']
        assert [int(x) for x in s['ALdims']] == [int(g) + p for g, p in zip(s['ALdims_G'], stored)], s['label']
        stored_G = [int(x) for x in s['ALdims_G']]
        assert stored_G == [G.get(s['label'], {}).get(al_string(d, N), 0) for d in divs], \
            "%s: stored ALdims_G is not the G newforms in divisor order" % s['label']
        # the populate script computed ALdims_P without looking at j
        assert stored == [old_lift_dim(k, N, d) for d in divs], "%s: stored ALdims_P is not the old formula" % s['label']
        if j != 0:
            new = [0] * len(divs)     # no Saito-Kurokawa lifts for j > 0
            stats['j > 0, set to 0'] += 1
        else:
            new = [int(Saito_Kurokawa_lift_dim(k, N, al=d)) for d in divs]
        assert all(x >= 0 for x in new), s['label']
        if j == 0 and s.get('old_cusp_P_dim') == 0:
            # only new lifts: each class is the classical newforms of that class with root number -1
            cls = Counter()
            for f in classical.get((2*k - 2, N), []):
                if f['fricke_eigenval'] == (-1)**k:
                    cls[prod(int(p) for p, e in f['atkin_lehner_eigenvals'] if e == -1)] += f['dim']
            assert new == [cls[d] for d in divs], "%s: %s, classical newforms give %s" % (s['label'], new, [cls[d] for d in divs])
            stats['checked against the classical newforms'] += 1
        new_G, new_P = sign_order(stored_G, N), sign_order(new, N)
        signs = sign_order([al_string(d, N) for d in divs], N)
        assert signs == sorted(signs, key=lambda x: x.replace('+', '0').replace('-', '1')), s['label']
        assert new_G == [G.get(s['label'], {}).get(x, 0) for x in signs], s['label']
        new_all = [g + p for g, p in zip(new_G, new_P)]
        if new_G != stored_G:
            stats['reordered'] += 1
        if new != stored:
            stats['ALdims_P corrected'] += 1
        if (new_G, new_P, new_all) != (stored_G, stored, [int(x) for x in s['ALdims']]):
            rows.append({'label': s['label'], 'ALdims_G': new_G, 'ALdims_P': new_P, 'ALdims': new_all})
            stats['change'] += 1
        else:
            stats['unchanged'] += 1
    print(dict(stats))
    for x in rows[:3]:
        print("   ", x)
    if not rows:
        print("nothing to correct")
        return
    fname = (datafile or "smf_aldims_p_fix") + "_newspaces.txt"
    write_update_file(fname, db.smf_newspaces, 'label', rows, ['ALdims_G', 'ALdims_P', 'ALdims'])
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newspaces.reload_revert().
    db.smf_newspaces.update_from_file(fname, label_col='label')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update file (default: ./smf_aldims_p_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
