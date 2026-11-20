"""
Make ALdims of the paramodular newspaces the Atkin-Lehner decomposition of the
new subspace.

ALdims = ALdims_G + ALdims_P, where ALdims_G counts the G newforms of each
class of Atkin-Lehner signs, but ALdims_P was Saito_Kurokawa_lift_dim(k, N, al=d),
which sums the new lifts of every level M | N with signs gcd(d, M).  An old lift
from level M was thus counted in each of the 2^omega(N/M) classes that differ
from it only at the primes of N/M, and the sum of ALdims was neither the new
nor the full cusp dimension (2.K.129.3.0.a: ALdims [0, 12, 4, 4], sum 20, with
new cusp dimension 12 and cusp dimension 16).  An old form has no well defined
signs at the primes of N/M in any case.

ALdims_P becomes Saito_Kurokawa_new_lift_dim(k, N, al=d), which is what the
populate script now uses.  Before writing anything it checks, in every space,
that

- ALdims_G is the dimension of the stored G newforms of each class,
- the new ALdims_P is the dimension of the stored P newforms of each class,
- the new ALdims sums to new_cusp_dim.

The classes are in sign-vector order (+ before -, the first prime the most
significant), as in mf_newspaces.ALdims and fix_aldims_p.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_aldims_new [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_aldims_new [--datafile PREFIX] --apply
"""
import argparse
from collections import Counter, defaultdict

from sage.all import divisors, is_squarefree, prime_divisors
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import Saito_Kurokawa_new_lift_dim, cmf_ALdims_index


def sign_index(s):
    return int(s.replace('+', '0').replace('-', '1'), 2)


def main(apply=False, datafile=None):
    spaces = list(db.smf_newspaces.search({'family': 'K', 'ALdims': {'$exists': True}},
                                          ['label', 'level', 'weight', 'ALdims_G', 'ALdims_P', 'ALdims', 'new_cusp_dim']))
    stored = defaultdict(lambda: defaultdict(Counter))
    for f in db.smf_newforms.search({'family': 'K', 'is_cuspidal': True, 'space_label': {'$in': [s['label'] for s in spaces]}},
                                    ['space_label', 'aut_rep_type', 'atkin_lehner_string', 'dim', 'multiplicity']):
        stored[f['space_label']][f['aut_rep_type']][f['atkin_lehner_string']] += f['dim'] * (f.get('multiplicity') or 1)
    print("%d paramodular newspaces with ALdims" % len(spaces))
    rows, stats = [], Counter()
    for s in spaces:
        N, (k, j) = s['level'], s['weight']
        n = len(prime_divisors(N))
        forms = stored[s['label']]
        assert set(forms) <= {'G', 'P'}, "%s: forms of type %s" % (s['label'], sorted(forms))
        assert all(len(x) == n for t in forms.values() for x in t), "%s: an Atkin-Lehner string of the wrong length" % s['label']
        G = [0] * 2**n
        P = [0] * 2**n
        for x, d in forms['G'].items():
            G[sign_index(x)] += d
        for x, d in forms['P'].items():
            P[sign_index(x)] += d
        assert G == [int(x) for x in s['ALdims_G']], "%s: ALdims_G %s, stored G newforms %s" % (s['label'], s['ALdims_G'], G)
        new_P = [0] * 2**n
        if j == 0:
            for d in divisors(N):
                if is_squarefree(d):
                    new_P[cmf_ALdims_index(d, N)] = int(Saito_Kurokawa_new_lift_dim(k, N, al=d))
        assert new_P == P, "%s: new lifts %s, stored P newforms %s" % (s['label'], new_P, P)
        new_all = [g + p for g, p in zip(G, new_P)]
        assert sum(new_all) == s['new_cusp_dim'], "%s: %s does not sum to %s" % (s['label'], new_all, s['new_cusp_dim'])
        if (new_P, new_all) != ([int(x) for x in s['ALdims_P']], [int(x) for x in s['ALdims']]):
            rows.append({'label': s['label'], 'ALdims_P': new_P, 'ALdims': new_all})
            stats['change'] += 1
        else:
            stats['unchanged'] += 1
    print(dict(stats))
    for x in rows[:3]:
        print("   ", x)
    if not rows:
        print("nothing to correct")
        return
    fname = (datafile or "smf_aldims_new_fix") + "_newspaces.txt"
    write_update_file(fname, db.smf_newspaces, 'label', rows, ['ALdims_P', 'ALdims'])
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newspaces.reload_revert().
    db.smf_newspaces.update_from_file(fname, label_col='label')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update file (default: ./smf_aldims_new_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
