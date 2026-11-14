"""
Correct the Atkin-Lehner signs of the paramodular Saito-Kurokawa lifts
(family K, aut_rep_type P, level > 1).

The Saito-Kurokawa lift F to K(N) of a newform f in S_{2k-2}(Gamma_0(N)) has,
at every p | N, the Atkin-Lehner eigenvalue of f: w_p(F) = eps_p(f) = w_p(f)
([DPRT] Dummigan, Pacetti, Rama, Tornaria, Quinary forms and paramodular forms,
Theorem 10.1(2), via [RS07, Proposition 5.5.8(i)]).  So the Fricke sign of F is
that of f, which is (-1)^k for f to have a lift.  The stored rows have the
populate default instead: + at every p and Fricke sign +1.

smf_newforms: atkin_lehner_eigenvals, atkin_lehner_string and fricke_eigenval
are set to those of the classical form.

Before writing anything it checks that every row stores the default, that the
classical form's primes are those of N and its Fricke sign is (-1)^k, and that
in each space whose lifts are all stored (new_cusp_P_dim = their total
dimension) the corrected signs give, class by class, the dimensions
Saito_Kurokawa_new_lift_dim(k, N, al=d) of the dimension formulas, which come
from mf_newspaces.ALdims and not from the signs.  (This needs the corrected
indexing of mf_newspaces.ALdims in DimensionFormulas; with the old one, the
classes of two primes are swapped in 5 of these spaces.)

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_sk_al_signs [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_sk_al_signs [--datafile PREFIX] --apply
"""
import argparse
from collections import Counter, defaultdict

from sage.all import divisors, is_squarefree, prime_divisors, prod
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import al_str_to_num
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import CMF_SPACES, Saito_Kurokawa_new_lift_dim


def classical_label(r):
    return '.'.join(r['related_objects'][0].split('/')[-4:])


def main(apply=False, datafile=None):
    cols = ['id', 'label', 'level', 'weight', 'dim', 'space_label', 'related_objects',
            'atkin_lehner_string', 'atkin_lehner_eigenvals', 'fricke_eigenval']
    rows = list(db.smf_newforms.search({'family': 'K', 'aut_rep_type': 'P', 'level': {'$gt': 1}}, cols))
    cl = {f['label']: f for f in db.mf_newforms.search({'label': {'$in': [classical_label(r) for r in rows]}},
                                                      ['label', 'level', 'weight', 'atkin_lehner_eigenvals', 'fricke_eigenval'])}
    print("%d paramodular SK rows of level > 1" % len(rows))
    out, stats = [], Counter()
    for r in rows:
        N, k = r['level'], r['weight'][0]
        ps = [int(p) for p in prime_divisors(N)]
        f = cl[classical_label(r)]
        assert (f['level'], int(f['weight'])) == (N, 2*k - 2), r['label']
        # the stored default
        assert r['atkin_lehner_string'] == '+' * len(ps) and r['fricke_eigenval'] == 1 and \
            sorted(map(list, r['atkin_lehner_eigenvals'])) == [[p, 1] for p in ps], \
            "%s: stored signs are not the populate default" % r['label']
        evs = sorted([int(p), int(e)] for p, e in f['atkin_lehner_eigenvals'])
        assert [p for p, _ in evs] == ps, "%s: primes of %s" % (r['label'], f['label'])
        fricke = int(prod(e for _, e in evs))
        assert fricke == f['fricke_eigenval'], r['label']
        if fricke != (-1)**k:
            # not a lift (2.K.31.4.0.a.b before fix_sk_invalid_lifts): left to that script
            print("    skipped %s: %s has Fricke sign %d, so it has no lift" % (r['label'], f['label'], fricke))
            stats['skipped, not a lift'] += 1
            r['skip'] = True
            continue
        new = {'id': r['id'], 'atkin_lehner_eigenvals': evs,
               'atkin_lehner_string': ''.join('+' if e == 1 else '-' for _, e in evs),
               'fricke_eigenval': fricke}
        r.update(new)
        stats['change'] += 1
        out.append(new)

    # the dimension formula, by Atkin-Lehner class, against the corrected signs
    by_space = defaultdict(list)
    for r in rows:
        by_space[r['space_label']].append(r)
    spaces = {s['label']: s for s in db.smf_newspaces.search({'label': {'$in': list(by_space)}},
                                                            ['label', 'level', 'weight', 'new_cusp_P_dim'])}
    ws = sorted(set(2*s['weight'][0] - 2 for s in spaces.values()))
    for f in db.mf_newspaces.search({'weight': {'$in': ws}, 'char_orbit_index': 1, 'level': {'$in': sorted(set(s['level'] for s in spaces.values()))}},
                                    ['weight', 'level', 'dim', 'plus_dim', 'ALdims']):
        CMF_SPACES[(f['weight'], f['level'])] = f
    for sp, rs in by_space.items():
        s = spaces[sp]
        N, k = s['level'], s['weight'][0]
        if any(r.get('skip') for r in rs):
            stats['spaces with a skipped row, not compared'] += 1
            continue
        if s['new_cusp_P_dim'] != sum(r['dim'] for r in rs):
            stats['spaces with lifts missing, not compared'] += 1
            continue
        divs = [d for d in divisors(N) if is_squarefree(d)]
        dims = Counter()
        for r in rs:
            dims[al_str_to_num(r['atkin_lehner_string'], N)] += r['dim']
        want = [Saito_Kurokawa_new_lift_dim(k, N, al=d) for d in divs]
        assert [dims[d] for d in divs] == want, "%s: corrected signs give %s, the dimension formula %s" % (sp, [dims[d] for d in divs], want)
        stats['spaces where the signs give the dimension formula'] += 1
    print(dict(stats))
    for x in out[:3]:
        print("   ", x)

    fname = (datafile or "smf_sk_al_fix") + "_newforms.txt"
    write_update_file(fname, db.smf_newforms, 'id', out, ['atkin_lehner_eigenvals', 'atkin_lehner_string', 'fricke_eigenval'])
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newforms.reload_revert().
    db.smf_newforms.update_from_file(fname, label_col='id')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update file (default: ./smf_sk_al_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
