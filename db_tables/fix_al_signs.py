"""
Correct atkin_lehner_eigenvals and fricke_eigenval of the paramodular newforms
computed from quinary orthogonal modular forms (families K, aut_rep_type G).

These were loaded from omf5_data, whose atkin_lehner_eigenvals encoded the
spinor-character conductor d rather than the paramodular Atkin-Lehner signs.
The stored atkin_lehner_string is correct for weight (3,0); for the single
forms of weight (3,2), level 19 and weight (4,0), level 31 the string was wrong
too (d = 1 at a prime level p_0 gives sign -1, see DPRT Theorem 10.1).

Usage (from the directory containing smf_lmfdb, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_al_signs          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_al_signs --apply
"""
import sys
from functools import reduce

from sage.all import prime_divisors
from lmfdb import db

STRING_FIXES = {'2.K.19.3.2.a.a': '-', '2.K.31.4.0.a.a': '-'}


def corrected(rec):
    al_str = STRING_FIXES.get(rec['label'], rec['atkin_lehner_string'])
    ps = prime_divisors(rec['level'])
    assert len(al_str) == len(ps), rec['label']
    evs = [[int(p), -1 if s == '-' else 1] for p, s in zip(ps, al_str)]
    fricke = reduce(lambda x, y: x * y, [e for _, e in evs], 1)
    return {'atkin_lehner_string': al_str, 'atkin_lehner_eigenvals': evs,
            'fricke_eigenval': fricke}


def main(apply=False):
    cols = ['label', 'level', 'atkin_lehner_string', 'atkin_lehner_eigenvals',
            'fricke_eigenval']
    changes = {}
    for rec in db.smf_newforms.search({'family': 'K', 'aut_rep_type': 'G'}, cols):
        if rec['level'] == 1:  # no bad primes, nothing to correct
            continue
        new = corrected(rec)
        if any(rec[k] != v for k, v in new.items()):
            changes[rec['label']] = new
    print("%d newforms need correcting" % len(changes))
    for label in list(changes)[:5]:
        print("  ", label, changes[label])
    if not apply:
        print("dry run; pass --apply to write")
        return
    for i, (label, new) in enumerate(changes.items()):
        db.smf_newforms.update({'label': label}, new, restat=False)
        if i % 500 == 0:
            print("updated", i)
    db.smf_newforms.stats.refresh_stats()


if __name__ == '__main__':
    main(apply='--apply' in sys.argv)
