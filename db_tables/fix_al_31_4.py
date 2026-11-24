"""
Correct the Atkin-Lehner sign of 2.K.31.4.0.a.a, the type G newform of weight (4,0) and level 31,
from -1 to +1, and the Atkin-Lehner dimensions of its space 2.K.31.4.0.a.

fix_al_signs.py set this sign to -1, on the grounds that a spinor-character conductor d = 1 at a
prime level gives -1 (DPRT Theorem 10.1). But this form has d = 31, and its sign is +1. Its
L-function has root number (-1)^k eps_31 with k = 4, and PARI's lfuncheckfeq on 199 coefficients
gives root number +1 (error 2^-57); root number -1 fails.

fix_al_signs.py is not rerun for this: it recomputes ALdims_G of every space with its own divisor
indexing, which no longer matches the stored one (its dry run proposes changes to 674 spaces).
This script changes only the rows below, and checks their current values first.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_al_31_4            # dry run
    sage -python -m smf_lmfdb.db_tables.fix_al_31_4 --apply
"""
import argparse
from lmfdb import db

NEWFORM = '2.K.31.4.0.a.a'
SPACE = '2.K.31.4.0.a'
NF_OLD = {'atkin_lehner_string': '-', 'atkin_lehner_eigenvals': [[31, -1]], 'fricke_eigenval': -1}
NF_NEW = {'atkin_lehner_string': '+', 'atkin_lehner_eigenvals': [[31, 1]], 'fricke_eigenval': 1}
# ALdims are indexed by the squarefree divisors of N: index 0 is sign +, index 1 is sign - at 31.
# The 5 Saito-Kurokawa lifts (ALdims_P = [5, 0]) are unchanged.
SP_OLD = {'ALdims_G': [0, 1], 'ALdims_P': [5, 0], 'ALdims': [5, 1]}
SP_NEW = {'ALdims_G': [1, 0], 'ALdims': [6, 0]}


def check(table, label, expected):
    rec = table.lucky({'label': label}, list(expected))
    if rec is None:
        raise ValueError("%s: no row" % label)
    for k, v in expected.items():
        if rec.get(k) != v:
            raise ValueError("%s: %s is %r, expected %r" % (label, k, rec.get(k), v))


def main(apply=False):
    check(db.smf_newforms, NEWFORM, NF_OLD)
    check(db.smf_newspaces, SPACE, SP_OLD)
    print("%s: %s -> %s" % (NEWFORM, NF_OLD, NF_NEW))
    print("%s: %s -> %s" % (SPACE, {k: SP_OLD[k] for k in SP_NEW}, SP_NEW))
    if not apply:
        print("dry run; pass --apply to write")
        return
    db.smf_newforms.update({'label': NEWFORM}, NF_NEW)
    db.smf_newspaces.update({'label': SPACE}, SP_NEW)
    check(db.smf_newforms, NEWFORM, NF_NEW)
    check(db.smf_newspaces, SPACE, dict(SP_NEW, ALdims_P=SP_OLD['ALdims_P']))
    print("written and read back")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    main(apply=parser.parse_args().apply)
