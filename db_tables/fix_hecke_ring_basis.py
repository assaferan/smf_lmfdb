"""
Rewrite the Hecke ring bases of smf_hecke_nf that do not start with 1.

The newform pages show eigenvalues and Dirichlet coefficients in a basis
1, beta_1, ..., beta_{d-1} of the Hecke ring: the coordinate on the first basis element
is printed as the constant term.  Some stored bases do not start with 1:
  - the paramodular orbits from omf5_data (Sage's basis of the order, as parse_omf5,
    recompute_hecke_rings.py and add_coordinates_14_20.py chose it);
  - the level-2 Yoshida lifts (the integral basis of the compositum, F.integral_basis()).
For these the pages are wrong: for example, a_1 of 2.K.181.3.0.a.a shows as a ten-term
expression.

For each such row this replaces the basis by the triangular basis of the same lattice
(the Hermite form with respect to 1, nu, ..., nu^(d-1); its first element is 1, and the
i-th element is (nu^i + lower terms) / den_i), and rewrites every coordinate column
(lambda_p, lambda_p_square[_0,_1,_2], an) in it.  Before writing anything it checks that:
  - the stored basis and its stored inverse are inverse to each other;
  - the new basis spans the same lattice (the change of basis is in GL_d(Z));
  - every eigenvalue is the same element of the field before and after;
  - a_1 has coordinates (1, 0, ..., 0) in the new basis.
The ring, its index and the field are unchanged, so smf_newforms is not touched.

The update is written to PREFIX_hecke_nf.txt (an update_from_file file keyed by id).

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_hecke_ring_basis [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_hecke_ring_basis [--datafile PREFIX] --apply
"""
import argparse
from collections import Counter

from sage.all import Matrix, QQ, ZZ, lcm, vector
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file

COORD_COLS = ['lambda_p', 'lambda_p_square', 'lambda_p_square_0', 'lambda_p_square_1', 'lambda_p_square_2', 'an']
BASIS_COLS = ['hecke_ring_numerators', 'hecke_ring_denominators',
              'hecke_ring_inverse_numerators', 'hecke_ring_inverse_denominators']


def rows_of(nums, dens):
    return Matrix(QQ, [[QQ(c) / QQ(den) for c in num] for num, den in zip(nums, dens)])


def triangular(B):
    """The basis, starting with 1, of the lattice spanned by the rows of B (coordinates in the
    power basis): the echelon form with respect to nu^(d-1), ..., nu, 1, reversed."""
    d = B.nrows()
    D = lcm([c.denominator() for c in B.list()])
    M = Matrix(ZZ, [[ZZ(c * D) for c in reversed(list(r))] for r in B.rows()]).echelon_form()
    T = Matrix(QQ, [list(reversed(list(r))) for r in reversed(M.rows())]) / D
    if T[0, 0] < 0:
        T.rescale_row(0, -1)
    assert T.row(0) == vector(QQ, [1] + [0] * (d - 1)), "1 is not a primitive element of the lattice"
    return T


def as_cols(M):
    """numerators and denominators of the rows of M, as stored."""
    dens = [int(r.denominator()) for r in M.rows()]
    nums = [[int(c * den) for c in r] for r, den in zip(M.rows(), dens)]
    return nums, dens


def is_one_first(h):
    d = len(h['field_poly']) - 1
    return h['hecke_ring_numerators'][0] == [1] + [0] * (d - 1) and h['hecke_ring_denominators'][0] == 1


def rebasis(h):
    """The new basis and coordinate columns of the smf_hecke_nf row h."""
    d = len(h['field_poly']) - 1
    B = rows_of(h['hecke_ring_numerators'], h['hecke_ring_denominators'])
    Binv = rows_of(h['hecke_ring_inverse_numerators'], h['hecke_ring_inverse_denominators'])
    assert B * Binv == 1, "%s: the stored inverse basis is not the inverse" % h['label']
    T = triangular(B)
    Tinv = T.inverse()
    # x = c B = c' T, so c' = c (B T^-1); C must be in GL_d(Z)
    C = B * Tinv
    assert all(c in ZZ for c in C.list()) and abs(C.det()) == 1, "%s: the new basis spans another lattice" % h['label']
    C = C.change_ring(ZZ)
    new = {}
    for col in COORD_COLS:
        vals = h.get(col)
        if vals is None:
            new[col] = None
            continue
        out = []
        for v in vals:
            if v is None or v == 'NULL':
                out.append(v)
                continue
            w = vector(ZZ, [ZZ(c) for c in v]) * C
            assert w * T == vector(QQ, [QQ(c) for c in v]) * B, h['label']
            out.append([int(c) for c in w])
        new[col] = out
    if new.get('an'):
        assert new['an'][0] == [1] + [0] * (d - 1), "%s: a_1 is not 1 in the new basis" % h['label']
    new['hecke_ring_numerators'], new['hecke_ring_denominators'] = as_cols(T)
    new['hecke_ring_inverse_numerators'], new['hecke_ring_inverse_denominators'] = as_cols(Tinv)
    return new


def main(apply=False, datafile=None):
    base = datafile or "smf_hecke_ring_basis"
    cols = ['id', 'label', 'family', 'level', 'field_poly', 'hecke_ring_power_basis'] + BASIS_COLS + COORD_COLS
    out, stats = [], Counter()
    for h in db.smf_hecke_nf.search({'hecke_ring_power_basis': False}, cols):
        if not h.get('hecke_ring_numerators') or len(h['field_poly']) <= 2 or is_one_first(h):
            continue
        new = rebasis(h)
        new['id'] = h['id']
        out.append(new)
        stats[(h['family'], h['level'] if h['level'] <= 2 else '>2')] += 1
        if len(out) % 250 == 0:
            print("%d rows" % len(out), flush=True)
    print("%d rows to rewrite: %s" % (len(out), dict(stats)))
    upd_cols = BASIS_COLS + COORD_COLS
    # update_from_file needs a value in every column: keep the NULL columns NULL
    write_update_file(base + "_hecke_nf.txt", db.smf_hecke_nf, 'id', out, upd_cols)
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table (undo with reload_revert())
    db.smf_hecke_nf.update_from_file(base + "_hecke_nf.txt", label_col='id')
    print("updated")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--datafile", help="prefix for the output file (default: ./smf_hecke_ring_basis)")
    parser.add_argument("--apply", action="store_true", help="write to the database")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
