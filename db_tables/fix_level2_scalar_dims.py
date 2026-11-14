"""
Correct the dimensions of the scalar-valued level-2 spaces (family P, weight (k, 0)).

smf_dims_degree_2_level_2 sent j = 0 through the vector-valued formulas, whose
generating series start at j = 2: Pick_j(0) reads L[j/2 - 1] = L[-1], the series
of j = 50.  So the 18 spaces 2.P.2.k.0.a (3 <= k <= 20) store the dimensions of
weight (k, 50) -- 2.P.2.8.0.a has cusp_dim 43974 -- except cusp_Y_dim (the
Yoshida part does use j = 0, and is 0).  The scalar-valued formulas
(All_List_Mult_Irrep_SV) were there all along; the generator now uses them.

Before writing anything it checks that the stored values are those of the j = 50
series, and that the new ones agree with what is known independently:
- total_dim = dim M_k(Gamma[2]) from Igusa's ring of theta constants:
  C(n+4, 4) - C(n, 4) for k = 2n, and dim M_{k-5} for k odd (chi_5 M_{k-5});
- for k even, total_dim - cusp_dim = 15 + 15 dim S_k(Gamma_0(4)) (Siegel and
  Klingen Eisenstein series at the 15 + 15 boundary components), and for k odd
  every form is a cusp form;
- the old_* parts (the trivial representation s[6] of S_6, i.e. level 1) are the
  dimensions stored for 2.S.1.k.0.a;
- the S_6 decomposition of the Saito-Kurokawa part is Conjecture 7.6 of
  Bergstrom-Faber-van der Geer (arXiv:0803.0917), from the classical newforms of
  weight 2k - 2 and levels 1, 2 (by Atkin-Lehner sign) and 4 in mf_newspaces.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.fix_level2_scalar_dims [--datafile PREFIX]          # dry run
    sage -python -m smf_lmfdb.db_tables.fix_level2_scalar_dims [--datafile PREFIX] --apply
"""
import argparse

from sage.all import binomial
from lmfdb import db
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file
import smf_lmfdb.db_tables.sage_functions as sf

ISOTYPES = ['total_dim', 'cusp_dim', 'eis_dim', 'eis_F_dim', 'eis_Q_dim', 'cusp_P_dim', 'cusp_Y_dim', 'cusp_G_dim']
COLS = [p + c for c in ISOTYPES for p in ('', 'old_', 'new_')]
# the order of the S_6 irreps in the multiplicity lists
IRREPS = ['6', '51', '42', '411', '33', '321', '3111', '222', '2211', '21111', '111111']


def entry(all_dims):
    """The smf_newspaces columns of smf_dims_degree_2_level_2 from a list of All_List_Mult_Irrep_*."""
    out = {}
    for c, (mults, dim) in zip(ISOTYPES, all_dims):
        out[c], out['old_' + c], out['new_' + c] = int(dim), int(mults[0]), int(dim - mults[0])
    return out


def old_all_dims(k):
    """What the generator computed for weight (k, 0): the j = 50 series, with the Yoshida part of j = 0."""
    L = sf.All_List_Mult_Irrep_VV(k, 50)
    Y = sf.List_Mult_Irrep_Yoshida_Lift_VV(k, 0)
    L[6] = [Y, sf.Dim_List(Y)]
    G = [L[1][0][x] - Y[x] for x in range(11)]
    L[7] = [G, sf.Dim_List(G)]
    return L


def igusa(k):
    if k < 0:
        return 0
    if k % 2:
        return igusa(k - 5)
    n = k // 2
    return binomial(n + 4, 4) - binomial(n, 4)


def classical(w):
    """{1: dim S_w(1), 4: dim S_w(4)^new, '+': , '-': dim S_w^{+/-}(2)^new}, from mf_newspaces."""
    out = {}
    for N in (1, 2, 4):
        r = db.mf_newspaces.lucky({'level': N, 'weight': w, 'char_orbit_index': 1}, ['dim', 'ALdims'])
        assert r is not None or w % 2, "no mf_newspaces row for level %d, weight %d" % (N, w)
        out[N] = r
    plus, minus = out[2]['ALdims'] if out[2]['dim'] else (0, 0)
    return {1: out[1]['dim'], 4: out[4]['dim'], '+': plus, '-': minus}


def sk_conjecture(k):
    """The multiplicities of Conjecture 7.6 (l = m, k = m + 3, classical weight 2k - 2)."""
    t = classical(2*k - 2)
    m = dict.fromkeys(IRREPS, 0)
    if k % 2 == 0:      # m odd
        m['42'] = t['+'] + t[1]
        m['222'] = t['-'] + t[1]
        m['6'] = t[1]
    else:               # m even
        m['33'] = t[4]
        m['51'] = t['-']
        m['111111'] = t['+']
    return [m[i] for i in IRREPS]


def main(apply=False, datafile=None):
    spaces = {s['label']: s for s in db.smf_newspaces.search({'level': 2, 'family': 'P'}, ['label', 'weight'] + COLS)
              if s['weight'][1] == 0}
    level1 = {s['weight'][0]: s for s in db.smf_newspaces.search({'level': 1, 'family': 'S'}, ['label', 'weight'] + ISOTYPES)
              if s['weight'][1] == 0}
    print("%d scalar-valued level-2 spaces" % len(spaces))
    rows = []
    for label, s in sorted(spaces.items(), key=lambda x: x[1]['weight']):
        k = s['weight'][0]
        old, L = entry(old_all_dims(k)), sf.All_List_Mult_Irrep_SV(k)
        new = entry(L)
        assert new == sf.smf_dims_degree_2_level_2(k, 0), label
        bad = [c for c in COLS if s[c] != old[c]]
        assert not bad, "%s: stored %s are not the j = 50 values %s" % (label, {c: s[c] for c in bad}, {c: old[c] for c in bad})
        # independent checks of the new values
        assert new['total_dim'] == igusa(k), "%s: total_dim %d, Igusa %d" % (label, new['total_dim'], igusa(k))
        if k % 2:
            assert new['total_dim'] == new['cusp_dim'], label
        else:
            dimS = classical_S_Gamma0_4(k)
            assert new['total_dim'] - new['cusp_dim'] == 15 + 15*dimS, \
                "%s: %d Eisenstein dimensions, 15 + 15 dim S_k(Gamma_0(4)) = %d" % (label, new['total_dim'] - new['cusp_dim'], 15 + 15*dimS)
        l1 = level1.get(k)
        for c in ['total_dim', 'cusp_dim', 'eis_F_dim', 'eis_Q_dim', 'cusp_P_dim', 'cusp_G_dim']:
            want = (l1[c] or 0) if l1 else 0
            assert new['old_' + c] == want, "%s: old_%s %d, level 1 stores %s" % (label, c, new['old_' + c], want)
        assert list(L[5][0]) == sk_conjecture(k), "%s: SK multiplicities %s, Conjecture 7.6 %s" % (label, L[5][0], sk_conjecture(k))
        rows.append(dict(label=label, **{c: new[c] for c in COLS}))
    for x in rows[:3]:
        print("   ", x)
    print("%d spaces to correct" % len(rows))
    fname = (datafile or "smf_level2_scalar_dims_fix") + "_newspaces.txt"
    write_update_file(fname, db.smf_newspaces, 'label', rows, COLS)
    if not apply:
        print("dry run; pass --apply to write")
        return
    # update_from_file builds a replacement table; undo with db.smf_newspaces.reload_revert().
    db.smf_newspaces.update_from_file(fname, label_col='label')


def classical_S_Gamma0_4(k):
    """dim S_k(Gamma_0(4)) = new(4) + 2 new(2) + 3 dim S_k(1)."""
    t = classical(k)
    return t[4] + 2*(t['+'] + t['-']) + 3*t[1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the corrections to the database")
    parser.add_argument("--datafile", help="prefix for the update file (default: ./smf_level2_scalar_dims_fix)")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
