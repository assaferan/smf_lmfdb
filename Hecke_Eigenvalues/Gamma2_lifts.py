"""
The lifts in the spaces of Siegel modular forms of degree 2 and principal level 2
(family P, level 2: Gamma[2] = Gamma_2[2]), with their multiplicities.

Sp(4, Z)/Gamma[2] = Sp(4, F_2) = S_6 acts on M_{k,j}(Gamma[2]) and commutes with
the Hecke operators at odd p, so an eigensystem occurs with multiplicity the
dimension of its S_6-isotypic part.  The new part (everything but the trivial
representation s[6], which is level 1) decomposes, after Bergstrom-Faber-van der
Geer (arXiv:0803.0917; Conjectures 7.4 and 7.6, proved by Roesner), as follows.
Here h, f, g are classical newforms of trivial character, and tau^{+/-} are the
level-2 newforms with w_2 = +1 / -1.

* F (Siegel Eisenstein series, j = 0, k >= 4 even): the level-1 eigensystem, at
  the 15 boundary components, 14 of them new.
* Q (Klingen-Eisenstein series, j even, k >= 4): 15 dim S_{j+k}(Gamma_0(4)), in
  which h of level 1, 2, 4 and weight j + k occurs 3, 2, 1 times; multiplicity
  45 - 1 (the level-1 copy is old), 30, 15.
* P (Saito-Kurokawa lifts, j = 0), Conj. 7.6, f of weight 2k - 2:
  k even: tau^+ in s[4,2] (9), tau^- in s[2^3] (5), level 1 in s[4,2] + s[2^3] (14);
  k odd:  tau^+ in s[1^6] (1), tau^- in s[5,1] (5), level 4 in s[3^2] (5).
* Y (Yoshida lifts, j > 0), Conj. 7.4, (f, g) of weights j + 2k - 2 and j + 2:
  level-2 pairs of opposite signs in s[1^6] (1), of equal signs in s[2^3] (5),
  level-4 pairs in s[2,1^4] (5).

For k = 2 (j > 0) the Yoshida part is Clery-van der Geer, Thm 1.1 (arXiv:1709.01748),
with f, g of the same weight j + 2 and the pairs unordered:
  level-2 pairs of opposite signs in s[1^6] (1);
  pairs of distinct eigenforms of level 2 and equal signs in s[2^3] (5), of dimension
  dim Lambda^2 S^+ + dim Lambda^2 S^-;
  pairs of distinct eigenforms of level 4 in s[2,1^4] (5), of dimension dim Lambda^2 S^new.
A pair of distinct eigenforms in one Galois orbit f (f and a conjugate) is a lift
of rational eigenvalues when dim f = 2; larger orbits are not handled.  Conj. 1.2
says that these lifts span S_{j,2}(Gamma[2]).

In each space the sum of dim * multiplicity over each type is the stored new_*_dim.

Conventions: the eigenvalues at p = 2 are the formulas evaluated at the classical
a_2, as in the rows stored before (Gamma[2] has no T(2), and the site does not show
them).  The Dirichlet series is the spin L-function of the automorphic
representation, which is a product of classical L-functions, Euler factor at 2
included:
  F: zeta(s) zeta(s-k+1) zeta(s-k+2) zeta(s-2k+3)  (that of level 1),
  Q: L(h, s) L(h, s-k+2),
  P: zeta(s-k+1) zeta(s-k+2) L(f, s),
  Y: L(f, s) L(g, s-k+2).
"""
import json
import os

from sage.all import Matrix, NumberField, PolynomialRing, QQ, ZZ, cyclotomic_polynomial, euler_phi, pari, prime_range
from sage.arith.misc import previous_prime

from lmfdb import db
from smf_lmfdb.Hecke_Eigenvalues.Klingen_Yoshida import (
    CLASSICAL_COLS, FIELD_COLS, HECKE_TYPES, YO, Classical, Lift, KE_valid, Y_lift,
    coordinates, dirichlet, integral_basis_one_first, ordered_eigenvalues, ring_basis, square_bound)
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import (
    SK_all_forms, SK_all_evs)
from smf_lmfdb.db_tables.common_populate import MAX_P

LEVEL = 2
TYPE_ORDER = 'FQPYG'   # the order of the types in a space, as at level 1
HOLO = 'ModularForm/GL2/Q/holomorphic/'

# polredabs is slow on the larger absolute fields of the split Yoshida lifts with
# older PARI (2.15): its results, which do not depend on the database, can be kept
# in the JSON file named by SMF_POLREDABS_CACHE, made where PARI is faster.
POLREDABS_CACHE = os.environ.get('SMF_POLREDABS_CACHE')
_polredabs = {}
if POLREDABS_CACHE and os.path.exists(POLREDABS_CACHE):
    with open(POLREDABS_CACHE) as fp:
        _polredabs = json.load(fp)


def polredabs(poly):
    """[P, a] for a rational polynomial poly: P is polredabs(poly), and a (coefficients,
    as strings, constant first) gives a root of poly as a polynomial in a root of P."""
    poly = PolynomialRing(QQ, 'x')(poly)
    key = json.dumps([str(c) for c in poly.list()])
    if key not in _polredabs:
        P, a = pari(poly).polredabs(1)
        _polredabs[key] = [[int(c) for c in P.Vec()[::-1]],
                           [str(QQ(c)) for c in a.lift().Vec()[::-1]]]
    return _polredabs[key]


def save_polredabs_cache():
    if POLREDABS_CACHE:
        with open(POLREDABS_CACHE, 'w') as fp:
            json.dump(_polredabs, fp)


def classical_query(N, w, sign=None):
    q = {'level': N, 'weight': w, 'char_order': 1}
    if sign is not None:
        q['atkin_lehner_string'] = '+' if sign == 1 else '-'
    return q


def classical(N, w, sign=None, cols=CLASSICAL_COLS + FIELD_COLS):
    if w < 2 or w % 2:
        return []
    cols = list(dict.fromkeys(cols))     # without repeats ('dim' is in both lists)
    return list(db.mf_newforms.search(classical_query(N, w, sign), cols))


def related(label):
    return HOLO + '/'.join(label.split('.'))


# The eigensystems

def eigensystems(k, j):
    """
    The lift eigensystems of the new part of M_{k,j}(Gamma[2]), as
    (type, multiplicity, data): data is ('F', k), ('Q', h), ('P', query) or ('Y', f, g),
    with h, f, g mf_newforms rows (for P, an mf_newforms query of one orbit).
    """
    E = []
    if j == 0 and k >= 4 and k % 2 == 0:
        E.append(('F', 14, ('F', k)))
    if KE_valid(k, j):
        for N, mult in ((1, 44), (2, 30), (4, 15)):
            for h in classical(N, j + k):
                E.append(('Q', mult, ('Q', h)))
    if j == 0 and k >= 3:
        w = 2*k - 2
        parts = [(2, 1, 9), (2, -1, 5), (1, None, 14)] if k % 2 == 0 else [(2, 1, 1), (2, -1, 5), (4, None, 5)]
        for N, sign, mult in parts:
            for f in classical(N, w, sign, ['label']):
                E.append(('P', mult, ('P', {'label': f['label']})))
    if j > 0 and k >= 3:
        w1, w2 = j + 2*k - 2, j + 2
        for s1 in (1, -1):
            for s2 in (1, -1):
                for f in classical(2, w1, s1):
                    for g in classical(2, w2, s2):
                        E.append(('Y', 5 if s1 == s2 else 1, ('Y', f, g)))
        for f in classical(4, w1):
            for g in classical(4, w2):
                E.append(('Y', 5, ('Y', f, g)))
    if j > 0 and k == 2:
        w = j + 2
        plus, minus = classical(2, w, 1), classical(2, w, -1)
        for f in minus:     # as Klingen_Yoshida.Y_pairs
            for g in plus:
                E.append(('Y', 1, ('Y', f, g)))
        for forms in (plus, minus, classical(4, w)):
            for i, f in enumerate(forms):
                if f['dim'] > 1:
                    E.append(('Y', 5, ('Y', f, f)))     # f with a conjugate
                for g in forms[i+1:]:
                    E.append(('Y', 5, ('Y', f, g)))
    return E


# The smf_newforms fields (as in KE_all_forms, SK_all_forms, Y_all_forms) and the
# smf_hecke_nf fields (as in KE_all_evs, SK_all_evs, Y_all_evs) of one eigensystem

def siegel_eisenstein(k):
    """The level-1 row of E_k and its smf_hecke_nf row; the eigenvalues are those of level 1."""
    f = db.smf_newforms.lucky({'space_label': '2.P.1.%d.0.a' % k, 'aut_rep_type': 'F'}, 2)
    assert f is not None, k
    h = db.smf_hecke_nf.lucky({'hecke_orbit_code': f['hecke_orbit_code']}, 2)
    row = {c: f[c] for c in FIELD_COLS + ['trace_' + ht for ht in HECKE_TYPES] + ['related_objects'] if c in f}
    row['traces'] = f['traces']
    row['is_cuspidal'] = False
    row['aut_rep_type'] = 'F'
    ev = None
    if h is not None:
        ev = {c: v for c, v in h.items() if c not in ('id', 'label', 'hecke_orbit_code', 'level', 'weight', 'family')}
    return row, ev


def klingen_eisenstein(k, j, h):
    c = Classical(h)
    lift = Lift('Q', k, j, 1, c.K, c.an, c.an, h)    # N = 1: no 'NULL' Dirichlet coefficients
    row = {col: h[col] for col in FIELD_COLS if col in h}
    row.update(lift.traces(MAX_P + 1))
    row['related_objects'] = [related(c.label)]
    row['is_cuspidal'] = False
    row['aut_rep_type'] = 'Q'
    ev = None
    if c.hnf is not None:
        prime_bound = 100
        ev = dict(c.hnf)
        K, B = ring_basis(ev)
        coords = coordinates(K, B)
        maxp, maxp_square = previous_prime(prime_bound), square_bound(prime_bound) - 1
        ev['maxp'], ev['maxp_square'] = maxp, maxp_square
        for ht in HECKE_TYPES:
            ps = prime_range((maxp if ht == 'lambda_p' else maxp_square) + 1)
            ev[ht] = [coords(lift.eigenvalues(p)[ht]) for p in ps]
        ev['an'] = [coords(x) for x in dirichlet(k, 1, lift.A, lift.B, maxp + 1)]
        for name in ['id', 'label', 'hecke_orbit_code', 'ap', 'weight', 'level']:
            ev.pop(name, None)
    return row, ev


def saito_kurokawa(k, query):
    rows = SK_all_forms(k, query, True, full=True)
    assert len(rows) == 1, (k, query)
    row = rows[0]
    evs = SK_all_evs(k, query, full=True)
    ev = evs[0] if evs else None
    if ev is not None:
        for name in ['id', 'level']:
            ev.pop(name, None)
    return row, ev


def yoshida(k, j, f, g):
    """[(row, ev)]: one Galois orbit, or several when f and g both have dimension > 1."""
    if f['label'] == g['label']:
        return [yoshida_conjugates(k, j, f)]
    if f['dim'] > 1 and g['dim'] > 1:
        return yoshida_split(k, j, f, g)
    return [yoshida_one(k, j, f, g)]


def yoshida_one(k, j, f, g):
    cf, cg, lift = Y_lift(k, j, f, g)
    lift.N = 1    # no 'NULL' Dirichlet coefficients
    row = {col: lift.fields[col] for col in FIELD_COLS if col in lift.fields}
    row['dim'] = lift.dim
    row.update(lift.traces(MAX_P + 1))
    row['related_objects'] = [related(c.label) for c in [cf, cg]]
    row['is_cuspidal'] = True
    row['aut_rep_type'] = 'Y'
    ev = None
    if cf.hnf is not None and cg.hnf is not None:
        F0, F1 = ring_basis(cf.hnf)[0], ring_basis(cg.hnf)[0]
        F, to_F1, to_F0, _ = F1.composite_fields(F0, both_maps=True)[0]
        ev = yoshida_ev(k, j, cf, F, [to_F0(F0(x)) for x in cf.an], [to_F1(F1(x)) for x in cg.an])
    return row, ev


def yoshida_conjugates(k, j, f):
    """
    The lift of (f, f^sigma) for f of dimension 2 and sigma the conjugation of its
    field K (k = 2): its eigenvalues and a_n are symmetric in the two, so rational.
    """
    cf = Classical(f)
    if cf.dim != 2:
        raise NotImplementedError("the lifts of pairs of conjugates of %s, of dimension %d" % (cf.label, cf.dim))
    K = cf.K
    sigma = [s for s in K.automorphisms() if s(K.gen()) != K.gen()][0]
    A = cf.an
    B = [sigma(x) for x in A]
    lift = Lift('Y', k, j, 1, QQ, A, B, None)    # trace(x, QQ) = ZZ(x): the values are rational elements of K
    rational = db.mf_newforms.lucky({'dim': 1, 'char_order': 1}, FIELD_COLS)
    row = {col: rational[col] for col in FIELD_COLS if col in rational}
    row.update(lift.traces(MAX_P + 1))
    row['related_objects'] = [related(cf.label)] * 2
    row['is_cuspidal'] = True
    row['aut_rep_type'] = 'Y'
    ev = None
    if cf.hnf is not None:
        Q1 = NumberField(PolynomialRing(QQ, 'x')([0, 1]), 'nu')
        ev = yoshida_ev(k, j, cf, Q1, A, B, to_F=lambda x: Q1(QQ(x)))
    return row, ev


def field_columns(L):
    """The smf_newforms field columns of an absolute field L given by a polredabs polynomial.
    The Hecke ring of the lift is not computed: its index columns are NULL."""
    poly = [int(c) for c in L.defining_polynomial().list()]
    d = L.degree()
    disc = ZZ(L.disc())
    def same_field(h):
        return polredabs(h)[0] == poly
    def real_cyclotomic(n):
        z = NumberField(cyclotomic_polynomial(n), 'z').gen()
        return (z + 1/z).minpoly()
    root = next((n for n in range(3, 4*d*d + 3) if euler_phi(n) == d and same_field(cyclotomic_polynomial(n))), 0)
    real = root == 0 and any(same_field(real_cyclotomic(n)) for n in range(5, 8*d*d + 5) if euler_phi(n) == 2*d)
    nf = db.nf_fields.lucky({'coeffs': poly}, 'label')
    return {'dim': d, 'relative_dim': d, 'field_poly': poly, 'is_polredabs': True, 'nf_label': nf,
            'field_disc': disc, 'field_disc_factorization': [[int(p), int(e)] for p, e in disc.factor()],
            'field_poly_is_cyclotomic': root > 0, 'field_poly_root_of_unity': root,
            'field_poly_is_real_cyclotomic': real,
            'hecke_ring_index': None, 'hecke_ring_index_proved': None, 'hecke_ring_index_factorization': None,
            'hecke_ring_generator_nbound': None}


def yoshida_split(k, j, f, g):
    """
    The Galois orbits of the Yoshida lift of (f, g), both of dimension > 1: one for
    each irreducible factor h of the polynomial of the field of g over the field K
    of f, with field K[y]/h (given by its polredabs polynomial).
    """
    cf, cg = Classical(f), Classical(g)
    K, Kg = cf.K, cg.K
    out = []
    for h, e in Kg.defining_polynomial().change_ring(K).factor():
        assert e == 1, (cf.label, cg.label)
        Lrel = K.extension(h, 'b')
        Labs = Lrel.absolute_field('c')
        from_abs, to_abs = Labs.structure()
        red, a = polredabs(Labs.defining_polynomial())
        L = NumberField(PolynomialRing(QQ, 'x')(red), 'nu')
        to_L = Labs.hom([L(PolynomialRing(QQ, 'x')([QQ(c) for c in a]))])
        b = to_L(to_abs(Lrel.gen()))
        A = [to_L(to_abs(Lrel(x))) for x in cf.an]
        B = [Kg(x).polynomial()(b) for x in cg.an]
        lift = Lift('Y', k, j, 1, L, A, B, None)    # N = 1: no 'NULL' Dirichlet coefficients
        row = field_columns(L)
        row.update(lift.traces(MAX_P + 1))
        row['related_objects'] = [related(c.label) for c in [cf, cg]]
        row['is_cuspidal'] = True
        row['aut_rep_type'] = 'Y'
        ev = None
        if cf.hnf is not None and cg.hnf is not None:
            ev = yoshida_ev(k, j, cf, L, A, B)
        out.append((row, ev))
    assert sum(r['dim'] for r, _ in out) == cf.dim * cg.dim, (cf.label, cg.label)
    return out


def yoshida_ev(k, j, cf, F, A, Bv, prime_bound=100, to_F=None):
    """As in Y_all_evs: the eigenvalues in the field F (containing the a_n A of f and Bv of g), in a basis
    of its maximal order starting with 1 (integral_basis_one_first).
    With to_F, A and Bv lie elsewhere and the eigenvalues and a_n are mapped to F by to_F."""
    maxp = previous_prime(prime_bound)
    lams = ordered_eigenvalues(k, j, A, Bv, maxp)
    if to_F is not None:
        lams = [to_F(x) for x in lams]
    basis = integral_basis_one_first(F, lams)
    mat = Matrix([list(b) for b in basis])
    ev = dict(cf.hnf)
    ev['hecke_ring_denominators'] = [row.denominator() for row in mat]
    ev['hecke_ring_numerators'] = [list(row.denominator() * row) for row in mat]
    ev['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
    ev['hecke_ring_inverse_numerators'] = [list(row.denominator() * row) for row in mat**(-1)]
    ev['field_poly'] = list(F.defining_polynomial())
    ev['hecke_ring_power_basis'] = False
    ev['hecke_ring_cyclotomic_generator'] = 0
    ev['hecke_ring_rank'] = F.degree()
    ev['maxp'] = ev['maxp_square'] = maxp
    coords = coordinates(F, basis)
    if to_F is not None:
        coords = (lambda c: lambda x: c(to_F(x)))(coords)
    for ht in HECKE_TYPES:
        ev[ht] = [coords(YO(k, j, p, A[p-1], Bv[p-1])[ht]) for p in prime_range(maxp + 1)]
    ev['an'] = [coords(x) for x in dirichlet(k, 1, A, Bv, maxp + 1)]
    for name in ['id', 'label', 'hecke_orbit_code', 'ap', 'weight', 'level']:
        ev.pop(name, None)
    return ev


def lift_data(k, j, data):
    """[(smf_newforms fields, smf_hecke_nf fields or None)]: the Galois orbits of one eigensystem."""
    typ = data[0]
    if typ == 'F':
        return [siegel_eisenstein(k)]
    if typ == 'Q':
        return [klingen_eisenstein(k, j, data[1])]
    if typ == 'P':
        return [saito_kurokawa(k, data[1])]
    return yoshida(k, j, data[1], data[2])


def sources(data):
    """The classical forms an eigensystem lifts, as stored in related_objects."""
    typ = data[0]
    if typ == 'F':
        return None
    if typ in 'PQ':
        return (related(data[1]['label']),)
    return (related(data[1]['label']), related(data[2]['label']))


def check_sum(E):
    """{type: sum of dim * multiplicity}, from the classical dimensions."""
    out = {}
    for typ, mult, data in E:
        if typ == 'F':
            d = 1
        elif typ == 'P':
            d = db.mf_newforms.lucky(data[1], 'dim')
        elif typ == 'Q':
            d = data[1]['dim']
        elif data[1]['label'] == data[2]['label']:
            d = data[1]['dim'] * (data[1]['dim'] - 1) // 2
        else:
            d = data[1]['dim'] * data[2]['dim']
        out[typ] = out.get(typ, 0) + ZZ(d) * mult
    return out
