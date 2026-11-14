"""
Klingen-Eisenstein series (aut_rep_type 'Q', level 1) and Yoshida lifts
(aut_rep_type 'Y', level 2, family P) of classical newforms.

At each embedding:

* the Klingen-Eisenstein series of a newform f of level 1 and weight j+k has
  the eigenvalues KE(k, j, p, a_p(f), a_{p^2}(f)) and spinor L-function
  L(f, s) L(f, s-k+2);
* the Yoshida lift of a pair (f, g) of newforms of level 2, of weights j+2k-2
  and j+2 and opposite Atkin-Lehner signs, has the eigenvalues
  YO(k, j, p, a_p(f), a_p(g)) and L-function L(f, s) L(g, s-k+2).

Both agree with these L-functions: the Euler factor at a good p is
1 - l t + p (l_1 + (p^2+1) l_2) t^2 - ..., with l_i = lambda_p_square_i.

A trace over a Galois orbit is the sum over its embeddings.  The KE_* formulas
are affine in a_{p^2}, so their constant terms count dim times; the YO_* ones
are not linear at all.  Everything is therefore computed at one embedding, in
the coefficient field, from the classical a_n (built from mf_hecke_nf.ap by the
Hecke recursion and checked against mf_hecke_nf.an), and traced.

Conventions, as for the other types: Dirichlet coefficients are 'NULL' at every
n not prime to the level.  For k = 2 the Yoshida square columns hold p times
lambda_p_square_0 and _1, and p^2 times _2 when j = 0 (the values have
denominators otherwise); this has always been so, and is kept.

A Yoshida pair (f, g) with both dimensions > 1 can split into several Galois
orbits, according to the factorization of the field of g over the field of f;
this is not implemented and raises NotImplementedError (there is no such pair
for k, j <= 20).
"""
from sage.all import Matrix, NumberField, PolynomialRing, QQ, ZZ, divisors, factor, floor, prime_range, vector
from sage.arith.misc import previous_prime

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import MAX_P

HECKE_TYPES = ['lambda_p', 'lambda_p_square', 'lambda_p_square_0', 'lambda_p_square_1', 'lambda_p_square_2']
NUM_COEFFS = MAX_P + 1      # length of smf_newforms.traces
RING_COLS = ['field_poly', 'hecke_ring_power_basis', 'hecke_ring_numerators', 'hecke_ring_denominators']
# columns of the classical orbit kept in the smf_newforms row of the lift
FIELD_COLS = ['dim', 'nf_label', 'hecke_ring_index', 'field_poly_is_cyclotomic', 'field_poly_root_of_unity',
              'field_poly_is_real_cyclotomic', 'hecke_ring_index_proved', 'hecke_ring_generator_nbound', 'field_disc',
              'field_disc_factorization', 'field_poly', 'hecke_ring_index_factorization', 'relative_dim', 'is_polredabs']
CLASSICAL_COLS = ['label', 'level', 'weight', 'dim', 'char_order', 'traces', 'hecke_orbit_code', 'atkin_lehner_string']


def KE(k, j, p, a, a2):
    """The Klingen-Eisenstein eigenvalues at one embedding (a = a_p, a2 = a_{p^2}; a2 = None for p^2 > NUM_COEFFS)."""
    p = QQ(p)
    ev = {'lambda_p': (1 + p**(k-2)) * a}
    if a2 is not None:
        ev.update({'lambda_p_square': (1 + p**(k-2) + p**(2*k-4)) * a2 + p**(2*k+j-4) * (p-1),
                   'lambda_p_square_0': (1 + (p-1) * p**(k-3) + p**(2*k-4)) * a2 + p**(k+j-2) * (p**(k-2) * (p-1) - p**(2*k-4) - 1),
                   'lambda_p_square_1': p**(k-3) * a2 + p**(k+j-2) * (1 - p**(k-4) + p**(2*k-4)),
                   'lambda_p_square_2': p**(2*k+j-6) + 0 * a})
    return ev


def YO(k, j, p, a, b):
    """The Yoshida eigenvalues at one embedding (a = a_p(f), b = a_p(g)), with the k = 2 scaling."""
    p = QQ(p)
    ev = {'lambda_p': a + p**(k-2) * b,
          'lambda_p_square': a**2 + p**(k-2) * a * b + p**(2*k-4) * b**2 - p**(2*k+j-4) * (1 + 2*p),
          'lambda_p_square_0': a**2 + p**(2*k-4) * b**2 + a * b * p**(k-3) * (p-1) - 2 * p**(2*k+j-4) * (1+p),
          'lambda_p_square_1': a * b * p**(k-3) + (p**2-1) * p**(2*k+j-6),
          'lambda_p_square_2': p**(2*k+j-6) + 0 * a}
    if k == 2:
        ev['lambda_p_square_0'] *= p
        ev['lambda_p_square_1'] *= p
        if j == 0:
            ev['lambda_p_square_2'] *= p**2
    return ev


def square_bound(prime_bound):
    """The primes p < square_bound(prime_bound) are those with p^2 < prime_bound (as before)."""
    return previous_prime(floor(prime_bound**(QQ(1)/2))) + 1


def ring_basis(h):
    """The field and Hecke ring basis of an mf_hecke_nf row."""
    K = NumberField(PolynomialRing(QQ, 'x')([int(c) for c in h['field_poly']]), 'nu')
    nu = K.gen()
    if h['hecke_ring_power_basis']:
        return K, [nu**i for i in range(K.degree())]
    return K, [sum(int(c) * nu**i for i, c in enumerate(num)) / int(den)
               for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])]


def coordinates(K, B):
    """x -> its (integral) coordinates in the basis B of K."""
    Minv = Matrix(QQ, [K(b).list() for b in B]).inverse()
    def coords(x):
        c = vector(QQ, K(x).list()) * Minv
        assert all(t in ZZ for t in c), "non-integral coordinates %s" % c
        return [int(t) for t in c]
    return coords


def trace(x, L):
    return ZZ(x) if L is QQ else ZZ(L(x).trace())


class Classical:
    """A classical newform orbit of trivial character: its field K and a_1, ..., a_NUM_COEFFS at one embedding."""
    def __init__(self, f):
        self.row = f
        self.label, self.level, self.weight, self.dim = f['label'], f['level'], int(f['weight']), f['dim']
        assert f['char_order'] == 1, "%s: character is not trivial" % self.label
        self.hnf = db.mf_hecke_nf.lucky({'hecke_orbit_code': f['hecke_orbit_code']})
        if self.dim == 1:
            self.K, self.an = QQ, [QQ(x) for x in f['traces'][:NUM_COEFFS]]
            return
        if self.hnf is None:
            raise ValueError("%s has dimension %d and no mf_hecke_nf coordinates" % (self.label, self.dim))
        self.K, B = ring_basis(self.hnf)
        elt = lambda c: sum(QQ(x) * b for x, b in zip(c, B))
        ap = dict(zip(prime_range(self.hnf['maxp'] + 1), [elt(c) for c in self.hnf['ap']]))
        an = [None, self.K(1)] + [None] * (NUM_COEFFS - 1)
        for n in range(2, NUM_COEFFS + 1):
            fa = factor(n)
            if len(fa) > 1:
                an[n] = self.K(1)
                for p, e in fa:
                    an[n] *= an[p**e]
            else:
                p, e = fa[0]
                chi = 0 if self.level % p == 0 else 1
                an[n] = ap[p] if e == 1 else ap[p] * an[p**(e-1)] - chi * p**(self.weight - 1) * an[p**(e-2)]
        self.an = an[1:]
        assert self.an[:len(self.hnf['an'])] == [elt(c) for c in self.hnf['an']], \
            "%s: a_n from the Hecke recursion differ from mf_hecke_nf.an" % self.label
        assert [trace(x, self.K) for x in self.an] == [ZZ(x) for x in f['traces'][:NUM_COEFFS]], self.label


def dirichlet(k, N, A, B, length=NUM_COEFFS):
    """a_n, n <= length, of L(f, s) L(g, s-k+2) at one embedding ('NULL' at n not prime to N)."""
    return ['NULL' if ZZ(n).gcd(N) != 1 else sum(A[d-1] * B[n//d - 1] * QQ(n//d)**(k-2) for d in divisors(n))
            for n in range(1, length + 1)]


class Lift:
    """The lift of one Galois orbit: its field L, the a_n of the two classical forms in L, and its eigenvalues."""
    def __init__(self, typ, k, j, N, L, A, B, fields):
        self.typ, self.k, self.j, self.N, self.L, self.A, self.B = typ, k, j, N, L, A, B
        self.fields = fields    # the mf_newforms row giving the coefficient field
        self.dim = 1 if L is QQ else L.degree()

    def eigenvalues(self, p):
        if self.typ == 'Q':
            return KE(self.k, self.j, p, self.A[p-1], self.A[p*p-1] if p*p <= NUM_COEFFS else None)
        return YO(self.k, self.j, p, self.A[p-1], self.B[p-1])

    def traces(self, prime_bound=MAX_P+1, num_coeffs=NUM_COEFFS):
        """The smf_newforms trace columns and `traces`."""
        sq = prime_bound if self.typ == 'Y' else square_bound(prime_bound)
        out = {}
        for ht in HECKE_TYPES:
            ps = prime_range(prime_bound if ht == 'lambda_p' else sq)
            out['trace_' + ht] = [trace(self.eigenvalues(p)[ht], self.L) for p in ps]
        out['traces'] = [c if c == 'NULL' else trace(c, self.L) for c in dirichlet(self.k, self.N, self.A, self.B, num_coeffs)]
        return out


# Klingen-Eisenstein series

def KE_query(k, j):
    return {'level': 1, 'weight': j + k, 'char_order': 1}

def KE_lifts(k, j, cols=CLASSICAL_COLS):
    for f in db.mf_newforms.search(KE_query(k, j), cols):
        c = Classical(f)
        yield c, Lift('Q', k, j, 1, c.K, c.an, c.an, f)

def KE_space_traces(k, j, e=0, prime_bound=MAX_P+1):
    """{ ht : { p : trace } } on the Klingen-Eisenstein series of weight (k, j) and level 1
    (none with character, e = 1, as in KE_num_forms' callers)."""
    sq = square_bound(prime_bound)
    L = {ht: {p: 0 for p in prime_range(prime_bound if ht == 'lambda_p' else sq)} for ht in HECKE_TYPES}
    if k < 4 or e == 1:
        return L
    for _, lift in KE_lifts(k, j):
        tr = lift.traces(prime_bound, num_coeffs=0)
        for ht in HECKE_TYPES:
            for p, v in zip(L[ht], tr['trace_' + ht]):
                L[ht][p] += v
    return L

def KE_num_forms(k, j):
    return 0 if k < 4 else db.mf_newforms.count(KE_query(k, j))

def KE_all_forms(k, j, prime_bound=MAX_P+1):
    """The smf_newforms rows of the Klingen-Eisenstein series of weight (k, j) and level 1."""
    if k < 4:
        return []
    forms = []
    for c, lift in KE_lifts(k, j, CLASSICAL_COLS + FIELD_COLS):
        orbit = {col: c.row[col] for col in FIELD_COLS if col in c.row}
        orbit.update(lift.traces(prime_bound))
        orbit['related_objects'] = ['ModularForm/GL2/Q/holomorphic/' + '/'.join(c.label.split('.'))]
        orbit['is_cuspidal'] = False
        orbit['aut_rep_type'] = 'Q'
        forms.append(orbit)
    return forms

def KE_all_evs(k, j, prime_bound=100):
    """The smf_hecke_nf rows: the eigenvalues and a_n at one embedding, in the Hecke ring basis of f."""
    if k < 4:
        return []
    evs = []
    maxp, maxp_square = previous_prime(prime_bound), square_bound(prime_bound) - 1
    for c, lift in KE_lifts(k, j):
        if c.hnf is None:
            continue
        ev = dict(c.hnf)
        K, B = ring_basis(ev)
        coords = coordinates(K, B)
        ev['maxp'], ev['maxp_square'] = maxp, maxp_square
        for ht in HECKE_TYPES:
            ps = prime_range((maxp if ht == 'lambda_p' else maxp_square) + 1)
            ev[ht] = [coords(lift.eigenvalues(p)[ht]) for p in ps]
        ev['an'] = [x if x == 'NULL' else coords(x) for x in dirichlet(k, 1, lift.A, lift.B, maxp + 1)]
        # We don't want to accidentally use the same fields for the Klingen-Eisenstein series
        for field_name in ['id', 'label', 'hecke_orbit_code', 'ap', 'weight']:
            ev.pop(field_name, None)
        evs.append(ev)
    return evs


# Yoshida lifts

def Y_pairs(k, j, cols=CLASSICAL_COLS):
    """
    The pairs (f, g) of level-2 newforms of weights j+2k-2 and j+2 and opposite
    Atkin-Lehner signs, in the order the rows have always been numbered.  For
    k = 2 the two weights agree and each pair is listed once.
    """
    ws = [j + 2, j + 2*k - 2]
    signs = ['+', '-']
    search = lambda w, s: list(db.mf_newforms.search({'level': 2, 'weight': w, 'char_order': 1, 'atkin_lehner_string': s}, cols))
    if ws[0] != ws[1]:
        res = [[search(w, s) for w in ws] for s in signs]
        return [(f, g) for i in [0, 1] for f in res[1-i][1] for g in res[i][0]]
    res = [search(ws[0], s) for s in signs]
    return [(f, g) for f in res[1] for g in res[0]]

def Y_lift(k, j, f, g):
    cf, cg = Classical(f), Classical(g)
    if cf.dim > 1 and cg.dim > 1:
        raise NotImplementedError("Yoshida lift of %s and %s: both have dimension > 1" % (cf.label, cg.label))
    if cg.dim == 1:
        L, A, B, fields = cf.K, cf.an, [cf.K(x) for x in cg.an] if cf.K is not QQ else cg.an, cf.row
    else:
        L, A, B, fields = cg.K, [cg.K(x) for x in cf.an], cg.an, cg.row
    return cf, cg, Lift('Y', k, j, 2, L, A, B, fields)

def Y_valid(k, j, e):
    return e > 0 and j % 2 == 0 and k >= 2

def Y_space_traces(k, j, e, prime_bound=MAX_P+1):
    """{ ht : { p : trace } } on the Yoshida lifts of weight (k, j) and level 2."""
    L = {ht: {p: 0 for p in prime_range(prime_bound)} for ht in HECKE_TYPES}
    if not Y_valid(k, j, e):
        return L
    for f, g in Y_pairs(k, j):
        tr = Y_lift(k, j, f, g)[2].traces(prime_bound, num_coeffs=0)
        for ht in HECKE_TYPES:
            for p, v in zip(L[ht], tr['trace_' + ht]):
                L[ht][p] += v
    return L

def Y_num_forms(k, j, e):
    if not Y_valid(k, j, e):
        return 0
    return len(Y_pairs(k, j, ['label']))

def Y_all_forms(k, j, e, prime_bound=MAX_P+1):
    """The smf_newforms rows of the Yoshida lifts of weight (k, j) and level 2."""
    if not Y_valid(k, j, e):
        return []
    forms = []
    for f, g in Y_pairs(k, j, CLASSICAL_COLS + FIELD_COLS):
        cf, cg, lift = Y_lift(k, j, f, g)
        orbit = {col: lift.fields[col] for col in FIELD_COLS if col in lift.fields}
        orbit['dim'] = lift.dim
        orbit.update(lift.traces(prime_bound))
        orbit['related_objects'] = ['ModularForm/GL2/Q/holomorphic/' + '/'.join(c.label.split('.')) for c in [cf, cg]]
        orbit['is_cuspidal'] = True
        orbit['aut_rep_type'] = 'Y'
        forms.append(orbit)
    return forms

def Y_all_evs(k, j, e, prime_bound=100):
    """
    The smf_hecke_nf rows: the eigenvalues and a_n at one embedding, in the
    integral basis of the composite of the fields of f and g (as before).
    """
    if not Y_valid(k, j, e):
        return []
    evs = []
    maxp = previous_prime(prime_bound)
    for f, g in Y_pairs(k, j):
        cf, cg, lift = Y_lift(k, j, f, g)
        if cf.hnf is None or cg.hnf is None:
            continue
        F0, F1 = ring_basis(cf.hnf)[0], ring_basis(cg.hnf)[0]
        F, to_F1, to_F0, _ = F1.composite_fields(F0, both_maps=True)[0]
        A = [to_F0(F0(x)) for x in cf.an]
        Bv = [to_F1(F1(x)) for x in cg.an]
        nu = F.gen()
        basis = F.integral_basis()
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
        for ht in HECKE_TYPES:
            ev[ht] = [coords(YO(k, j, p, A[p-1], Bv[p-1])[ht]) for p in prime_range(maxp + 1)]
        ev['an'] = [x if x == 'NULL' else coords(x) for x in dirichlet(k, 2, A, Bv, maxp + 1)]
        # We don't want to accidentally use the same fields for the Yoshida lift
        for field_name in ['id', 'label', 'hecke_orbit_code', 'ap', 'weight']:
            ev.pop(field_name, None)
        evs.append(ev)
    return evs
