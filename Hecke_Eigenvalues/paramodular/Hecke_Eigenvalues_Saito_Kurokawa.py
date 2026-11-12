"""
Saito-Kurokawa lifts F of classical newforms f of weight w = 2k-2 and trivial
character.

At each embedding the Hecke eigenvalues of F are the SK_* polynomials (of
degree <= 2) in a_p(f), and L(F, s) = zeta(s-k+1) zeta(s-k+2) L(f, s).  The
traces over a Galois orbit are therefore sums over the embeddings, which need
Tr(a_p) and Tr(a_p^2), not the SK_* formulas applied to Tr(a_p).

A lift needs f to have root number (-1)^(w/2) * (Fricke sign) = -1, so the
Fricke sign of f must be (-1)^k (at level 1 the Fricke sign is +1, so k is even).

Conventions (as for the G rows): Dirichlet coefficients are 'NULL' at every n
with a prime factor dividing the level; a trace is 'NULL' when the classical
data does not determine Tr(a_p^2).  At p | N the eigenvalues are the SK_*
formulas, as at the good primes.
"""
from sage.all import Matrix, NumberField, PolynomialRing, QQ, ZZ, divisors, prime_range
from sage.arith.misc import previous_prime

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import MAX_P


def ms_label(k,N):
    return '.'.join([str(N), str(k), "a"])

def SK_lambda_p(k,p,a_p):
    return p**(k-2)+ p**(k-1)+a_p

def SK_lambda_p_square(k,p,a_p):
    return a_p**2+(p+1)*p**(k-2)*a_p+p**(2*k-2)

def SK_lambda_p_square_0(k,p,a_p):
    return a_p**2+p**(k-3)*(p**2-1)*(a_p+p**(k-1))

def SK_lambda_p_square_1(k,p,a_p):
    return p**(k-3)*(p+1)*a_p+p**(2*k-6)*(p**2-1)

def SK_lambda_p_square_2(k,p,a_p):
    return p**(2*k-6)

SK_FUNCS = {'lambda_p' : SK_lambda_p,
            'lambda_p_square' : SK_lambda_p_square,
            'lambda_p_square_0' : SK_lambda_p_square_0,
            'lambda_p_square_1' : SK_lambda_p_square_1,
            'lambda_p_square_2' : SK_lambda_p_square_2}
HECKE_TYPES = list(SK_FUNCS)

CLASSICAL_FORM_COLS = ['label', 'level', 'weight', 'dim', 'nf_label', 'hecke_ring_index', 'field_poly_is_cyclotomic',
                       'field_poly_root_of_unity', 'field_poly_is_real_cyclotomic', 'hecke_ring_index_proved',
                       'hecke_ring_generator_nbound', 'field_disc', 'field_disc_factorization', 'field_poly',
                       'hecke_ring_index_factorization', 'relative_dim', 'traces', 'is_polredabs', 'hecke_orbit_code']
HECKE_RING_COLS = ['field_poly', 'hecke_ring_power_basis', 'hecke_ring_numerators', 'hecke_ring_denominators']

def SK_query(k, N):
    """The classical newforms whose lifts are the SK forms of weight (k,0) and paramodular level N."""
    return {'level' : N, 'weight' : 2*k-2, 'char_order' : 1, 'fricke_eigenval' : (-1)**k}

def SK_coefficients(F, k, p):
    """c0, c1, c2 with F(k, p, a) = c0 + c1 a + c2 a^2."""
    c0, s, m = F(k, p, 0), F(k, p, 1), F(k, p, -1)
    return c0, (s - m) // 2, (s + m - 2 * c0) // 2

def SK_dirichlet(k, N, a, length):
    """
    Coefficients 1..length of zeta(s-k+1) zeta(s-k+2) L(f, s), 'NULL' at n not prime
    to N; a[n-1] is a_n(f), or its trace over the orbit (the map is linear).
    """
    return ['NULL' if ZZ(n).gcd(N) != 1 else
            sum(d1**(k-1) * sum(d2**(k-2) * a[n//d1//d2 - 1] for d2 in divisors(n//d1)) for d1 in divisors(n))
            for n in range(1, length+1)]

def hecke_ring_basis(h):
    """The Hecke ring basis of an mf_hecke_nf row, as elements of its field."""
    K = NumberField(PolynomialRing(QQ, 'x')([int(c) for c in h['field_poly']]), 'nu')
    nu = K.gen()
    if h['hecke_ring_power_basis']:
        return [nu**i for i in range(K.degree())]
    return [sum(int(c) * nu**i for i, c in enumerate(num)) / int(den)
            for num, den in zip(h['hecke_ring_numerators'], h['hecke_ring_denominators'])]

def trace_of_squares(f, primes):
    """{p : Tr(a_p(f)^2)} for the p in primes for which the classical data determines it."""
    N, w, d, tr = f['level'], int(f['weight']), f['dim'], f['traces']
    t2 = {}
    for p in primes:
        # a_{p^2} = a_p^2 - p^(w-1) for p prime to N (trivial character), a_{p^2} = a_p^2 otherwise
        if p*p <= len(tr):
            t2[p] = tr[p*p-1] + (d * p**(w-1) if N % p else 0)
    if all(p in t2 for p in primes):
        return t2
    h = db.mf_hecke_nf.lucky({'hecke_orbit_code' : f['hecke_orbit_code']}, ['ap', 'maxp'] + HECKE_RING_COLS)
    if h is not None:
        B = hecke_ring_basis(h)
        for p, c in zip(prime_range(h['maxp']+1), h['ap']):
            if p in primes and p not in t2:
                a = sum(QQ(x) * b for x, b in zip(c, B))
                t2[p] = ZZ((a*a).trace())
    for x in db.mf_hecke_charpolys.search({'hecke_orbit_code' : f['hecke_orbit_code']}, ['p', 'charpoly_factorization']):
        if x['p'] in primes and x['p'] not in t2:
            s = 0
            for coeffs, mult in x['charpoly_factorization']:  # monic, coefficients from low to high degree
                e = len(coeffs) - 1
                s += mult * (coeffs[e-1]**2 - 2 * (coeffs[e-2] if e >= 2 else 0))
            t2[x['p']] = ZZ(s)
    return t2

def SK_orbit_traces(k, f, prime_bound=MAX_P+1, num_coeffs=MAX_P+1):
    """
    The smf_newforms trace columns of the lift of the classical orbit f: the
    traces of the SK eigenvalues for p < prime_bound, and `traces`, the first
    num_coeffs Dirichlet coefficients.
    """
    primes = prime_range(prime_bound)
    d, tr = f['dim'], f['traces']
    t2 = trace_of_squares(f, primes)
    out = {}
    for ht, F in SK_FUNCS.items():
        vals = []
        for p in primes:
            c0, c1, c2 = SK_coefficients(F, k, p)
            if c2 and p not in t2:
                vals.append('NULL')
            else:
                vals.append(d*c0 + c1*tr[p-1] + (c2*t2[p] if c2 else 0))
        out['trace_' + ht] = vals
    out['traces'] = SK_dirichlet(k, f['level'], tr, num_coeffs)
    return out

def SK_all_forms(k, query, is_cuspidal, prime_bound=MAX_P+1):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = []
    for orbit in db.mf_newforms.search(query, CLASSICAL_FORM_COLS):
        orbit.update(SK_orbit_traces(k, orbit, prime_bound))
        orbit['related_objects'] = [ 'ModularForm/GL2/Q/holomorphic/' + '/'.join(orbit['label'].split('.'))]
        orbit['is_cuspidal'] = is_cuspidal
        orbit['aut_rep_type'] = 'P'
        # We don't want to accidentally use the label or the data of the classical form
        # for the lift ('traces' is now that of the lift)
        for field_name in ['label', 'level', 'weight', 'hecke_orbit_code']:
            dummy = orbit.pop(field_name)
        forms.append(orbit)
    return forms

def SK_space_traces(k, query, prime_bound=200):
    """
    Traces of the SK eigenvalues on the space of lifts of the forms in query,
    { ht : { p : trace } }: the sums over the orbits ('NULL' if one of them is).
    """
    primes = prime_range(prime_bound)
    L = { ht : {p : 0 for p in primes} for ht in HECKE_TYPES }
    for f in db.mf_newforms.search(query, ['level', 'weight', 'dim', 'traces', 'hecke_orbit_code']):
        tr = SK_orbit_traces(k, f, prime_bound, num_coeffs=0)
        for ht in HECKE_TYPES:
            for p, v in zip(primes, tr['trace_' + ht]):
                L[ht][p] = 'NULL' if 'NULL' in (v, L[ht][p]) else L[ht][p] + v
    return L

def SK_all_evs(k, query, prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf: the SK eigenvalues and Dirichlet
    coefficients at each embedding, in the Hecke ring basis of the classical form.
    '''
    # 'an' of mf_hecke_nf only goes up to 100
    evs = []
    maxp = previous_prime(prime_bound)
    for orbit in db.mf_newforms.search(query, ['hecke_orbit_code', 'level']):
        ev = db.mf_hecke_nf.lucky({'hecke_orbit_code' : orbit['hecke_orbit_code']})
        if not ev:
            continue
        B = hecke_ring_basis(ev)
        K = B[0].parent()
        Minv = Matrix(QQ, [list(K(b)) for b in B]).inverse()
        def elt(c):
            return sum(QQ(x) * b for x, b in zip(c, B))
        def coords(x):
            c = list((Matrix(QQ, [list(K(x))]) * Minv)[0])
            assert all(t.is_integer() for t in c), (ev['label'], x)
            return [int(t) for t in c]
        ap = dict(zip(prime_range(ev['maxp']+1), [elt(c) for c in ev['ap']]))
        an = [elt(c) for c in ev['an']]
        ev['maxp'] = ev['maxp_square'] = maxp
        for ht, F in SK_FUNCS.items():
            ev[ht] = [coords(F(k, p, ap[p])) for p in prime_range(maxp+1)]
        ev['an'] = [x if x == 'NULL' else coords(x) for x in SK_dirichlet(k, orbit['level'], an, maxp+1)]
        # We don't want to accidentally use the same fields for the Saito-Kurokawa lift
        for field_name in ['label', 'hecke_orbit_code', 'ap', 'weight']:
            dummy = ev.pop(field_name)
        evs.append(ev)
    return evs

# Paramodular level N

def Hecke_Traces_Eigenvalues_Saito_Kurokawa(k,j,N,prime_bound=200):
    """
    Traces of the Hecke eigenvalues on the space of SK lifts of paramodular level N,
    { ht : { p : trace } }.
    """
    if j > 0:
        return { ht : {p : 0 for p in prime_range(prime_bound)} for ht in HECKE_TYPES }
    return SK_space_traces(k, SK_query(k, N), prime_bound)

def Hecke_Eigenvalues_PM_SK_new_num_forms(k,j,N):
    if (j > 0):
        return 0
    return db.mf_newforms.count(SK_query(k, N))

# We probably don't need this one, but just in case
def Hecke_Eigenvalues_PM_SK_num_forms(k,j,N):
    return sum([Hecke_Eigenvalues_PM_SK_new_num_forms(k,j,M) for M in divisors(N)])

def Hecke_Eigenvalues_PM_SK_all_forms(k,j,N,prime_bound=MAX_P+1):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    if (j > 0):
        return []
    return SK_all_forms(k, SK_query(k, N), True, prime_bound)

def Hecke_Eigenvalues_PM_SK_all_evs(k,j,N,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    if (j > 0):
        return []
    return SK_all_evs(k, SK_query(k, N), prime_bound)
