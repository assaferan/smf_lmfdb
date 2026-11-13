"""
The eigenvalues of T(p^2) of the paramodular newforms of omf5_data.

omf5_data stores, for p <= 7 (weight (3,0); p <= 29 in the other weights), the eigenvalues of T(p,2) in 'lambda_p_square' and
their traces in 'trace_lambda_p_square' (see its README), while smf_newforms and
smf_hecke_nf store eigenvalues of T(p^2):

    T(p^2) = T(p)^2 - p^(mu-1) - (T(p,2) + p^2 + 1) p^(mu-2),    mu = j + 2k - 3.

For an orbit with eigenvalue coordinates this is computed at each embedding.  The
trace of T(p^2) over an orbit also needs the sum of lambda_p^2 over its embeddings,
which the traces of T(p) do not determine; for the orbits stored without
coordinates it is read off the Hecke matrices of omf5_data (hecke_mats, T(p,1) for
p < 10 on each space of spinor conductor d, weight (3,0)): the spaces are split
into blocks by a random combination of the T(p,1), and an orbit is the block with
its dimension and traces of T(p,1).
"""
import os
import random
import re
from functools import lru_cache

from sage.all import Matrix, QQ, ZZ, prime_range
from smf_lmfdb.db_tables.nf_elt import get_nf_basis

OMF5_DIR = "smf_lmfdb/Hecke_Eigenvalues/paramodular/omf5_data"
SQUARE_PRIMES = [2, 3, 5, 7]   # the primes with T(p,2) in omf5_data, and with Hecke matrices
PRIMES = prime_range(100)


def mu(k, j):
    return j + 2*k - 3


def T_p_square(ap, bp2, p, k, j):
    """The eigenvalue of T(p^2) from those of T(p) and T(p,2)."""
    m = mu(k, j)
    return ap**2 - p**(m-1) - (bp2 + p**2 + 1)*p**(m-2)


def trace_T_p_square(sum_sq, trace_bp2, dim, p, k, j):
    """The trace of T(p^2) over an orbit of dimension dim, from the sum of the
    squares of the eigenvalues of T(p) and the trace of T(p,2)."""
    m = mu(k, j)
    return sum_sq - dim*p**(m-1) - (trace_bp2 + dim*(p**2 + 1))*p**(m-2)


def _load(fname):
    with open(fname) as f:
        r = f.read()
    return eval(re.sub(r'(\d)L', r'\1', r))


def orbit_basis(orbit):
    """The basis of the Hecke field in which the eigenvalues are written."""
    return get_nf_basis(orbit)[0]


def orbit_elements(orbit, key, basis):
    return [None if x in ('NULL', None) or len(x) == 0 else sum(QQ(c)*b for c, b in zip(x, basis))
            for x in orbit[key]]


def has_coordinates(orbit):
    return bool(orbit.get('lambda_p')) and bool(orbit.get('lambda_p_square'))


def sums_of_squares_from_coordinates(orbit, N):
    """At every good prime with T(p,2) (p <= 7 in weight (3,0), more in the other weights)."""
    basis = orbit_basis(orbit)
    ap = orbit_elements(orbit, 'lambda_p', basis)
    n = len(orbit['trace_lambda_p_square'])
    return {PRIMES[i]: ZZ((ap[i]**2).trace()) for i in range(n) if N % PRIMES[i] and ap[i] is not None}


@lru_cache(maxsize=4)
def hecke_matrix_blocks(N, folder=OMF5_DIR, seed=1):
    """
    The blocks of the spaces of weight (3,0) and level N under the T(p,1), p < 10
    not dividing N: a list of (d, dim, {p: tr T(p,1)}, {p: tr T(p,1)^2}, irreducible),
    where irreducible says that the block is a single Galois orbit of eigensystems.
    """
    ps = [p for p in SQUARE_PRIMES if N % p]
    if not ps:
        return ()
    mats = {p: _load(os.path.join(folder, "hecke_mats", "hecke_mat_%d_%d.dat" % (N, p))) for p in ps}
    rnd = random.Random(seed)
    out = []
    for d in mats[ps[0]]:
        Ms = {p: Matrix(QQ, mats[p][d]) for p in ps if mats[p][d]}
        if not Ms:
            continue
        for p in Ms:
            for q in Ms:
                if p < q and Ms[p]*Ms[q] != Ms[q]*Ms[p]:
                    raise ValueError("hecke_mat_%d: T(%d,1) and T(%d,1) do not commute on d = %s" % (N, p, q, d))
        T = sum(rnd.randint(1, 50)*Ms[p] for p in Ms)
        for V, irreducible in T.decomposition():
            R = {p: Ms[p].restrict(V) for p in Ms}
            out.append((d, V.dimension(), {p: ZZ(R[p].trace()) for p in Ms}, {p: ZZ((R[p]**2).trace()) for p in Ms},
                        bool(irreducible)))
    return tuple(out)


def sums_of_squares_from_matrices(orbit, N, folder=OMF5_DIR):
    """
    The sum over the orbit of lambda_p^2 for p < 10 not dividing N, from the Hecke
    matrices: the orbit must be the only block with its dimension and traces, and
    that block a single Galois orbit.  Otherwise (few good primes) ValueError.

    It also needs two good primes below 10.  Against the orbits with coordinates
    (5612 that it places) it is right at every one with two or more, but with a
    single good prime the dimension and one trace place an orbit in the wrong
    block 6 times out of 649, so those orbits are refused.
    """
    dim = len(orbit['field_poly']) - 1
    tr = {p: orbit['trace_lambda_p'][PRIMES.index(p)] for p in SQUARE_PRIMES if N % p}
    if not tr:
        return {}
    if len(tr) < 2:
        raise ValueError("level %d: one good prime below 10 cannot place an orbit" % N)
    found = [b for b in hecke_matrix_blocks(N, folder) if b[1] == dim and b[2] == tr]
    if len(found) != 1 or not found[0][4]:
        raise ValueError("level %d: %d blocks of dimension %d with the traces %s of the orbit%s"
                         % (N, len(found), dim, tr, "" if len(found) != 1 else ", not a single Galois orbit"))
    return dict(found[0][3])


def sums_of_squares(orbit, N, k, j, folder=OMF5_DIR):
    """{p: sum of lambda_p^2 over the orbit}, p in SQUARE_PRIMES not dividing N."""
    if has_coordinates(orbit):
        return sums_of_squares_from_coordinates(orbit, N)
    if (k, j) != (3, 0):
        raise ValueError("no eigenvalue coordinates and no Hecke matrices for weight %s" % ((k, j),))
    return sums_of_squares_from_matrices(orbit, N, folder)


def T_p_square_traces(orbit, N, k, j, sum_sq):
    """The traces of T(p^2) at the primes of orbit['trace_lambda_p_square'] (T(p,2) in omf5_data)."""
    dim = len(orbit['field_poly']) - 1
    out = []
    for i, t in enumerate(orbit['trace_lambda_p_square']):
        p = PRIMES[i]
        out.append('NULL' if t in ('NULL', None) or p not in sum_sq else trace_T_p_square(sum_sq[p], t, dim, p, k, j))
    return out


def T_p_square_elements(orbit, k, j):
    """The eigenvalues of T(p^2), as elements of the Hecke field, where omf5_data has T(p,2)."""
    basis = orbit_basis(orbit)
    ap = orbit_elements(orbit, 'lambda_p', basis)
    bp2 = orbit_elements(orbit, 'lambda_p_square', basis)
    return [None if b is None else T_p_square(ap[i], b, PRIMES[i], k, j) for i, b in enumerate(bp2)]


def orbit_T_p_square_traces(orbit, N, k, j, folder=OMF5_DIR):
    """
    (traces of T(p^2), sums of squares) for an omf5_data orbit, in place of its
    traces of T(p,2).  Where the sum of squares cannot be found (an orbit without
    coordinates at a level with a single good prime below 10), the trace is 'NULL'.
    """
    try:
        sq = sums_of_squares(orbit, N, k, j, folder)
    except ValueError:
        sq = {}
    return T_p_square_traces(orbit, N, k, j, sq), sq
