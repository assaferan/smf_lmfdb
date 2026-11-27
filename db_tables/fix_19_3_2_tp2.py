"""
Fix T(p^2) and a_n of 2.K.19.3.2.a.a (weight (3,2), mu = 5), as fix_31_4_tp2 did for
2.K.31.4.0.a.a.  The stored T(p^2) at the good primes p was computed as
T^2 - p^4 - p^3 T2 - p^3 (p^2 + 1) instead of T^2 - p^4 - p T2 - p^3 (p^2 + 1), T2 = omf5's
T(p,2).  The two agree only when T2 = 0 (here p = 2).  Invert the stored value to get T2
(checked against omf5 3564131 -square runs for p <= 29), recompute T(p^2), and recompute a_n
at the good primes from the Euler factor
1 - a_p x + b x^2 - p^mu a_p x^3 + p^(2 mu) x^4, b = p T2 + p^3 (p^2 + 1).
The values at p = 19 (lambda_19 = -950 and its T(p^2) entry) and a_{19^k} are kept.

Evidence: with the stored a_n the spin L-function fails glfunc's certification (ERR_RH_ERROR,
ERR_DBL_ZERO); with the corrected factors it certifies (4-19-1.1-c3c5-0-0, root number 1).

    sage -python -m smf_lmfdb.db_tables.fix_19_3_2_tp2 [--apply]

Run fill_trace_tables afterwards: smf_hecke_traces and smf_hecke_newspace_traces hold Tr(a_n).
"""
import sys
from sage.all import prime_range, factor, ComplexField, PolynomialRing, Integer
from lmfdb import db

L, SP, N, MU = '2.K.19.3.2.a.a', '2.K.19.3.2.a', 19, 5
# omf5 3564131, -weight=1,1 -disc=19 -p=p -square (lava ~/omf5b/pilotw/1,1/19/s<p>.txt)
OMF5_T2 = {2: 0, 3: -54, 5: -576, 7: 220, 11: -20436, 13: -31598, 17: 35730, 23: -24324, 29: -422310}
nf = db.smf_newforms.lookup(L)
h = db.smf_hecke_nf.lookup(L)
sp = db.smf_newspaces.lookup(SP, ['id', 'cusp_G_lambda_p_square', 'new_cusp_G_dim'])
assert sp['new_cusp_G_dim'] == 1
ps = prime_range(100)
T = dict(zip(ps, nf['trace_lambda_p']))
sq = [p for p in ps if p <= h['maxp_square']]
old = dict(zip(sq, nf['trace_lambda_p_square']))
assert [[x] for x in nf['trace_lambda_p_square']] == h['lambda_p_square'] == [[x] for x in sp['cusp_G_lambda_p_square']]
assert [[x] for x in nf['traces']] == h['an']
good = [p for p in sq if p != N]
T2 = {}
for p in good:
    num = T[p]**2 - p**4 - p**3 * (p**2 + 1) - old[p]
    assert num % p**3 == 0, p
    T2[p] = num // p**3
    assert T2[p] == OMF5_T2[p], (p, T2[p], OMF5_T2[p])
new = dict(old)
for p in good:
    new[p] = T[p]**2 - p**4 - p * T2[p] - p**3 * (p**2 + 1)
an_old = nf['traces']
M = len(an_old)


def local(p, kmax):
    if p == N:
        return [an_old[p**k - 1] if p**k <= M else None for k in range(kmax + 1)]
    ap = T[p]
    assert an_old[p - 1] == ap, p
    b = None
    if kmax >= 2:
        assert p in new, ("need T(p^2) for", p)
        b = ap**2 - (new[p] + p**(MU - 1))
    a = [1, ap]
    for k in range(2, kmax + 1):
        v = ap * a[k-1] - b * a[k-2]
        if k >= 3:
            v += p**MU * ap * a[k-3]
        if k >= 4:
            v -= p**(2*MU) * a[k-4]
        a.append(v)
    return a


cache = {}
an_new = []
for n in range(1, M + 1):
    v = 1
    for p, e in (factor(n) if n > 1 else []):
        if (p, e) not in cache:
            cache[(p, e)] = local(int(p), int(e))[int(e)]
        v *= cache[(p, e)]
    an_new.append(int(v))
# entries with no good p^2 | n must be unchanged
for n in range(1, M + 1):
    if all(p == N or e == 1 for p, e in (factor(n) if n > 1 else [])):
        assert an_new[n-1] == an_old[n-1], n
# Ramanujan on the corrected local factors (unitary normalization: roots of the reversed
# polynomial have absolute value p^(mu/2), i.e. those in x have absolute value p^(-mu/2))
C = ComplexField(600)
Cx = PolynomialRing(C, 'x')
for p in good:
    b = T[p]**2 - (new[p] + p**(MU - 1))
    pol = Cx([1, -T[p], b, -p**MU * T[p], p**(2*MU)])
    assert all(abs(abs(r) - C(p)**(-Integer(MU)/2)) < C(2)**-100 for r in pol.roots(multiplicities=False)), p
changed = [n for n in range(1, M + 1) if an_new[n-1] != an_old[n-1]]
nl = [new[p] for p in sq]
print('T2 (inverted):', T2)
print('T(p^2) old:', [old[p] for p in sq])
print('T(p^2) new:', nl)
print('a_n changed at %d of %d n: %s; a_9 %d -> %d, a_25 %d -> %d'
      % (len(changed), M, changed, an_old[8], an_new[8], an_old[24], an_new[24]))
if '--apply' in sys.argv:
    db.smf_newforms.update({'label': L}, {'trace_lambda_p_square': nl, 'traces': an_new}, restat=False)
    db.smf_hecke_nf.update({'label': L}, {'lambda_p_square': [[x] for x in nl], 'an': [[x] for x in an_new]}, restat=False)
    db.smf_newspaces.update({'label': SP}, {'cusp_G_lambda_p_square': nl}, restat=False)
    nf2, h2 = db.smf_newforms.lookup(L), db.smf_hecke_nf.lookup(L)
    sp2 = db.smf_newspaces.lookup(SP, ['cusp_G_lambda_p_square'])
    assert nf2['trace_lambda_p_square'] == nl and nf2['traces'] == an_new
    assert h2['lambda_p_square'] == [[x] for x in nl] and h2['an'] == [[x] for x in an_new]
    assert sp2['cusp_G_lambda_p_square'] == nl
    print('applied, read back ok')
