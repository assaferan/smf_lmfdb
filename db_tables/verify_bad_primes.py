"""
Independent check of a dry run of add_bad_primes.py (PREFIX_values.json), against the db.
Does not import add_bad_primes, bad_primes_ring or bad_primes_local.

For each smf_hecke_nf row in the output (with its smf_newforms row):
  basis      the rows start with 1, the inverse is the inverse, every coordinate is an integer;
             if the ring did not change, the basis is the stored one;
  good       at every n prime to N and every good p, lambda_p, lambda_p_square, lambda_p_square_i
             and a_n are the stored field elements;
  bad        at p | N, for the weight (k, j), a = k - 3 + j/2, w = 2k + j - 3 (a = 0, w = 3 for (3, 0)):
                       p || N:  mu = -eps p^a lambda - p^(2a) - p^(2a+1), a_p = c = lambda + eps p^a,
                                 a_{p^2}, a_{p^3} from the Euler factor
                                 (1 - (c + eps p^(a+1)) x + p^w x^2)(1 + eps p^(a+1) x);
                       p^2 | N: a_p = lambda, a_{p^2} = lambda^2 - p^w - p mu;
                       Ramanujan at p || N: |c + eps p^(a+1)| <= 2 p^(w/2) (exact, 1000-bit roots);
             T(p^2, i) stays NULL at p | N;
  multiplicativity  a_{mn} = a_m a_n for coprime m, n (all n stored);
  traces     smf_newforms traces[n] = Tr(a_n), trace_lambda_p[i] = Tr(lambda_p), trace_lambda_p_square;
  ring       the coordinates of the lambda_j (j = p, p^2, all p) generate the new ring: the Z-module
             of the closure of their coordinates has index 1; nbound is minimal for that, with the
             good p^2 generator being lambda_p_square... (only checked to be an upper bound: the
             converter uses omf5's T(p,2), which smf_hecke_nf does not store);
             for degree <= 8: [O_K : H] = the stored index (disc(H) = index^2 disc(K)).
For each smf_newforms row: the values at the good primes / coprime n are the stored ones.
For each smf_newspaces row: at p | N, cusp_G/P_lambda_p(_square) is the sum over the forms of the
type (new values where given, stored otherwise), or NULL when one of them has none.

Usage (from the lmfdb root):  sage -python -m smf_lmfdb.db_tables.verify_bad_primes PREFIX
"""
import json, sys
from collections import Counter
from sage.all import (Matrix, QQ, ZZ, vector, lcm, prime_range, gcd, factor, NumberField, PolynomialRing,
                      ComplexField, RealField, prod)
from lmfdb import db

PS = prime_range(1000)
CF, R = ComplexField(1000), RealField(1000)
Qx = PolynomialRing(QQ, 'x')
prefix = sys.argv[1]
out = json.load(open(prefix + "_values.json"))
stats, bad = Counter(), []


def fail(*a):
    bad.append(a)
    stats['FAIL ' + a[0]] += 1


def zspan(vs):
    M = Matrix(QQ, vs)
    den = lcm([c.denominator() for c in M.list()] + [1])
    return [r / den for r in (M * den).change_ring(ZZ).echelon_form().rows() if r]


def basis(nums, dens, K):
    if K is None:
        return [QQ(n[0]) / d for n, d in zip(nums, dens)]
    b = K.gen()
    return [sum(QQ(c) * b**i for i, c in enumerate(n)) / d for n, d in zip(nums, dens)]


def el(v, B):
    return None if v in (None, 'NULL') else sum(QQ(c) * x for c, x in zip(v, B))


labels = list(out['hecke_nf'])
for k, hl in enumerate(labels):
    new = out['hecke_nf'][hl]
    h = db.smf_hecke_nf.lookup(hl)
    f = db.smf_newforms.lucky({'hecke_orbit_code': h['hecke_orbit_code']})
    nfn = out['newforms'].get(f['label'])
    if nfn is None:
        fail('no newforms row', hl); continue
    D = f['level']; d = f['dim']
    a = f['weight'][0] - 3 + f['weight'][1] // 2; w = 2 * a + 3
    poly = Qx([QQ(c) for c in h['field_poly']])
    K = None if poly.degree() == 1 else NumberField(poly, 'b')
    Bo = basis(h['hecke_ring_numerators'], h['hecke_ring_denominators'], K)
    Bn = basis(new.get('hecke_ring_numerators', h['hecke_ring_numerators']), new.get('hecke_ring_denominators', h['hecke_ring_denominators']), K)
    changed = new.get('hecke_ring_numerators', h['hecke_ring_numerators']) != h['hecke_ring_numerators'] or new.get('hecke_ring_denominators', h['hecke_ring_denominators']) != h['hecke_ring_denominators']
    if K is not None:
        Mn = Matrix(QQ, [list(x) for x in Bn])
        if Mn.row(0) != vector(QQ, [1] + [0] * (d - 1)): fail('basis does not start with 1', hl)
        if changed:
            Mi = Matrix(QQ, [[QQ(c) / den for c in n] for n, den in zip(new['hecke_ring_inverse_numerators'], new['hecke_ring_inverse_denominators'])])
            if Mn * Mi != 1: fail('inverse basis', hl)
            # the new ring contains the old
            Mo = Matrix(QQ, [list(x) for x in Bo])
            if not all(c in ZZ for c in (Mo * Mn.inverse()).list()): fail('new ring does not contain the old', hl)
    eps = dict(f['atkin_lehner_eigenvals'])
    badp = [p for p in PS if D % p == 0]
    for col in ('lambda_p', 'lambda_p_square', 'lambda_p_square_0', 'lambda_p_square_1', 'lambda_p_square_2', 'an'):
        if new.get(col) is None: continue
        for i, v in enumerate(new[col]):
            if v not in (None, 'NULL') and not all(isinstance(c, int) for c in v): fail('non-integer coordinate', hl, col, i)
    E = lambda col, i: el(new[col][i], Bn)
    O = lambda col, i: el(h[col][i], Bo)
    # good primes / coprime n: the stored elements
    for col in ('lambda_p', 'lambda_p_square', 'lambda_p_square_0', 'lambda_p_square_1', 'lambda_p_square_2'):
        if h.get(col) is None or col not in new: continue
        for i in range(len(h[col])):
            p = PS[i]
            if D % p == 0:
                if col.startswith('lambda_p_square_') and new[col][i] not in (None, 'NULL'): fail('T(p^2,i) not NULL at bad p', hl, col, p)
                continue
            if E(col, i) != O(col, i): fail('good value changed', hl, col, p)
    an = [E('an', n - 1) for n in range(1, len(new['an']) + 1)]
    for n in range(1, len(an) + 1):
        if gcd(n, D) == 1 and an[n - 1] != O('an', n - 1): fail('coprime a_n changed', hl, n)
    if an[0] != 1: fail('a_1', hl)
    # bad primes
    for p in badp:
        i = PS.index(p)
        if i >= len(new['lambda_p']):
            stats['bad p beyond maxp'] += 1; continue
        lam = E('lambda_p', i)
        mu = E('lambda_p_square', i) if i < len(new['lambda_p_square'] or []) else None
        if lam is not None and mu is None and (D // p) % p:
            mu = -eps[p] * p**a * lam - p**(2 * a) - p**(2 * a + 1)   # beyond maxp_square: implied by lambda at p || N
            stats['bad p beyond maxp_square (p||N)'] += 1
        elif mu is None and lam is not None:
            stats['bad p beyond maxp_square (p^2|N): Euler factor not checked'] += 1
            if p <= len(an) and an[p - 1] != lam: fail('a_p != lambda (p^2|N)', hl, p)
            continue
        if lam is None or mu is None: fail('bad value missing', hl, p); continue
        e_ = eps[p]
        if (D // p) % p:                      # p || N
            c = lam + e_ * p**a
            if mu != -e_ * p**a * lam - p**(2 * a) - p**(2 * a + 1): fail('mu relation (p||N)', hl, p)
            ap, ap2 = c, None
            # 1/((1 - (c + e p^(a+1)) x + p^w x^2)(1 + e p^(a+1) x)): a_p = c, ... computed by recursion
            b = [1, -c, p**w - (c + e_ * p**(a + 1)) * e_ * p**(a + 1), p**w * e_ * p**(a + 1)]
            a_ = [1]
            for kk in range(1, 4):
                a_.append(-sum(b[j] * a_[kk - j] for j in range(1, min(kk, 3) + 1)))
            y = c + e_ * p**(a + 1)
            vals = [CF(y)] if K is None else [CF(r) for r, _ in y.minpoly().roots(CF)] if y not in QQ else [CF(y)]
            if f['aut_rep_type'] == 'G':
                if max(abs(v) for v in vals) > 2 * R(p) ** (R(w) / 2) * (1 + R(10)**-50): fail('Ramanujan at p||N (G)', hl, p)
                else: stats['Ramanujan ok (G, p||N)'] += 1
        else:
            b = [1, -lam, p**w + p * mu, 0]
            a_ = [1]
            for kk in range(1, 4):
                a_.append(-sum(b[j] * a_[kk - j] for j in range(1, min(kk, 3) + 1)))
        for kk in range(1, 4):
            if p**kk <= len(an) and an[p**kk - 1] != a_[kk]: fail('a_{p^k} from the Euler factor', hl, p, kk)
    # multiplicativity
    for n in range(2, len(an) + 1):
        fa = factor(n)
        if len(fa) > 1:
            q, e2 = fa[0]; m = q**e2
            if an[n - 1] != an[m - 1] * an[n // m - 1]: fail('multiplicativity', hl, n); break
    # traces
    tr = (lambda x: d * QQ(x)) if K is None else (lambda x: x.trace() if x not in QQ else d * QQ(x))
    for n in range(1, min(len(an), len(nfn['traces'])) + 1):
        if nfn['traces'][n - 1] is not None and an[n - 1] is not None and nfn['traces'][n - 1] != tr(an[n - 1]): fail('trace of a_n', hl, n); break
    for i in range(min(len(new['lambda_p']), len(nfn['trace_lambda_p']))):
        if nfn['trace_lambda_p'][i] is not None and new['lambda_p'][i] not in (None, 'NULL') and nfn['trace_lambda_p'][i] != tr(E('lambda_p', i)): fail('trace lambda_p', hl, PS[i]); break
    # good-prime traces in smf_newforms unchanged
    for i, p in enumerate(PS[:len(f['trace_lambda_p'])]):
        if D % p and nfn['trace_lambda_p'][i] != f['trace_lambda_p'][i]: fail('stored good trace changed', f['label'], p); break
    # ring generated by the eigenvalues, nbound an upper bound, index
    if K is not None and changed:
        Mi = Mn.inverse()
        co = lambda x: vector(QQ, list(K(x))) * Mi
        gens = sorted([(PS[i], E('lambda_p', i)) for i in range(len(new['lambda_p']))] +
                      [(PS[i] ** 2, E('lambda_p_square', i)) for i in range(len(new['lambda_p_square'] or []))], key=lambda t: t[0])
        nb = nfn['hecke_ring_generator_nbound']
        def gen_ok(bound):
            M = zspan([co(K(1))] + [co(x) for j, x in gens if x is not None and j <= bound])
            while True:
                Bm = [sum(c * x for c, x in zip(r, Bn)) for r in M]
                M2 = zspan([co(x * y) for ii, x in enumerate(Bm) for y in Bm[ii:]] + M)
                if M2 == M: break
                M = M2
            return len(M) == d and abs(Matrix(QQ, M).det()) == 1
        if gen_ok(nb): stats['nbound ok (with T(p^2))'] += 1
        elif gen_ok(10**9): stats['generated with all lambda_j, not by j <= nbound with T(p^2) (nbound uses T(p,2))'] += 1
        else: fail('the eigenvalues do not generate the new ring', hl)
        idx = prod([ZZ(p) ** e2 for p, e2 in nfn['hecke_ring_index_factorization']]) if nfn['hecke_ring_index_factorization'] else ZZ(1)
        if d <= 8:
            discH = Matrix(QQ, [[(x * y).trace() for y in Bn] for x in Bn]).det()
            if discH != idx**2 * K.discriminant(): fail('index vs disc(K)', hl, idx)
            else: stats['index checked against disc(K)'] += 1
    stats['rows checked'] += 1
    stats['ring changed' if changed else 'ring unchanged'] += 1
    if k % 500 == 0: print(k, dict(stats), flush=True)

# newspace sums at p | N
for sp, row in out.get('newspaces', {}).items():
    D = int(sp.split('.')[2])
    fs = list(db.smf_newforms.search({'space_label': sp, 'aut_rep_type': {'$in': ['G', 'P']}}, ['label', 'aut_rep_type', 'trace_lambda_p', 'trace_lambda_p_square']))
    for t in ('G', 'P'):
        for col, fcol in (('cusp_%s_lambda_p' % t, 'trace_lambda_p'), ('cusp_%s_lambda_p_square' % t, 'trace_lambda_p_square')):
            for i, x in enumerate(row.get(col) or []):
                if D % PS[i]: continue
                vals = []
                for f in fs:
                    if f['aut_rep_type'] != t: continue
                    src = out['newforms'][f['label']][fcol] if f['label'] in out['newforms'] else (f[fcol] or [])
                    vals.append(src[i] if i < len(src) else None)
                want = sum(vals) if all(v is not None for v in vals) else None
                if x != want: fail('newspace sum', sp, col, PS[i], x, want)
                else: stats['newspace sums ok'] += 1

print(dict(stats))
for b_ in bad[:30]: print("  ", b_)
json.dump([list(map(str, b_)) for b_ in bad], open(prefix + "_verify_fail.json", "w"))
