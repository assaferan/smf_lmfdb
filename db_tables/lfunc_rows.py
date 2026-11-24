"""
Rows of the LMFDB L-function tables (lfunc_lfunctions, lfunc_instances, lfunc_plot)
for the spinor L-functions of Siegel paramodular newforms with rational Hecke
eigenvalues, from a certified computation with Platt-Costa's glfunc
(edgarcosta/lfunctions) through the driver smf_driver.c (kept next to this file;
on lava it is built as ~/lfun/smf_driver.exe).

Two steps:

  1. `input`: build the glfunc input line from omf5's per-prime output,
         label:4:N:w:[mus]:[[1,b1,b2,b3,b4], ...]    (one factor per prime, in order)
     and a sidecar <out>.factors.json with the exact local factors (null for an
     unknown coefficient).
  2. `rows`: from the driver's JSON and the sidecar, write the rows, as JSON and as
     tab-separated upload files (header: column names, column types, blank line).

Conventions, for weight (k, j), w = 2k + j - 3 (the motivic weight):
  good p:  1 - l x + p^(w-2) (b + p^2 + 1) x^2 - p^w l x^3 + p^(2w) x^4,
           l = T(p) and b = T(p,2) eigenvalues from omf5 (b only needed when
           p^2 <= the coefficient bound; above it the x^2 coefficient does not affect
           a_n for n <= bound, and glfunc is given 0 there)
  p || N:  (1 - (c + eps p) x + p^3 x^2)(1 + eps p x)          [bad_primes_local, case 1]
  p^2 | N: 1 - lam x + (p^3 + p mu) x^2                        [bad_primes_local, case 2]
  gamma:   Gamma_C(s) Gamma_C(s - (k - 2)), i.e. Gamma_R shifts mus = [0, 1, 2-k, 3-k]
           (arithmetic normalisation), nu_real_doubled = [j + 1, w] (analytic).

Usage (on lava; the driver runs on lava):
  python3 lfunc_rows.py input --dir ~/omf5b/timing/ev61 --level 61 --weight 3 0 \
      --field-poly '2  7 1' --bad 61:1:144:-1 --bound 1084 --name 61.a --out in61_full.txt
  ~/lfun/smf_driver.exe in61_full.txt e61.json --self-dual
  python3 lfunc_rows.py rows --driver e61.json --factors in61_full.txt.factors.json \
      --smf-label 2.K.61.3.0.a.a --weight 3 0 --aut-rep-type G --outdir rows61
"""
import argparse
import ast
import json
import math
import os
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# small arithmetic helpers (pure python, so the script runs without sage)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def primes_upto(n):
    return [p for p in range(2, n + 1) if is_prime(p)]


def prime_factors(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d += 1
    if n > 1:
        out.append((n, 1))
    return out


def trim(f):
    """Drop trailing zero coefficients (LMFDB stores bad factors that way)."""
    f = list(f)
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def inverse_series(f, n):
    """Coefficients 0..n-1 of 1/f for f with f[0] = 1 (integers)."""
    out = [1] + [0] * (n - 1)
    for i in range(1, n):
        out[i] = -sum(f[m] * out[i - m] for m in range(1, min(i, len(f) - 1) + 1))
    return out


def dirichlet_an(factors, nmax):
    """a_1..a_nmax (index 0 unused) from {p: ascending local factor}."""
    a = [0, 1] + [None] * (nmax - 1)
    pp = {}
    for p in primes_upto(nmax):
        r = 1
        while p ** (r + 1) <= nmax:
            r += 1
        f = factors[p]
        if any(c is None for c in f[:r + 1]):
            raise ValueError("a_%d needs the unknown coefficients of the factor at %d" % (p ** r, p))
        pp[p] = inverse_series(f, r + 1)
    for n in range(2, nmax + 1):
        val = 1
        for p, e in prime_factors(n):
            val *= pp[p][e]
        a[n] = val
    return a


def digamma(x):
    """psi(x) for real x > 0, to double precision."""
    r = 0.0
    while x < 12:
        r -= 1 / x
        x += 1
    f = 1 / (x * x)
    return r + math.log(x) - 0.5 / x - f * (1 / 12 - f * (1 / 120 - f * (1 / 252 - f * (1 / 240 - f / 132))))


def analytic_conductor(N, mu_real, nu_real_doubled):
    """N * prod_R exp(psi((1/2 + mu)/2))/pi * prod_C (exp(psi(1/2 + nu))/(2 pi))^2, at the centre,
    with nu = nu_real_doubled/2 (analytic normalisation); this reproduces the LMFDB values."""
    c = float(N)
    for mu in mu_real:
        c *= math.exp(digamma((0.5 + mu) / 2)) / math.pi
    for nu2 in nu_real_doubled:
        c *= (math.exp(digamma(0.5 + nu2 / 2)) / (2 * math.pi)) ** 2
    return c


# ---------------------------------------------------------------------------
# local factors


def good_factor(p, l, b, k, j):
    """Spinor factor at a good prime, b = T(p,2) eigenvalue or None if unknown."""
    w = 2 * k + j - 3
    b2 = None if b is None else p ** (w - 2) * (b + p ** 2 + 1)
    return [1, -l, b2, -p ** w * l, p ** (2 * w)]


def bad_factor(case, p, eps=None, c=None, lam=None, mu=None):
    """Factor at p | N, as in bad_primes_local.local_from_omf5 (weight (3,0))."""
    if case == 1:
        # (1 - (c + eps p) x + p^3 x^2)(1 + eps p x)
        return trim([1, -c, p ** 3 - (c + eps * p) * eps * p, p ** 4 * eps])
    if case == 2:
        mu = Fraction(mu)
        b2 = p ** 3 + p * mu
        assert b2.denominator == 1
        return trim([1, -lam, int(b2)])
    raise ValueError(case)


# ---------------------------------------------------------------------------
# step 1: glfunc input from omf5 output


def read_omf5(path, field_poly):
    """(lambda_p, lambda_p_square or None) of the rational form with this field_poly
    in an omf5 -format=GG output file (a python-literal dict {disc: [forms]})."""
    with open(path) as F:
        txt = F.read()
    txt = txt[txt.index('{'):txt.rindex('}') + 1]
    data = ast.literal_eval(txt)
    for forms in data.values():
        for f in forms:
            if f['field_poly'] == field_poly:
                lp = int(f['lambda_p'][0])
                lp2 = int(f['lambda_p_square'][0]) if f['lambda_p_square'] else None
                return lp, lp2
    raise ValueError("no form with field_poly %r in %s" % (field_poly, path))


def omf5_file(directory, p):
    """The newest omf5 output for p: q<p>.txt (full run) before p<p>.txt (pilot)."""
    for name in ['q%d.txt' % p, 'p%d.txt' % p]:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    return None


def parse_bad(spec):
    """'p:1:c:eps' (p || N) or 'p:2:lam:mu' (p^2 | N) or 'p:f0,f1,...' (explicit factor)."""
    parts = spec.split(':')
    p = int(parts[0])
    if len(parts) == 2:
        return p, trim([int(c) for c in parts[1].split(',')])
    case = int(parts[1])
    if case == 1:
        return p, bad_factor(1, p, eps=int(parts[3]), c=int(parts[2]))
    return p, bad_factor(2, p, lam=int(parts[2]), mu=Fraction(parts[3]))


def gamma_data(k, j):
    w = 2 * k + j - 3
    mus = sorted([0, 1, 2 - k, 3 - k])           # Gamma_R shifts, arithmetic normalisation
    gamma_factors = [[], sorted([0, 2 - k])]      # Gamma_C shifts, arithmetic normalisation
    nu_real_doubled = sorted([j + 1, w])          # 2 * Gamma_C shifts, analytic normalisation
    return w, mus, gamma_factors, nu_real_doubled


def spec_line(name, N, k, j, factors, nprimes=None):
    w, mus, _, _ = gamma_data(k, j)
    ps = sorted(factors)
    if nprimes is not None:
        ps = ps[:nprimes]
    assert ps == primes_upto(ps[-1]), "missing primes"
    rows = []
    for p in ps:
        f = [0 if c is None else c for c in factors[p]]
        f = (f + [0] * 5)[:5]
        rows.append('[%s]' % ','.join(map(str, f)))
    return '%s:4:%d:%d:[%s]:[%s]' % (name, N, w, ','.join(map(str, mus)), ','.join(rows))


def cmd_input(args):
    k, j = args.weight
    N = args.level
    bad = dict(parse_bad(s) for s in args.bad)
    assert sorted(bad) == [p for p, _ in prime_factors(N)], "need a --bad for each p | N"
    factors, missing = {}, []
    for p in primes_upto(args.bound):
        if p in bad:
            factors[p] = bad[p]
            continue
        path = omf5_file(args.dir, p)
        if path is None:
            missing.append(p)
            continue
        l, b = read_omf5(path, args.field_poly)
        if p * p <= args.bound and b is None:
            raise ValueError("T(p,2) needed at p = %d (p^2 <= %d) but missing in %s" % (p, args.bound, path))
        factors[p] = good_factor(p, l, b, k, j)
    if missing:
        sys.stderr.write("missing omf5 output for %d primes: %s\n" % (len(missing), missing[:20]))
        if not args.allow_missing:
            sys.exit(1)
        # keep the initial segment of consecutive primes
        first = min(missing)
        factors = {p: f for p, f in factors.items() if p < first}
    line = spec_line(args.name, N, k, j, factors)
    with open(args.out, 'w') as F:
        F.write(line + '\n')
    with open(args.out + '.factors.json', 'w') as F:
        json.dump({'level': N, 'weight': [k, j], 'bound': args.bound,
                   'factors': {str(p): f for p, f in sorted(factors.items())}}, F)
    sys.stderr.write("wrote %s (%d primes, up to %d)\n" % (args.out, len(factors), max(factors)))


def factors_from_spec(path, bound, bad_primes):
    """Fallback when no sidecar exists: read the factors from a glfunc input line,
    marking the x^2 coefficient of good p with p^2 > bound as unknown."""
    with open(path) as F:
        line = F.readline().strip()
    efs = json.loads(line.split(':')[5])
    out = {}
    for p, f in zip(primes_upto(10 ** 5), efs):
        f = list(f)
        if p in bad_primes:
            f = trim(f)
        elif p * p > bound:
            f[2] = None
        out[p] = f
    return out


# ---------------------------------------------------------------------------
# step 2: rows

def ball_digits(mid, rad, maxdigits=30):
    """mid (a decimal string) truncated to the digits its radius supports, at most maxdigits."""
    if rad <= 0:
        good = maxdigits
    else:
        good = max(1, min(maxdigits, int(-math.log10(rad)) - 1))
    # good = digits after the decimal point
    if '.' not in mid:
        return mid
    ip, fp = mid.split('.')
    if 'e' in fp or 'E' in fp:
        return mid
    return ip + '.' + fp[:good]


def sig_digits(mid, n=30):
    """mid with n significant digits (as in the CMF z1, z2, z3 columns)."""
    neg = mid.startswith('-')
    s = mid.lstrip('-')
    ip, _, fp = s.partition('.')
    digits = len(ip.lstrip('0'))
    keep = max(0, n - digits)
    return ('-' if neg else '') + ip + ('.' + fp[:keep] if keep and fp else '')


class Samples:
    """Z(t) on the grid t_i = i h, i = 0..n-1, as returned by Lfunc_plot_data."""

    def __init__(self, spacing, values):
        self.h = spacing
        self.y = values
        self.T = spacing * (len(values) - 1)

    def _stencil(self, x):
        n = len(self.y)
        i = int(math.floor(x / self.h)) - 1
        i = min(max(i, 0), n - 4)
        return [(self.h * m, self.y[m]) for m in range(i, i + 4)]

    def value(self, x):
        pts = self._stencil(x)
        s = 0.0
        for a, (xa, ya) in enumerate(pts):
            L = 1.0
            for b, (xb, _) in enumerate(pts):
                if a != b:
                    L *= (x - xb) / (xa - xb)
            s += ya * L
        return s

    def deriv(self, x):
        pts = self._stencil(x)
        s = 0.0
        for a, (xa, ya) in enumerate(pts):
            dL = 0.0
            for c, (xc, _) in enumerate(pts):
                if c == a:
                    continue
                term = 1.0 / (xa - xc)
                for b, (xb, _) in enumerate(pts):
                    if b != a and b != c:
                        term *= (x - xb) / (xa - xb)
                dL += term
            s += ya * dL
        return s


def hermite(control, x):
    import bisect
    xs = [c[0] for c in control]
    i = min(max(bisect.bisect_right(xs, x), 1), len(xs) - 1)
    (x0, y0, d0), (x1, y1, d1) = control[i - 1], control[i]
    h = x1 - x0
    t = (x - x0) / h if h else 1.0
    return ((2 * t ** 3 - 3 * t ** 2 + 1) * y0 + (t ** 3 - 2 * t ** 2 + t) * h * d0
            + (3 * t ** 2 - 2 * t ** 3) * y1 + (t ** 3 - t ** 2) * h * d1)


def plot_row(S, zeros, sign, rank, tol_rel=1e-3):
    """lfunc_plot row (PR LMFDB/lmfdb#7272) of a self-dual L-function: control points
    (x, Z(x), Z'(x)) for cubic Hermite interpolation, x in [-T, T], Z > 0 just right of 0."""
    flip = 1.0
    # sign of Z just right of 0 (rank r: Z ~ c t^r)
    probe = next((y for y in S.y[1:] if y != 0), 1.0)
    if probe < 0:
        flip = -1.0
    y = [flip * v for v in S.y]
    S = Samples(S.h, y)
    zs = [z for z in zeros if 0 < z < S.T]
    pos_zeros = [(z, S.deriv(z)) for z in zs]
    # extrema between consecutive zeros (and before the first, after the last)
    bounds = [0.0] + zs + [S.T]
    extrema = []
    for a, b in zip(bounds, bounds[1:]):
        i0 = int(math.ceil(a / S.h)) + 1
        i1 = int(math.floor(b / S.h)) - 1
        if i1 - i0 < 2:
            continue
        im = max(range(i0, i1 + 1), key=lambda i: abs(y[i]))
        if im in (i0, i1):          # monotone on this stretch (e.g. near 0 or T)
            continue
        ym, y0, y1 = y[im], y[im - 1], y[im + 1]
        den = y0 - 2 * ym + y1
        dx = 0.5 * S.h * (y0 - y1) / den if den else 0.0
        xe = im * S.h + dx
        extrema.append((xe, S.value(xe)))
    control = ([(x, 0.0, d) for x, d in pos_zeros] + [(x, v, 0.0) for x, v in extrema]
               + [(0.0, y[0], 0.0 if rank % 2 == 0 else S.deriv(0.0)), (S.T, y[-1], S.deriv(S.T))])
    extra = [(0.0, y[0], control[-2][2]), (S.T, y[-1], control[-1][2])]
    tol = tol_rel * max(abs(v) for v in y)
    for _ in range(200):
        control.sort()
        errs = [(abs(hermite(control, i * S.h) - y[i]), i) for i in range(len(y))]
        err, i = max(errs)
        if err <= tol:
            break
        pt = (i * S.h, y[i], S.deriv(i * S.h))
        control.append(pt)
        extra.append(pt)
    # mirror: Z(-x) = sign * Z(x) for a self-dual L-function
    def mirror_pts(pts):
        return [(-x, sign * v, -sign * d) for x, v, d in pts if x > 0]
    extra = sorted(mirror_pts(extra) + extra)
    zeros_rows = sorted([[-x, -sign * d] for x, d in pos_zeros] + [[x, d] for x, d in pos_zeros])
    if rank % 2 == 1:
        zeros_rows.append([0.0, S.deriv(0.0)])
        zeros_rows.sort()
        extra = [e for e in extra if e[0] != 0.0]
    extrema_rows = sorted([[-x, sign * v] for x, v in extrema] + [[x, v] for x, v in extrema])
    return {
        'zeros': zeros_rows,
        'extrema': extrema_rows,
        'extra_control': [list(e) for e in extra],
        'central_value': y[0] if rank == 0 else None,
        'sign_flipped': flip < 0,
    }


LFUNC_TYPES = {
    'A10': 'numeric', 'A2': 'numeric', 'A3': 'numeric', 'A4': 'numeric', 'A5': 'numeric',
    'A6': 'numeric', 'A7': 'numeric', 'A8': 'numeric', 'A9': 'numeric', 'Lhash': 'text',
    'a10': 'jsonb', 'a2': 'jsonb', 'a3': 'jsonb', 'a4': 'jsonb', 'a5': 'jsonb', 'a6': 'jsonb',
    'a7': 'jsonb', 'a8': 'jsonb', 'a9': 'jsonb', 'accuracy': 'smallint', 'algebraic': 'boolean',
    'analytic_conductor': 'double precision', 'analytic_normalization': 'numeric',
    'bad_lfactors': 'jsonb', 'bad_primes': 'bigint[]', 'central_character': 'text',
    'coeff_info': 'jsonb', 'coefficient_field': 'text', 'conductor': 'numeric',
    'conductor_radical': 'integer', 'conjugate': 'text', 'credit': 'text', 'degree': 'smallint',
    'dirichlet_coefficients': 'jsonb', 'euler_factors': 'jsonb',
    'euler_factors_factorization': 'jsonb', 'gamma_factors': 'jsonb', 'group': 'text',
    'index': 'smallint', 'label': 'text', 'leading_term': 'text', 'load_key': 'text',
    'motivic_weight': 'smallint', 'mu_imag': 'numeric[]', 'mu_real': 'smallint[]',
    'nu_imag': 'numeric[]', 'nu_real_doubled': 'smallint[]', 'order_of_vanishing': 'smallint',
    'origin': 'text', 'plot_delta': 'numeric', 'plot_values': 'jsonb', 'positive_zeros': 'jsonb',
    'precision': 'smallint', 'prelabel': 'text', 'primitive': 'boolean', 'rational': 'boolean',
    'root_analytic_conductor': 'double precision', 'root_angle': 'double precision',
    'root_number': 'text', 'self_dual': 'boolean', 'sign_arg': 'numeric',
    'spectral_label': 'text', 'st_group': 'text', 'symmetry_type': 'text',
    'trace_hash': 'bigint', 'types': 'jsonb', 'values': 'jsonb', 'z1': 'numeric',
    'z2': 'numeric', 'z3': 'numeric',
}
INSTANCE_TYPES = {'url': 'text', 'Lhash': 'text', 'type': 'text', 'Lhash_array': 'text[]',
                  'factors': 'text[]', 'label': 'text'}
PLOT_TYPES = {'Lhash': 'text', 'label': 'text', 'dual_Lhash': 'text', 'dual_label': 'text',
              'dual_yflip': 'boolean', 'central_value': 'double precision',
              'zeros': 'double precision[]', 'extrema': 'double precision[]',
              'extra_control': 'double precision[]'}


def conductor_str(N):
    """N as in L-function labels: m'e'e when N = m^e with e > 1 maximal (13e2, 2e13), else N."""
    for e in range(int(math.log2(N)) if N > 1 else 1, 1, -1):
        m = round(N ** (1.0 / e))
        for mm in (m - 1, m, m + 1):
            if mm > 1 and mm ** e == N:
                return '%de%d' % (mm, e)
    return str(N)


def spectral(nu_real_doubled):
    """Spectral label of an algebraic L-function with only Gamma_C factors (c1e2-0, c1c3-0),
    in the normal form of lmfdb.lfunctions.main.parse_spectral."""
    from collections import Counter
    cnt = Counter(nu_real_doubled)
    g = 0
    for v in cnt.values():
        g = math.gcd(g, v)
    s = ''.join('c%d' % x for x in sorted(cnt) for _ in range(cnt[x] // g))
    return s + ('' if g == 1 else 'e%d' % g) + '-0'


def smf_url(smf_label):
    deg, fam, N, k, j, chi, orb = smf_label.split('.')
    return 'ModularForm/GSp/Q/%s/%s/%s/%s.%s/%s/%s' % (deg, fam, N, k, j, chi, orb)


def lfunction_rows(drv, factors, smf_label, k, j, aut_rep_type='G', index=0,
                   accuracy=100, max_zeros=None, plot_stride=8, load_key='SMF',
                   instance_type='SMF', st_group=None, allow_uncertified=False):
    N = drv['conductor']
    w, mus, gamma_factors, nu2 = gamma_data(k, j)
    assert drv['degree'] == 4 and drv['motivic_weight'] == w and drv['mus'] == mus, \
        "driver run does not match weight (%d,%d)" % (k, j)
    if not drv['certified'] and not allow_uncertified:
        raise ValueError("driver output not certified: %s (nmax %d of %d)"
                         % (drv['errors'], drv['nmax_used'], drv['nmax_required']))
    rank = drv['rank']
    sign = round(float(drv['sign']['re']['mid']))
    assert sign in (1, -1) and abs(float(drv['sign']['im']['mid'])) < 1e-10
    # (-1)^k times the Fricke sign; checked by the caller when known

    # zeros: those known to `accuracy` bits (the LMFDB rounds them to 2^-accuracy)
    zeros = [z for z in drv['zeros'] if z['rad'] < 2.0 ** -(accuracy + 1)]
    if max_zeros:
        zeros = zeros[:max_zeros]
    if len(zeros) < 3:
        raise ValueError("only %d zeros known to %d bits" % (len(zeros), accuracy))
    # CMF convention: round(first positive zero * 2^100); G2Q rows use trace_hash instead
    assert zeros[0] is drv['zeros'][0], "first zero not known to %d bits" % accuracy
    Lhash = zeros[0]['x2e100']

    bad_primes = [p for p, _ in prime_factors(N)]
    ps100 = primes_upto(100)
    euler = [factors[p] for p in ps100]
    an = dirichlet_an(factors, 10)
    complete = all(c is not None for f in euler for c in f)
    spectral_label = spectral(nu2)
    prelabel = '4-%s-1.1-%s' % (conductor_str(N), spectral_label)
    label = '%s-%d' % (prelabel, index)
    ac = analytic_conductor(N, [], nu2)
    taylor = drv['taylor']
    leading = ball_digits(taylor['mid'], taylor['rad'])
    efac = None
    if complete:
        efac = euler_factorizations(euler)
    S = Samples(drv['plot']['spacing'], drv['plot']['values'])
    origin = smf_url(smf_label)
    row = {
        'Lhash': Lhash,
        'label': label,
        'prelabel': prelabel,
        'index': index,
        'origin': origin,
        'load_key': load_key,
        'types': [instance_type],
        'degree': 4,
        'conductor': N,
        'conductor_radical': math.prod(bad_primes),
        'bad_primes': bad_primes,
        'central_character': '%d.1' % N,
        'motivic_weight': w,
        'analytic_normalization': str(Fraction(w, 2)) if w % 2 == 0 else '%.1f' % (w / 2),
        'gamma_factors': gamma_factors,
        'mu_real': [], 'mu_imag': [],
        'nu_real_doubled': nu2, 'nu_imag': [0, 0],
        'spectral_label': spectral_label,
        'analytic_conductor': ac,
        'root_analytic_conductor': ac ** 0.25,
        'self_dual': True, 'conjugate': None,
        'algebraic': True, 'rational': True,
        'primitive': aut_rep_type == 'G',
        'coefficient_field': '1.1.1.1',
        'group': 'GSp4',
        'symmetry_type': None,
        'st_group': st_group if st_group is not None else ('USp(4)' if aut_rep_type == 'G' else None),
        'root_number': str(sign),
        'root_angle': 0.0 if sign == 1 else 0.5,
        'sign_arg': 0.0 if sign == 1 else 0.5,
        'order_of_vanishing': rank,
        'leading_term': leading,
        'positive_zeros': [z['mid'] for z in zeros],
        'accuracy': accuracy,
        'precision': None,
        'z1': sig_digits(zeros[0]['mid']), 'z2': sig_digits(zeros[1]['mid']),
        'z3': sig_digits(zeros[2]['mid']),
        'euler_factors': euler,
        'euler_factors_factorization': efac,
        'bad_lfactors': [[p, factors[p]] for p in bad_primes],
        'plot_delta': S.h * plot_stride,
        'plot_values': S.y[::plot_stride],
        'trace_hash': None,
        # lmfdb's an_from_data cannot expand factors with unknown (null) coefficients, so
        # when T(p,2) is missing for some p < 100 we store a_1..a_100 instead
        'dirichlet_coefficients': None if complete else dirichlet_an(factors, 100)[1:],
        'coeff_info': None,
        'credit': None,
        'values': None,
    }
    for n in range(2, 11):
        row['a%d' % n] = an[n]
        row['A%d' % n] = an[n]
    assert sorted(row) == sorted(LFUNC_TYPES), set(row) ^ set(LFUNC_TYPES)
    instance = {'url': origin, 'Lhash': Lhash, 'type': instance_type, 'Lhash_array': [Lhash],
                'factors': [label], 'label': label}
    pr = plot_row(S, [float(z['mid']) for z in drv['zeros']], sign, rank)
    flipped = pr.pop('sign_flipped')
    plot = dict(pr, Lhash=Lhash, label=label, dual_Lhash=None, dual_label=None, dual_yflip=None)
    meta = {
        'certified': drv['certified'], 'errors': drv['errors'],
        'nmax_required': drv['nmax_required'], 'nmax_used': drv['nmax_used'],
        'zeros_found': len(drv['zeros']), 'zeros_stored': len(zeros),
        'max_zero_radius_stored': max(z['rad'] for z in zeros),
        'taylor': taylor, 'plot_range': S.T, 'plot_sign_flipped': flipped,
        'euler_factors_complete': complete,
    }
    return row, instance, plot, meta


def euler_factorizations(euler):
    """[[[factor, multiplicity], ...] per prime] over ZZ, as in the CMF rows (needs sage)."""
    try:
        from sage.all import PolynomialRing, ZZ
    except ImportError:
        return None
    R = PolynomialRing(ZZ, 'x')
    out = []
    for f in euler:
        F = R(f).factor()
        # normalise each factor to constant term 1, as in the CMF rows
        gs = [(g if g[0] > 0 else -g, e) for g, e in F]
        out.append([[[int(c) for c in g.list()], int(e)] for g, e in
                    sorted(gs, key=lambda ge: (ge[0].degree(), ge[0].list()))])
    return out


def encode(value, typ):
    """A value as postgres COPY text (psycodict upload files)."""
    if value is None:
        return r'\N'
    if typ == 'jsonb':
        return json.dumps(value)
    if typ == 'boolean':
        return 't' if value else 'f'
    if typ.endswith('[]'):
        def arr(v):
            if isinstance(v, list):
                return '{' + ','.join(arr(x) for x in v) + '}'
            return '"%s"' % v if isinstance(v, str) else repr(v)
        return arr(value)
    if isinstance(value, float):
        return repr(value)
    return str(value)


def write_upload(path, rows, types):
    cols = sorted(types)
    with open(path, 'w') as F:
        F.write('\t'.join(cols) + '\n')
        F.write('\t'.join(types[c] for c in cols) + '\n\n')
        for r in rows:
            F.write('\t'.join(encode(r[c], types[c]) for c in cols) + '\n')


def cmd_rows(args):
    k, j = args.weight
    with open(args.driver) as F:
        drivers = [json.loads(line) for line in F if line.strip()]
    if args.factors:
        with open(args.factors) as F:
            fj = json.load(F)
        factors = {int(p): f for p, f in fj['factors'].items()}
    else:
        drv = drivers[0]
        factors = factors_from_spec(args.spec, args.spec_bound or drv['nmax_used'],
                                    [p for p, _ in prime_factors(drv['conductor'])])
    lrows, irows, prows, metas = [], [], [], []
    for drv in drivers:
        row, inst, plot, meta = lfunction_rows(
            drv, factors, args.smf_label, k, j, aut_rep_type=args.aut_rep_type,
            index=args.index, accuracy=args.accuracy, max_zeros=args.max_zeros,
            load_key=args.load_key, allow_uncertified=args.allow_uncertified)
        if args.fricke is not None:
            assert int(row['root_number']) == (-1) ** k * args.fricke, "root number != (-1)^k * Fricke"
        lrows.append(row)
        irows.append(inst)
        prows.append(plot)
        metas.append(meta)
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, 'rows.json'), 'w') as F:
        json.dump({'meta': metas, 'lfunc_lfunctions': lrows, 'lfunc_instances': irows,
                   'lfunc_plot': prows}, F, indent=1)
    write_upload(os.path.join(args.outdir, 'lfunc_lfunctions.txt'), lrows, LFUNC_TYPES)
    write_upload(os.path.join(args.outdir, 'lfunc_instances.txt'), irows, INSTANCE_TYPES)
    write_upload(os.path.join(args.outdir, 'lfunc_plot.txt'), prows, PLOT_TYPES)
    for m, r in zip(metas, lrows):
        sys.stderr.write("%s  Lhash %s  rank %d  eps %s  L* %s  z1 %s  certified %s %s\n"
                         % (r['label'], r['Lhash'], r['order_of_vanishing'], r['root_number'],
                            r['leading_term'], r['z1'], m['certified'], m['errors']))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('input', help='glfunc input line from omf5 output files')
    a.add_argument('--dir', required=True)
    a.add_argument('--level', type=int, required=True)
    a.add_argument('--weight', type=int, nargs=2, required=True)
    a.add_argument('--field-poly', required=True, help="omf5 field_poly of the rational form, e.g. '2  7 1'")
    a.add_argument('--bad', action='append', default=[], help="p:1:c:eps | p:2:lam:mu | p:f0,f1,..")
    a.add_argument('--bound', type=int, required=True, help='coefficient bound (Lfunc_nmax)')
    a.add_argument('--name', required=True)
    a.add_argument('--out', required=True)
    a.add_argument('--allow-missing', action='store_true')
    a.set_defaults(func=cmd_input)
    r = sub.add_parser('rows', help='LMFDB rows from smf_driver output')
    r.add_argument('--driver', required=True, help='smf_driver JSON output')
    r.add_argument('--factors', help='sidecar <input>.factors.json written by `input`')
    r.add_argument('--spec', help='glfunc input line (fallback when there is no sidecar)')
    r.add_argument('--spec-bound', type=int, help='bound the spec was built for (default: nmax_used)')
    r.add_argument('--smf-label', required=True)
    r.add_argument('--weight', type=int, nargs=2, required=True)
    r.add_argument('--aut-rep-type', default='G')
    r.add_argument('--fricke', type=int, help='Fricke sign, to check the root number')
    r.add_argument('--index', type=int, default=0)
    r.add_argument('--accuracy', type=int, default=100)
    r.add_argument('--max-zeros', type=int)
    r.add_argument('--load-key', default='SMF')
    r.add_argument('--allow-uncertified', action='store_true')
    r.add_argument('--outdir', required=True)
    r.set_defaults(func=cmd_rows)
    args = ap.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
