"""
glfunc inputs for the spinor L-functions of the rational paramodular G newforms, from the
database and the 40-bit omf5 sweep (T(p) at the good primes 97 < p <= Lfunc_nmax, a -prec=100
run p0.txt per level, -square runs s<p>.txt at 11 <= p <= 97, and T(p) at the bad p > 97).

Two steps:

  1. `targets` (sage, reads the db): one JSON record per rational G newform of the weights
     given, with its stored lambda_p, T(p^2), lambda and mu at p | N, Atkin-Lehner signs and a_n.
         sage -python -m smf_lmfdb.db_tables.sweep_lfun targets --out targets.json
  2. `inputs` (python3, no sage or db): per newform, the glfunc input line and the sidecar
     <line>.factors.json, as `lfunc_rows input` writes them (same conventions, same functions).
         python3 sweep_lfun.py inputs --targets targets.json --sweep ~/sweep40/out \
             --nmax ~/lfun/nmax_spin.exe --outdir inputs [LABEL ...]

  3. on lava, the driver and `lfunc_rows rows` for every input, at 40 bits (the sweep's bound):
         bash sweep_lfun_rows.sh inputs rows        # rows/<name>/rows/rows.json, rows/summary.txt
  4. load_lfunctions.py (dry run, --apply, --readback) on the rows/<name>/rows directories.
Steps 1-3 reproduce the loaded rows of 2.K.19.3.2.a.a and 2.K.31.4.0.a.a byte for byte.

Sources of the local factor at p <= nmax (w = 2k + j - 3, a = (w - 3) / 2):
  * p | N, p < 100: the stored lambda, mu and Atkin-Lehner sign eps (bad_primes_local's dictionary):
      p || N:  c = lambda + eps p^a, (1 - (c + eps p^(a+1)) x + p^w x^2)(1 + eps p^(a+1) x);
      p^2 | N: 1 - lambda x + (p^w + p mu) x^2, mu = T_{p^2} / 2 from s<p>.txt when not stored;
  * p | N, p > 97: omf5's printed c in <dir>/p<p>.txt (p || N only);
  * good p < 100: the stored lambda_p and, when stored, T(p^2), with the x^2 coefficient
    lambda_p^2 - T(p^2) - p^(w-1); otherwise T(p,2) from s<p>.txt; otherwise unknown;
  * good p > 97: T(p) from <dir>/p<p>.txt, T(p,2) from s<p>.txt if there is one.
The sweep form is the unique one of degree 1 and type G in p0.txt whose lambda_p at the good
p < 100 are the stored ones; it is the same (cond, pos) in every run of the level (same binary,
same genus; the field_poly is checked).

Checks, any failure skips the form: stored T(p^2) = the s<p>.txt value wherever both exist;
mu = -eps p^a c - p^(2a+1) at p || N; |lambda_p| <= 4 p^(w/2) at the good p; the a_n from the
factors equal the stored a_n at every n they determine; all primes up to nmax present.
"""
import argparse
import ast
import json
import os
import subprocess
import sys
from fractions import Fraction

try:
    from smf_lmfdb.db_tables.lfunc_rows import (bad_factor, good_factor, spec_line, primes_upto,
                                                prime_factors, dirichlet_an)
except ImportError:
    from lfunc_rows import bad_factor, good_factor, spec_line, primes_upto, prime_factors, dirichlet_an

OMF5_WEIGHT = {(3, 0): '0', (4, 0): '1', (3, 2): '1,1'}
PRIMES100 = primes_upto(100)


# ---------------------------------------------------------------------------
# step 1: targets from the db

def cmd_targets(args):
    from lmfdb import db
    out = []
    for wt in args.weights:
        k, j = map(int, wt.split(','))
        q = {'family': 'K', 'weight': [k, j], 'dim': 1, 'aut_rep_type': 'G'}
        if (k, j) == (3, 0) and args.max_level_30:
            q['level'] = {'$lt': args.max_level_30}
        for f in db.smf_newforms.search(q, ['label', 'level', 'weight', 'trace_lambda_p', 'trace_lambda_p_square',
                                            'atkin_lehner_eigenvals', 'fricke_eigenval', 'traces']):
            lp = dict(zip(PRIMES100, f['trace_lambda_p']))
            lp2 = dict(zip(PRIMES100, f.get('trace_lambda_p_square') or []))
            out.append({'label': f['label'], 'level': f['level'], 'weight': f['weight'],
                        'eps': {str(p): e for p, e in f['atkin_lehner_eigenvals']},
                        'fricke': f.get('fricke_eigenval'),
                        'lambda_p': {str(p): v for p, v in lp.items() if v is not None},
                        'lambda_p_square': {str(p): v for p, v in lp2.items() if v is not None},
                        'traces': f.get('traces')})
    out.sort(key=lambda r: (r['weight'], r['level'], r['label']))
    with open(args.out, 'w') as F:
        json.dump(out, F, default=int)          # the db returns sage Integers
    sys.stderr.write("wrote %d targets to %s\n" % (len(out), args.out))


# ---------------------------------------------------------------------------
# step 2: inputs

def read_run(path):
    """{(cond, pos): form} of an omf5 output file (the python-literal dict at its end)."""
    with open(path) as F:
        txt = F.read()
    txt = txt[txt.index('{'):txt.rindex('}') + 1]
    data = ast.literal_eval(txt)
    return {(cond, pos): f for cond, L in data.items() for pos, f in enumerate(L)}


def degree(field_poly):
    return int(field_poly.split()[0]) - 1


class Skip(Exception):
    pass


def find_dir(roots, W, N):
    for r in roots:
        d = os.path.join(r, W, str(N))
        if os.path.exists(os.path.join(d, 'p0.txt')):
            return d
    return None


def run_file(d, prefix, p):
    path = os.path.join(d, '%s%d.txt' % (prefix, p))
    return path if os.path.exists(path) else None


def build(t, roots, nmax_exe, stats):
    N = t['level']
    k, j = t['weight']
    w = 2 * k + j - 3
    a = (w - 3) // 2
    W = OMF5_WEIGHT[(k, j)]
    bound = int(subprocess.check_output([nmax_exe, str(N), str(k), str(j), '40']).split()[-1])
    d = find_dir(roots, W, N)
    if d is None:
        raise Skip('no p0.txt')
    lp = {int(p): v for p, v in t['lambda_p'].items()}
    lp2 = {int(p): v for p, v in t['lambda_p_square'].items()}
    eps = {int(p): e for p, e in t['eps'].items()}
    bad = [p for p, _ in prime_factors(N)]
    good100 = [p for p in PRIMES100 if N % p]
    # the sweep form
    run0 = read_run(os.path.join(d, 'p0.txt'))
    hits = [key for key, f in run0.items() if degree(f['field_poly']) == 1 and f['aut_rep_type'] == 'G'
            and [int(x) for x in f['lambda_p'][:len(good100)]] == [lp[p] for p in good100]]
    if len(hits) != 1:
        raise Skip('%d matches in p0.txt' % len(hits))
    key = hits[0]
    fpoly = run0[key]['field_poly']

    def value(path):
        f = read_run(path).get(key)
        if f is None or f['field_poly'] != fpoly:
            raise Skip('form %s not at the same place in %s' % (key, path))
        return (int(f['lambda_p'][0]),
                int(f['lambda_p_square'][0]) if f.get('lambda_p_square') else None)

    factors, missing = {}, []
    for p in primes_upto(bound):
        if p in bad:
            if p < 100:
                lam, mu = lp[p], lp2.get(p)
                if (N // p) % p:
                    c = lam + eps[p] * p ** a
                    if mu is not None and mu != -eps[p] * p ** a * c - p ** (2 * a + 1):
                        raise Skip('mu at %d is not -eps p^a c - p^(2a+1)' % p)
                    factors[p] = bad_factor(1, p, eps=eps[p], c=c, w=w)
                else:
                    if mu is None:
                        # not stored (p beyond the level's T(p^2) range): mu = T_{p^2} / 2 from s<p>.txt
                        sq = run_file(d, 's', p)
                        if sq is None:
                            raise Skip('no mu at %d (p^2 | N) and no s%d.txt' % (p, p))
                        tp, tp2 = value(sq)
                        if tp != lam or tp2 is None:
                            raise Skip('s%d.txt: T_p %d, stored lambda %d' % (p, tp, lam))
                        mu = Fraction(tp2, 2)
                    factors[p] = bad_factor(2, p, lam=lam, mu=mu, w=w)
            else:
                if (N // p) % p == 0:
                    raise Skip('p^2 | N at %d > 97' % p)
                path = run_file(d, 'p', p)
                if path is None:
                    missing.append(p)
                    continue
                factors[p] = bad_factor(1, p, eps=eps[p], c=value(path)[0], w=w)
            continue
        sq = run_file(d, 's', p)
        T2 = value(sq)[1] if sq else None
        if p < 100:
            l = lp[p]
            b2 = l ** 2 - lp2[p] - p ** (w - 1) if p in lp2 else None
            if T2 is not None:
                b2s = p * T2 + p ** (w - 2) * (p ** 2 + 1)
                if b2 is not None and b2 != b2s:
                    raise Skip('stored T(p^2) at %d disagrees with s%d.txt' % (p, p))
                b2 = b2s
                stats['T(p,2) checked against the db'] += p in lp2
        else:
            path = run_file(d, 'p', p)
            if path is None:
                missing.append(p)
                continue
            l = value(path)[0]
            b2 = None if T2 is None else p * T2 + p ** (w - 2) * (p ** 2 + 1)
        if l * l > 16 * p ** w:
            raise Skip('|lambda_%d| > 4 p^(w/2)' % p)
        f = good_factor(p, l, None, k, j)
        f[2] = b2
        if p * p <= bound and b2 is None:
            raise Skip('x^2 coefficient needed at %d (p^2 <= %d)' % (p, bound))
        factors[p] = f
    if missing:
        raise Skip('%d primes missing, first %s' % (len(missing), missing[:5]))
    # the stored a_n
    an = t.get('traces') or []
    M = len(an)
    while M > 1:
        try:
            calc = dirichlet_an(factors, M)
            break
        except ValueError:
            M -= 1
    n_checked = 0
    for n in range(1, M + 1):
        if an[n - 1] is None:
            continue
        if calc[n] != an[n - 1]:
            raise Skip('a_%d: factors give %d, stored %d' % (n, calc[n], an[n - 1]))
        n_checked += 1
    stats['a_n checked'] += n_checked
    complete = sum(1 for p in factors if None not in factors[p])
    return bound, factors, key, n_checked, complete


def cmd_inputs(args):
    from collections import Counter
    targets = json.load(open(args.targets))
    if args.labels:
        targets = [t for t in targets if t['label'] in args.labels]
    os.makedirs(args.outdir, exist_ok=True)
    stats, report = Counter(), []
    for t in targets:
        try:
            bound, factors, key, nck, complete = build(t, args.sweep, args.nmax, stats)
        except Skip as e:
            stats['skipped'] += 1
            report.append('%s  skipped: %s' % (t['label'], e))
            continue
        k, j = t['weight']
        name = t['label'][len('2.K.'):]
        path = os.path.join(args.outdir, 'in_%s.txt' % name)
        with open(path, 'w') as F:
            F.write(spec_line(name, t['level'], k, j, factors) + '\n')
        with open(path + '.factors.json', 'w') as F:
            json.dump({'level': t['level'], 'weight': [k, j], 'bound': bound, 'smf_label': t['label'],
                       'fricke': t['fricke'], 'factors': {str(p): f for p, f in sorted(factors.items())}}, F)
        stats['written'] += 1
        report.append('%s  ok: bound %d, %d primes (%d complete factors), sweep form %s, %d a_n checked'
                      % (t['label'], bound, len(factors), complete, key, nck))
    with open(os.path.join(args.outdir, 'report.txt'), 'w') as F:
        F.write('\n'.join(report) + '\n')
    print('\n'.join(report[-20:]))
    print(dict(stats))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('targets')
    a.add_argument('--out', required=True)
    a.add_argument('--weights', nargs='+', default=['3,0', '4,0', '3,2'])
    a.add_argument('--max-level-30', type=int, default=500, help='levels of weight (3, 0) below this')
    a.set_defaults(func=cmd_targets)
    b = sub.add_parser('inputs')
    b.add_argument('--targets', required=True)
    b.add_argument('--sweep', action='append', required=True, help='root holding <W>/<N>/p<p>.txt (repeatable)')
    b.add_argument('--nmax', required=True, help='nmax_spin executable: N k j bits -> ... nmax')
    b.add_argument('--outdir', required=True)
    b.add_argument('labels', nargs='*')
    b.set_defaults(func=cmd_inputs)
    args = ap.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
