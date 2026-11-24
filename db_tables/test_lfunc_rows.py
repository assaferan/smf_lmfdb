"""
Checks of lfunc_rows.py, runnable with plain python3 (or pytest) from the repo root:

    python3 db_tables/test_lfunc_rows.py

The genus-2 checks run lfunc_rows on smf_driver output for the genus-2 curves 169.a
(non-primitive, = CMF 13.2.e.a) and 587.a (primitive, rank 1), viewed as paramodular
forms of weight (2, 0), and compare with the LMFDB rows saved in
lfunc_example/lmfdb_reference_genus2.json (read from grace on 2026-11-24).
"""
import bisect
import json
import os
import sys
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lfunc_rows as LR  # noqa: E402

EX = os.path.join(HERE, 'lfunc_example')
getcontext().prec = 60


def load(name):
    with open(os.path.join(EX, name)) as F:
        return json.load(F)


def test_formulas():
    # 2.K.61.3.0.a.a: T(2) = -7, T(2,2) = 7; the bad factor at 61 from c = 144, eps = -1
    assert LR.good_factor(2, -7, 7, 3, 0) == [1, 7, 24, 56, 64]
    assert LR.bad_factor(1, 61, eps=-1, c=144) == [1, -144, 232044, -13845841]
    assert LR.gamma_data(3, 0) == (3, [-1, 0, 0, 1], [[], [-1, 0]], [1, 3])
    assert LR.gamma_data(2, 0) == (1, [0, 0, 1, 1], [[], [0, 0]], [1, 1])
    assert LR.spectral([1, 3]) == 'c1c3-0' and LR.spectral([1, 1]) == 'c1e2-0'
    assert [LR.conductor_str(n) for n in (61, 169, 8192, 1152, 1225)] == \
        ['61', '13e2', '2e13', '1152', '35e2']
    assert LR.smf_url('2.K.61.3.0.a.a') == 'ModularForm/GSp/Q/2/K/61/3.0/a/a'
    # the 61 input line agrees with the formulas where T(p,2) is known (p^2 <= 230)
    fac = LR.factors_from_spec(os.path.join(EX, 'in61_bound230.txt'), 232, [61])
    assert fac[61] == LR.bad_factor(1, 61, eps=-1, c=144)
    for p in (2, 3, 5, 7, 11, 13):
        f = fac[p]
        l = -f[1]
        b = f[2] // p - p * p - 1
        assert LR.good_factor(p, l, b, 3, 0) == f
    # a_p = T(p) eigenvalues of 2.K.61.3.0.a.a (smf_newforms.trace_lambda_p)
    an = LR.dirichlet_an(fac, 10)
    assert [an[p] for p in (2, 3, 5, 7)] == [-7, -3, 3, -9]
    # analytic conductor, against the LMFDB values
    ref = load('lmfdb_reference_genus2.json')
    for key, N in (('CMF_169', 169), ('G2Q_587', 587)):
        ac = LR.analytic_conductor(N, [], [1, 1])
        assert abs(ac / ref[key]['analytic_conductor'] - 1) < 1e-12, (key, ac)


def z_plot_points(rec, dual=False, rmax=30, step=0.02):
    """lmfdb.lfunctions.Lfunctionutilities.z_plot_points from LMFDB/lmfdb#7272."""
    control = ([(x, 0.0, d) for x, d in rec["zeros"]]
               + [(x, y, 0.0) for x, y in rec["extrema"]]
               + [tuple(c) for c in rec["extra_control"]])
    control.sort()
    xs = [c[0] for c in control]
    r = min(-xs[0], xs[-1], rmax)
    n = max(1, int(2 * r / step))
    grid = sorted(set([-r + 2 * r * i / n for i in range(n + 1)] + [x for x in xs if -r < x < r]))
    points = []
    for x in grid:
        i = min(max(bisect.bisect_right(xs, x), 1), len(xs) - 1)
        (x0, y0, d0), (x1, y1, d1) = control[i - 1], control[i]
        h = x1 - x0
        t = (x - x0) / h if h else 1.0
        points.append((x, (2*t**3 - 3*t**2 + 1) * y0 + (t**3 - 2*t**2 + t) * h * d0
                       + (3*t**2 - 2*t**3) * y1 + (t**3 - t**2) * h * d1))
    return points


def compare_genus2(name, refkeys, plotkey, verbose=True):
    drv = load('driver_g2_%s.json' % name)
    N = drv['conductor']
    bad = [p for p, _ in LR.prime_factors(N)]
    fac = LR.factors_from_spec(os.path.join(EX, 'in_g2_%s.txt' % name), 10 ** 6, bad)
    row, inst, plot, meta = LR.lfunction_rows(drv, fac, '2.K.%d.2.0.a.a' % N, 2, 0,
                                              allow_uncertified=True)
    ref = load('lmfdb_reference_genus2.json')
    out = {'name': name, 'driver_errors': drv['errors'], 'nmax': drv['nmax_required']}
    ours = [Decimal(z['mid']) for z in drv['zeros']]
    for key in refkeys:
        R = ref[key]
        if R.get('label'):
            assert row['label'] == R['label'], (row['label'], R['label'])
        assert row['spectral_label'] == R.get('spectral_label', row['spectral_label'])
        assert row['root_number'] == R['root_number']
        assert row['order_of_vanishing'] == R['order_of_vanishing']
        if R.get('leading_term'):
            d = abs(Decimal(row['leading_term']) - Decimal(R['leading_term']))
            out[key + ' leading_term diff'] = float(d)
            assert d < Decimal('1e-11')
        # every LMFDB zero is one of ours
        worst = Decimal(0)
        for z in R['positive_zeros']:
            zz = Decimal(str(z))
            worst = max(worst, min(abs(zz - o) for o in ours))
        out[key + ' zeros'] = len(R['positive_zeros'])
        out[key + ' worst zero diff'] = float(worst)
        tol = Decimal('1e-28') if isinstance(R['positive_zeros'][0], str) else Decimal('1e-13')
        assert worst < tol, (key, worst)
        # and we find no zero they do not have, below their last one
        top = max(Decimal(str(z)) for z in R['positive_zeros'])
        assert sum(1 for o in ours if o <= top + Decimal('1e-9')) == len(R['positive_zeros']), key
        assert [f for f in row['euler_factors']][:len(R['euler_factors'])] == R['euler_factors'][:25]
        assert row['bad_lfactors'] == R['bad_lfactors']
        if isinstance(R.get('a2'), int):
            assert [row['a%d' % n] for n in range(2, 11)] == [R['a%d' % n] for n in range(2, 11)]
        # plot_values, where both exist (glfunc's plot stops at t = 16 in degree 4)
        S = LR.Samples(drv['plot']['spacing'], drv['plot']['values'])
        delta = float(R['plot_delta'])
        diffs = [abs(S.value(i * delta) - v) for i, v in enumerate(R['plot_values'])
                 if i * delta <= S.T]
        out[key + ' plot_values max diff (t<=16)'] = max(diffs)
        # CMF rows of dimension > 1 multiply spline interpolations of their factors' plots
        assert max(diffs) < (1e-6 if key.startswith('G2Q') else 5e-3), (key, max(diffs))
    if plotkey and ref.get(plotkey):
        ours_pts = dict(z_plot_points(plot, rmax=15))
        theirs = z_plot_points(ref[plotkey], rmax=15)
        S = LR.Samples(drv['plot']['spacing'], drv['plot']['values'])
        scale = max(abs(v) for v in S.y)
        mx = max(abs(y - (S.value(abs(x)) * (1 if x >= 0 else int(row['root_number']))))
                 for x, y in theirs)
        mo = max(abs(y - (S.value(abs(x)) * (1 if x >= 0 else int(row['root_number']))))
                 for x, y in ours_pts.items())
        out['lfunc_plot: LMFDB curve vs glfunc samples, max diff on [-15,15]'] = mx
        out['lfunc_plot: our curve vs glfunc samples, max diff on [-15,15]'] = mo
        out['max |Z| on [0,16]'] = scale
        assert mo < 2e-3 * scale
    out['our label'] = row['label']
    out['our Lhash (round(z1*2^100))'] = row['Lhash']
    out['analytic_conductor'] = row['analytic_conductor']
    if verbose:
        for k, v in out.items():
            print('  %s: %s' % (k, v))
    return out


def test_genus2_169():
    compare_genus2('169a', ['G2Q_169', 'CMF_169'], 'plot_169', verbose=False)


def test_genus2_587():
    compare_genus2('587a', ['G2Q_587'], 'plot_587', verbose=False)


def test_61_uncertified_rows():
    drv = load('driver_61a_uncertified_bound230.json')
    fac = LR.factors_from_spec(os.path.join(EX, 'in61_bound230.txt'), drv['nmax_used'], [61])
    try:
        LR.lfunction_rows(drv, fac, '2.K.61.3.0.a.a', 3, 0)
        assert False, "uncertified run accepted"
    except ValueError:
        pass
    row, inst, plot, meta = LR.lfunction_rows(drv, fac, '2.K.61.3.0.a.a', 3, 0, accuracy=40,
                                              allow_uncertified=True)
    assert row['label'] == '4-61-1.1-c1c3-0-0'
    assert row['root_number'] == str((-1) ** 3 * -1)          # (-1)^k * Fricke
    assert inst['url'] == 'ModularForm/GSp/Q/2/K/61/3.0/a/a' and inst['type'] == 'SMF'
    assert row['euler_factors'][17] == [1, -144, 232044, -13845841]   # p = 61
    assert row['euler_factors'][6][2] is None                         # p = 17, T(17,2) unknown


if __name__ == '__main__':
    test_formulas()
    print('formulas: ok')
    print('169.a:')
    compare_genus2('169a', ['G2Q_169', 'CMF_169'], 'plot_169')
    print('587.a:')
    compare_genus2('587a', ['G2Q_587'], 'plot_587')
    test_61_uncertified_rows()
    print('61.a uncertified rows: ok')
