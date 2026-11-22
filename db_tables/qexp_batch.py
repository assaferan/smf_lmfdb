"""
Fourier coefficients of whole spaces S_k(K(N)) from Brandon Williams' eigenform files, with
qexp_from_eigenforms, keeping only the spaces all of whose db newforms get coefficients.

Stage 1, one file per call (run them in parallel; each process opens its own db connection):

    sage -python -m smf_lmfdb.db_tables.qexp_batch file OUTDIR FILE

writes OUTDIR/FILE.json: the matched orbit and its normalized coefficients, or the error (an
oldform at a composite level matches no orbit of the same dimension, for instance).

Stage 2:

    sage -python -m smf_lmfdb.db_tables.qexp_batch spaces OUTDIR PREFIX

For each space (N, k) with files: complete if each of its db newforms is matched by exactly one
file, with unit content and coefficients consistent on every class.  Files that match nothing
are reported (expected: oldforms).  smf_qexp_reduction is per level, shared by the weights, so
the classes of a level are made once, from the box of its files (all of the same precision), and
compared with the rows the db already has for that level: if it has them they must be the same.
Orbits that already have coefficients in the db are skipped.  Writes PREFIX.json (reduction and
coefficient rows to insert, in the format of qexp_from_eigenforms) and PREFIX_spaces.json (the
status of each space).  Nothing is written to the db.
"""
import glob
import json
import os
import re
import sys
import traceback
from collections import Counter, defaultdict

from lmfdb import db

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qexp_from_eigenforms import process_file, classes, normalize, space_eigenforms, INT_MAX  # noqa: E402


def one_file(outdir, fname):
    out = os.path.join(outdir, os.path.basename(fname) + '.json')
    try:
        N, k, f, co, prec, info = process_file(fname, Counter())
        res = {'ok': True, 'level': N, 'weight': k, 'prec': prec, 'label': f['label'],
               'hecke_orbit_code': f['hecke_orbit_code'], 'info': info,
               'co': [[list(m), [int(x) for x in v]] for m, v in co.items()]}
    except Exception as e:
        N, k = map(int, re.search(r'level_(\d+)_weight_(\d+)_label_\d+(?:_wrong_label)?\.txt$', fname).groups())
        res = {'ok': False, 'level': N, 'weight': k, 'file': os.path.basename(fname),
               'error': '%s: %s' % (type(e).__name__, e), 'trace': traceback.format_exc()[-1500:]}
    json.dump(res, open(out, 'w'))


def one_level(outdir, N, k, edir):
    """All the files of S_k(K(N)) together (qexp_from_eigenforms.space_eigenforms): writes
    OUTDIR/level_N_weight_k_space.json (the diagnostics) and one result per eigenform found, as
    one_file does."""
    # the files named *_wrong_label.txt are used too: only the span of the files matters here
    fnames = sorted(glob.glob(os.path.join(edir, 'level_%d_weight_%d_label_*.txt' % (N, k))))
    fnames = [f for f in fnames if re.search(r'_label_\d+(_wrong_label)?\.txt$', f)]
    base = os.path.join(outdir, 'level_%d_weight_%d_space' % (N, k))
    try:
        diag, forms = space_eigenforms(fnames)
    except Exception as e:
        json.dump({'level': N, 'weight': k, 'error': '%s: %s' % (type(e).__name__, e), 'trace': traceback.format_exc()[-1500:]}, open(base + '.json', 'w'))
        return
    diag['eigenforms'] = []
    for j, (K, d) in enumerate(forms):
        name = 'level_%d_weight_%d_space_%d' % (N, k, j)
        try:
            N1, k1, f, co, prec, info = normalize(N, k, K, d, diag['prec'], name)
            res = {'ok': True, 'level': N, 'weight': k, 'prec': prec, 'label': f['label'],
                   'hecke_orbit_code': f['hecke_orbit_code'], 'info': info,
                   'co': [[list(m), [int(x) for x in v]] for m, v in co.items()]}
        except Exception as e:
            res = {'ok': False, 'level': N, 'weight': k, 'file': name, 'degree': int(K.degree()) if hasattr(K, 'degree') else 1,
                   'error': '%s: %s' % (type(e).__name__, e), 'trace': traceback.format_exc()[-1500:]}
        diag['eigenforms'].append([name, res.get('label'), res['ok']])
        json.dump(res, open(os.path.join(outdir, name + '.json'), 'w'))
    json.dump(diag, open(base + '.json', 'w'))


def spaces(outdir, prefix):
    res = [json.load(open(fn)) for fn in sorted(glob.glob(os.path.join(outdir, '*.json'))) if not fn.endswith('_space.json')]
    by_space = defaultdict(list)
    for r in res:
        by_space[(r['level'], r['weight'])].append(r)
    status, complete = [], {}
    for (N, k), rs in sorted(by_space.items()):
        forms = list(db.smf_newforms.search({'space_label': '2.K.%d.%d.0.a' % (N, k)}, ['label', 'hecke_orbit_code', 'dim']))
        ok = [r for r in rs if r['ok']]
        matched = Counter(r['hecke_orbit_code'] for r in ok)
        st = {'level': N, 'weight': k, 'db_newforms': len(forms), 'files': len(rs), 'matched': len(ok),
              'unmatched_files': [[r['file'], r['error'][:200]] for r in rs if not r['ok']],
              'not_unit_content': [[r['label'], r['info'].get('multiplier_index'), r['info'].get('c')] for r in ok if r['info']['content_index'] != 1],
              'not_normalized': [r['label'] for r in ok if not r['info'].get('generator_found', r['info']['content_index'] == 1)],
              'precs': sorted({r['prec'] for r in ok})}
        missing = [f['label'] for f in forms if f['hecke_orbit_code'] not in matched]
        dup = [c for c, n in matched.items() if n > 1]
        st['missing'] = missing
        if not forms:
            st['status'] = 'no db newforms'
        elif missing or dup or st['not_normalized'] or len(st['precs']) != 1:
            st['status'] = 'incomplete' if missing else 'problem'
        else:
            st['status'] = 'complete'
            complete[(N, k)] = ok
        status.append(st)
    # per level: the classes, shared by the weights
    by_level = defaultdict(dict)
    for (N, k), ok in complete.items():
        by_level[N][k] = ok
    red_lines, coeff_lines, summary = [], [], []
    for N, wts in sorted(by_level.items()):
        precs = {r['prec'] for ok in wts.values() for r in ok}
        if len(precs) != 1:
            summary.append({'level': N, 'error': 'precisions %s' % sorted(precs)})
            continue
        prec = precs.pop()
        reps_ref, lsum = None, {'level': N, 'weights': sorted(wts), 'prec': prec}
        values_all, bad_all = {}, Counter()
        for k, ok in sorted(wts.items()):
            orbits = {r['hecke_orbit_code']: {tuple(m): v for m, v in r['co']} for r in ok}
            cls, reps, values, bad = classes(N, prec, k, orbits)
            if reps_ref is None:
                reps_ref = reps
            elif reps != reps_ref:
                lsum['error'] = 'the classes of the weights differ'
            values_all.update(values)
            bad_all.update(bad)
        if bad_all:
            lsum['error'] = 'sign-inconsistent classes %s' % dict(bad_all)
        if 'error' in lsum:
            summary.append(lsum)
            continue
        disc = lambda key: 4 * reps_ref[key][0] * reps_ref[key][2] - reps_ref[key][1]**2
        keys = sorted(reps_ref, key=lambda key: (disc(key), key))
        rows = [{'level': N, 'family': 'K', 'index': i, 'is_minimal': True, 'qf_legendre': list(key[0]),
                 'qf_tag': list(key[1]), 'qf_orbit_rep': [reps_ref[key][2], reps_ref[key][1], reps_ref[key][0]]}
                for i, key in enumerate(keys)]
        have = {(tuple(r['qf_legendre']), tuple(r['qf_tag'])): (r['index'], r['qf_orbit_rep'])
                for r in db.smf_qexp_reduction.search({'level': N, 'family': 'K'}, ['qf_legendre', 'qf_tag', 'index', 'qf_orbit_rep'])}
        if have:
            mine = {(tuple(r['qf_legendre']), tuple(r['qf_tag'])): (r['index'], r['qf_orbit_rep']) for r in rows}
            if mine != have:
                lsum['error'] = 'the db has other classes at this level (%d rows, %d here)' % (len(have), len(mine))
                summary.append(lsum)
                continue
            lsum['reduction'] = 'already in the db'
        else:
            red_lines += rows
            lsum['reduction'] = len(rows)
        n_co = over = 0
        for code, vals in values_all.items():
            if db.smf_qexp_coeffs.exists({'hecke_orbit_code': code}):
                lsum.setdefault('already_in_db', []).append(code)
                continue
            for key in keys:
                if key in vals:
                    c = [int(x) for x in vals[key]]
                    over += any(abs(x) > INT_MAX for x in c)
                    coeff_lines.append({'hecke_orbit_code': code, 'qf_legendre': list(key[0]), 'qf_tag': list(key[1]), 'coeff': c})
                    n_co += 1
        lsum['coeff_rows'], lsum['overflow'] = n_co, over
        if over:
            lsum['error'] = 'bigint overflow'
        summary.append(lsum)
        print(json.dumps(lsum), flush=True)
    bad_levels = {s['level'] for s in summary if 'error' in s}
    hecke_nf = []
    for (N, k), ok in complete.items():
        if N in bad_levels:
            continue
        for r in ok:
            i = r['info']
            h = db.smf_hecke_nf.lucky({'hecke_orbit_code': r['hecke_orbit_code']}, ['label', 'hecke_ring_numerators', 'hecke_ring_denominators', 'field_poly'])
            if i.get('order'):
                nums, dens = i['order']
            elif len(h['field_poly']) == 2:
                nums, dens = [[1]], [1]
            else:
                nums, dens = h['hecke_ring_numerators'], h['hecke_ring_denominators']
            hecke_nf.append({'label': h['label'], 'qexp_prec': r['prec'], 'qexp_order_index': int(i.get('multiplier_index', 1)),
                             'qexp_scale': int(i.get('c', 1)), 'qexp_order_numerators': nums, 'qexp_order_denominators': dens})
    red_lines = [r for r in red_lines if r['level'] not in bad_levels]
    lvl = {c: N for N, wts in by_level.items() for ok in wts.values() for c in (r['hecke_orbit_code'] for r in ok)}
    coeff_lines = [r for r in coeff_lines if lvl[r['hecke_orbit_code']] not in bad_levels]
    json.dump({'spaces': status, 'levels': summary}, open(prefix + '_spaces.json', 'w'), indent=1)
    json.dump({'levels': summary, 'reduction': red_lines, 'coeffs': coeff_lines, 'hecke_nf': hecke_nf}, open(prefix + '.json', 'w'))
    print(json.dumps(Counter(s['status'] for s in status)))
    print("levels with errors: %s" % sorted(bad_levels))
    print("wrote %s.json: %d reduction rows, %d coefficient rows, %d hecke_nf rows" % (prefix, len(red_lines), len(coeff_lines), len(hecke_nf)))


if __name__ == '__main__':
    if sys.argv[1] == 'file':
        one_file(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'level':
        one_level(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif sys.argv[1] == 'spaces':
        spaces(sys.argv[2], sys.argv[3])
