"""
Load the L-function rows written by lfunc_rows.py (rows.json) into lfunc_lfunctions,
lfunc_instances, lfunc_plot and lfunc_search, and read them back.

lfunc_rows.py writes no lfunc_search row. search_row() derives it from the
lfunc_lfunctions and lfunc_instances rows:
- the columns shared with lfunc_lfunctions are copied, except central_character,
  which lfunc_search gives as 1.1 for the trivial character;
- euler_factors is cut to the primes below 100 and each factor is padded with zeros
  to degree + 1 entries; euler2, ..., euler97 hold the same factors;
- dirichlet_coefficients is a_1, ..., a_100, computed from those factors;
- instance_types and instance_urls list the instances, and is_instance_X is true when
  X is one of them. There is no is_instance_SMF column, so an SMF row is found
  through instance_types.
--check rebuilds the stored search rows of four genus-2 L-functions with this rule
and compares them with the database.

The rows go in with insert_many, not copy_from: the LMFDB space check budgets a
copy_from into lfunc_lfunctions (about 300 GB) as a rewrite of the whole table, and
refuses it. lfunc_lfunctions and lfunc_search share ids, and insert_many gives the
new rows the next ids of each, so --apply checks that they agree.

Usage (from the lmfdb root):
    sage -python -m smf_lmfdb.db_tables.load_lfunctions --check
    sage -python -m smf_lmfdb.db_tables.load_lfunctions --outdir D DIR...            # dry run
    sage -python -m smf_lmfdb.db_tables.load_lfunctions --outdir D DIR... --apply
    sage -python -m smf_lmfdb.db_tables.load_lfunctions --readback DIR...
Each DIR holds a rows.json. The dry run checks the rows, refuses labels, Lhashes or
urls already in the database, and writes upload files for inspection to D.
"""
import argparse
import json
import os
import sys
from decimal import Decimal

from lmfdb import db
from smf_lmfdb.db_tables.lfunc_rows import encode, dirichlet_an

P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
TABLES = ['lfunc_lfunctions', 'lfunc_instances', 'lfunc_plot', 'lfunc_search']
CHECK_URLS = ['Genus2Curve/Q/587/a', 'Genus2Curve/Q/169/a', 'Genus2Curve/Q/249/a', 'Genus2Curve/Q/277/a']


def types_of(t):
    return {c: v for c, v in db[t].col_type.items() if c != 'id'}


def search_row(L, insts):
    """The lfunc_search row of the lfunc_lfunctions row L with instances insts."""
    types = types_of('lfunc_search')
    S = {c: L[c] for c in types if c in L}
    deg = L['degree']
    ef = [list(f) + [0] * (deg + 1 - len(f)) for f in L['euler_factors'][:len(P)]]
    S['euler_factors'] = ef
    for i, p in enumerate(P):
        S['euler%d' % p] = ef[i] if i < len(ef) else None
    assert L['central_character'].endswith('.1'), L['central_character']   # trivial character only
    S['central_character'] = '1.1'
    if L.get('dirichlet_coefficients'):
        S['dirichlet_coefficients'] = L['dirichlet_coefficients'][:100]
    else:
        S['dirichlet_coefficients'] = dirichlet_an({p: ef[i] for i, p in enumerate(P)}, 100)[1:]
    S['instance_types'] = [i['type'] for i in insts]
    S['instance_urls'] = [i['url'] for i in insts]
    for c in types:
        if c.startswith('is_instance_'):
            S[c] = c[len('is_instance_'):] in S['instance_types']
    assert set(S) == set(types), (set(types) - set(S), set(S) - set(types))
    return S


def typed(v, typ):
    """A rows.json value as the Python type of its postgres column."""
    if v is None or typ == 'jsonb':
        return v
    base = typ[:-2] if typ.endswith('[]') else typ

    def one(x):
        if isinstance(x, list):
            return [one(y) for y in x]
        if x is None:
            return None
        if base == 'numeric':
            return Decimal(str(x))
        if base in ('smallint', 'integer', 'bigint'):
            assert int(x) == x, (x, typ)
            return int(x)
        if base in ('double precision', 'real'):
            return float(x)
        if base == 'boolean':
            assert isinstance(x, bool), (x, typ)
            return x
        if base == 'text':
            return str(x)
        raise ValueError('unhandled type %s' % typ)
    return one(v)


def same(a, b):
    """Equality of a typed expected value and a value read from the database."""
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(same(a[k], b[k]) for k in a)
    if a is None or b is None:
        return a is b
    if isinstance(a, (bool, str)) or isinstance(b, (bool, str)):
        return type(a) is type(b) and a == b
    try:
        x, y = Decimal(str(a)), Decimal(str(b))
    except Exception:
        return str(a) == str(b)
    return x == y or abs(float(x) - float(y)) <= 1e-15 * max(1, abs(float(x)))


def read_rows(dirs):
    """{table: [rows]} from the rows.json files, with the derived lfunc_search rows."""
    rows = {t: [] for t in TABLES}
    for d in dirs:
        with open(os.path.join(d, 'rows.json')) as F:
            r = json.load(F)
        assert all(m['certified'] and not m['errors'] for m in r['meta']), (d, 'not certified')
        for L in r['lfunc_lfunctions']:
            assert all(c is not None for f in L['euler_factors'] for c in f), (d, 'incomplete Euler factors')
        for t in TABLES[:3]:
            rows[t] += r[t]
        for L in r['lfunc_lfunctions']:
            rows['lfunc_search'].append(search_row(L, [i for i in r['lfunc_instances'] if i['label'] == L['label']]))
    return rows


def check():
    bad = 0
    for url in CHECK_URLS:
        lab = db.lfunc_instances.lucky({'url': url}, 'label')
        S = search_row(db.lfunc_lfunctions.lucky({'label': lab}), list(db.lfunc_instances.search({'label': lab})))
        stored = db.lfunc_search.lucky({'label': lab})
        diff = []
        for c in types_of('lfunc_search'):
            a, b = S[c], stored.get(c)
            if c in ('instance_types', 'instance_urls'):
                a, b = sorted(a), sorted(b)
            if not same(a, b):
                diff.append((c, str(a)[:60], str(b)[:60]))
        print(url, lab, 'same' if not diff else diff)
        bad += bool(diff)
    return bad


def readback(dirs):
    rows, bad = read_rows(dirs), 0
    for L in rows['lfunc_lfunctions']:
        lab = L['label']
        for t in TABLES:
            types = types_of(t)
            exp = [r for r in rows[t] if r['label'] == lab]
            got = list(db[t].search({'label': lab}, sort=[]))
            if len(got) != len(exp):
                bad += 1
                print('COUNT', lab, t, len(got), len(exp))
                continue
            for e, g in zip(exp, got):
                diff = [c for c in types if not same(typed(e.get(c), types[c]), g.get(c))]
                if diff:
                    bad += 1
                    print('DIFF', lab, t, [(c, str(e.get(c))[:50], str(g.get(c))[:50]) for c in diff])
        print(lab, 'read back')
    print('mismatches', bad)
    return bad


def write_upload(path, rows, types):
    cols = sorted(types)
    with open(path, 'w') as F:
        F.write('\t'.join(cols) + '\n')
        F.write('\t'.join(types[c] for c in cols) + '\n\n')
        for r in rows:
            F.write('\t'.join(encode(r.get(c), types[c]) for c in cols) + '\n')


def load(dirs, outdir, apply):
    rows = read_rows(dirs)
    labels = [L['label'] for L in rows['lfunc_lfunctions']]
    hashes = [L['Lhash'] for L in rows['lfunc_lfunctions']]
    urls = [i['url'] for i in rows['lfunc_instances']]
    assert len(set(labels)) == len(labels) and len(set(hashes)) == len(hashes)
    clash = [('label', t, x) for x in labels for t in TABLES if db[t].exists({'label': x})]
    clash += [('Lhash', t, x) for x in hashes for t in TABLES[:3] if db[t].exists({'Lhash': x})]
    clash += [('url', 'lfunc_instances', x) for x in urls if db.lfunc_instances.exists({'url': x})]
    if clash:
        print('already in the database:', clash)
        return 1
    os.makedirs(outdir, exist_ok=True)
    for t in TABLES:
        path = os.path.join(outdir, t + '.txt')
        write_upload(path, rows[t], types_of(t))
        print('wrote %s (%d rows)' % (path, len(rows[t])))
    for S in rows['lfunc_search']:
        print(S['label'], 'a_1..a_6', S['dirichlet_coefficients'][:6], S['instance_types'], S['instance_urls'])
    if not apply:
        print('dry run; pass --apply to load')
        return 0
    before = {t: db[t].max_id() for t in TABLES}
    for t in TABLES:
        types = types_of(t)
        db[t].insert_many([{c: typed(r.get(c), types[c]) for c in types} for r in rows[t]], restat=False)
    print('max ids before', before)
    for x in labels:
        ids = {t: db[t].lucky({'label': x}, 'id') for t in TABLES}
        print(x, ids)
        assert ids['lfunc_lfunctions'] == ids['lfunc_search'], 'lfunc_search id differs from lfunc_lfunctions'
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--check', action='store_true', help='check search_row() against stored genus-2 rows')
    ap.add_argument('--readback', action='store_true', help='compare the loaded rows with rows.json')
    ap.add_argument('--apply', action='store_true', help='insert the rows')
    ap.add_argument('--outdir', help='directory for the upload files of the dry run')
    ap.add_argument('dirs', nargs='*', help='directories with a rows.json')
    a = ap.parse_args()
    if a.check:
        sys.exit(check())
    if a.readback:
        sys.exit(readback(a.dirs))
    if not a.outdir:
        ap.error('--outdir is required')
    sys.exit(load(a.dirs, a.outdir, a.apply))


if __name__ == '__main__':
    main()
