"""
Add the column Nk2 = N(2k+j-3)^2(j+1)^2 (N the level, (k, j) the weight) to smf_newforms and
smf_newspaces, the analogue of Nk^2 for classical modular forms (the analytic conductor of the
spin L-function is about Nk2 / (256 pi^4)), used by the search pages
(search box and sort order) to bound the size of the forms.  New rows get it from
common_populate.entry_add_common_columns.

Dry run: writes PREFIX_{newforms,newspaces}.txt (id|Nk2 update files) and checks them.
--apply: adds the column where it is missing (or updates its description), loads the files with update_from_file, and reads
every row back.

    sage -python -m smf_lmfdb.db_tables.add_nk2 [--datafile PREFIX] [--apply]
"""
import argparse

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import Nk2
from smf_lmfdb.db_tables.fix_sk_traces import write_update_file

TABLES = {'newforms': db.smf_newforms, 'newspaces': db.smf_newspaces}
DESC = 'N(2k+j-3)^2(j+1)^2, for the level N and the weight (k,j): the analogue of Nk^2 for classical modular forms, proportional to the analytic conductor of the spin L-function up to a bounded factor'


def rows_of(table):
    rows = []
    for r in table.search({}, ['id', 'level', 'weight', 'degree']):
        assert r['degree'] == 2 and len(r['weight']) == 2, r
        rows.append({'id': r['id'], 'Nk2': Nk2(r['level'], r['weight'])})
    return rows


def main(apply=False, datafile='smf_nk2'):
    files = {}
    for name, table in TABLES.items():
        rows = rows_of(table)
        files[name] = (table, rows)
        vals = [r['Nk2'] for r in rows]
        print("%s: %d rows, Nk2 from %d to %d; column present: %s"
              % (table.search_table, len(rows), min(vals), max(vals), 'Nk2' in table.search_cols))
        assert max(vals) < 2**31
    if not apply:
        for name, (table, rows) in files.items():
            if 'Nk2' in table.col_type:
                write_update_file('%s_%s.txt' % (datafile, name), table, 'id', rows, ['Nk2'])
            else:
                with open('%s_%s.txt' % (datafile, name), 'w') as F:
                    F.write("id|Nk2\nbigint|integer\n\n")
                    for r in rows:
                        F.write("%d|%d\n" % (r['id'], r['Nk2']))
                print("wrote %s_%s.txt (%d rows)" % (datafile, name, len(rows)))
        print("dry run; pass --apply to write")
        return
    for name, (table, rows) in files.items():
        if 'Nk2' not in table.search_cols:
            table.add_column('Nk2', 'integer', description=DESC)
        elif table.column_description('Nk2') != DESC:
            table.column_description('Nk2', DESC)
        fname = '%s_%s.txt' % (datafile, name)
        write_update_file(fname, table, 'id', rows, ['Nk2'])
        table.update_from_file(fname, label_col='id')
        got = {r['id']: r.get('Nk2') for r in table.search({}, ['id', 'Nk2'])}
        bad = [r['id'] for r in rows if got.get(r['id']) != r['Nk2']]
        print("%s: read back %d rows, %d mismatches" % (table.search_table, len(rows), len(bad)))
        assert not bad, bad[:10]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datafile", default="smf_nk2")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    main(apply=args.apply, datafile=args.datafile)
