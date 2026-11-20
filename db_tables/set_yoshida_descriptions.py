"""
Say in the column descriptions of smf_hecke_nf and smf_newforms that the
Yoshida lifts of weight k = 2 store the eigenvalues (and traces) of
p T_{p^2,0} and p T_{p^2,1}, and, when j = 0, of p^2 T_{p^2,2}: those of
T_{p^2,i} are not integral there.  The newform pages say the same.

Usage (from the lmfdb root, with write access to the db):
    sage -python -m smf_lmfdb.db_tables.set_yoshida_descriptions            # print the changes
    sage -python -m smf_lmfdb.db_tables.set_yoshida_descriptions --apply
"""
import argparse

from lmfdb import db


def note(i):
    scale = "p^2" if i == 2 else "p"
    when = " and j = 0" if i == 2 else ""
    return ("; for the Yoshida lifts of weight k = 2%s these are the eigenvalues of %s T_{p^2,%s}, which are integral"
            % (when, scale, i))


def main(apply=False):
    for table, fmt in [(db.smf_hecke_nf, "lambda_p_square_%s"), (db.smf_newforms, "trace_lambda_p_square_%s")]:
        desc = table.column_description()
        for i in range(3):
            col = fmt % i
            old = desc[col]
            if "Yoshida" in old:
                print("%s.%s already says it" % (table.search_table, col))
                continue
            new = old + note(i).replace("the eigenvalues", "the traces" if table is db.smf_newforms else "the eigenvalues")
            print("%s.%s:\n    %s\n -> %s" % (table.search_table, col, old, new))
            if apply:
                table.column_description(col, new)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--apply", action="store_true", help="write the descriptions to the database")
    main(apply=parser.parse_args().apply)
