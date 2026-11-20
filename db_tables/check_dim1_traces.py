"""
Check the population path of the paramodular forms of type G and dimension 1:
the traces made by smf_newforms_populate.entry_add_columns from the output of
Hecke_Eigenforms_paramodular must satisfy a_{p^2} = T(p^2) + p^(mu-1),
mu = j + 2k - 3, at the good primes.

trace_lambda_p_square holds T(p^2) once Hecke_Eigenforms_paramodular has
converted omf5_data's T(p,2); entry_add_columns used to convert it a second
time. For the form of weight (4,0) and level 31, T(2) = -15 and T(4) = 33, so
a_4 = 33 + 2^4 = 49 (the double conversion gave -79).

Run from the root of the lmfdb checkout:
    sage -python -m smf_lmfdb.db_tables.check_dim1_traces [k j N ...]
"""
import sys

from sage.all import prime_range

from smf_lmfdb.db_tables.common_populate import common_entry_values
from smf_lmfdb.db_tables.smf_newforms_populate import entry_add_columns
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import Hecke_Eigenforms_paramodular


def check(k, j, N):
    mu = j + 2*k - 3
    checked = wrong = 0
    for i, form in enumerate(Hecke_Eigenforms_paramodular(k, j, N)):
        if form.get('aut_rep_type') != 'G' or form['dim'] != 1:
            continue
        e = common_entry_values(k, j, N, 'K')
        e.update(form)
        e = entry_add_columns(e, {'id': 0, 'num_forms': i + 1})
        for idx, p in enumerate(prime_range(len(e['traces']) + 1)):
            if p*p > len(e['traces']):
                break
            t2 = e['trace_lambda_p_square'][idx] if idx < len(e['trace_lambda_p_square']) else 'NULL'
            if N % p == 0 or t2 == 'NULL':
                continue
            checked += 1
            if e['traces'][p*p - 1] != t2 + p**(mu - 1):
                wrong += 1
                print("wrong: level %s weight (%s,%s) orbit %s: a_%s = %s, expected %s"
                      % (N, k, j, i + 1, p*p, e['traces'][p*p - 1], t2 + p**(mu - 1)))
    print("weight (%s,%s) level %s: %s coefficients a_{p^2} checked, %s wrong" % (k, j, N, checked, wrong))
    return wrong == 0


if __name__ == '__main__':
    args = [int(a) for a in sys.argv[1:]] or [4, 0, 31]
    triples = [args[i:i+3] for i in range(0, len(args), 3)]
    sys.exit(0 if all(check(*t) for t in triples) else 1)
