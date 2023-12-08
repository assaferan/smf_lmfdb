from sage.all import (prime_range, nth_prime, is_square)
from smf_lmfdb.db_tables.common_populate import MAX_P, table_reload_plain
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code

from lmfdb import db

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    return e

def get_space_entries(g, F, N, k, j):
    entries = []
    e = {'hecke_orbit_code' : make_orbit_code(g, F, N, k, j, 1, 1)}
    M = db.smf_newspaces.lucky(e, ['cusp_G_lambda_p'])
    aps = M['cusp_G_lambda_p']
    if (len(aps) == 0):
        max_p = 0
    else:
        max_p = min(MAX_P, nth_prime(len(aps)))
    for n,p in enumerate(prime_range(1,max_p+1)):
        e['n'] = p
        e['trace_an'] = aps[n-1]
        entries.append(e.copy())
    return entries
    
def create_entries(triple_list):
    entries = []
    for triple in triple_list:
        k,j,N = triple
        if (N == 1):
            for F in ['K','S','P']:
                entries += get_space_entries(2, F, N, k, j)
        # manually adding (3,0,N)
        if (not is_square(N)) and (k == 3) and (j == 0) and (N < 1000):
            entries += get_space_entries(2, 'K', N, k, j)
       # else:
           # For now we don't yet have traces for the higher level spaces
           # if (N == 2) and (j % 2 == 0) and (k >= 3):
           #    entries += get_space_entries(2, 'P', N, k, j)
           # entries += get_space_entries(2, 'K', N, k, j)
    return entries

def populate_smf_hecke_newspace_traces(triple_list):
    table = db.smf_hecke_newspace_traces
    aux_fname = "smf_lmfdb/db_tables/smf_hecke_newspace_traces_table.dat"
    entries = create_entries(triple_list)
    table_reload_plain(table, entries, entry_add_columns, aux_fname, "hecke_newspace_traces")
    return
