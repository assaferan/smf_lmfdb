from smf_lmfdb.db_tables.common_populate import MAX_P, table_reload
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.db_tables.sage_functions import Get_All_Hecke_Eigenvalues_Up_To

from lmfdb import db

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    return e

def get_space_entries(g, F, N, k, j):
    entries = []
    e = {'hecke_orbit_code' : make_orbit_code(g, F, N, k, j, 1, 1)}
    # !! TODO : For now, we take cusp_Y until we have cusp_G
    M = db.smf_newspaces.lucky(e, ['cusp_Y_lambda_p', 'cusp_Y_lambda_p_square'])
    if 'cusp_Y_lambda_p' not in M:
        print((k, j, N))
    aps = M['cusp_Y_lambda_p']
    aps2 = M['cusp_Y_lambda_p_square']
    an = Get_All_Hecke_Eigenvalues_Up_To(MAX_P+1, aps, aps2, (k,j))
    for n in range(1,MAX_P+1):
        e['n'] = n
        e['trace_an'] = an[n-1]
        entries.append(e)
    return entries
    
def create_entries(triple_list):
    entries = []
    for triple in triple_list:
       k,j,N = triple
       if (N == 1):
           entries += get_space_entries(2, 'P', N, k, j)
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
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
