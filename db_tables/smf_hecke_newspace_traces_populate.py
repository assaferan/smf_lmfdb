from smf_lmfdb.db_tables.common_populate import MAX_P
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code

from lmfdb import db

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    return e

def get_space_entries(g, F, N, k, j):
    entries = []
    e = {'hecke_orbit_code' : make_orbit_code(g, F, N, k, j, 0, 0)}
    M = db.smf_newspaces.lucky(e, ['cusp_G_lambda_p', 'cusp_G_lambda_p_square'])
    aps = M['cusp_G_lambda_p']
    aps2 = M['cusp_G_lambda_p_square']
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
           for F in ['K', 'S', 'P']:
               entries += get_space_entries(2, F, N, k, j)
       else:
           if (N == 2) and (j % 2 == 0) and (k >= 3):
               entries += get_space_entries(2, 'P', N, k, j)
            entries += get_space_entries(2, 'K', N, k, j)
    return entries

def populate_smf_hecke_newspace_traces(triple_list):
    table = db.smf_hecke_newspace_traces
    aux_fname = "smf_hecke_newspace_traces_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
