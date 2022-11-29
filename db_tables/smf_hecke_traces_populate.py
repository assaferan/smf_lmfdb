from smf_lmfdb.db_tables.common_populate import MAX_P
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code

from lmfdb import db

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    return e

def get_form_entries(g, F, N, k, j, orb):
    entries = []
    e = {'hecke_orbit_code' : make_orbit_code(g, F, N, k, j, 0, orb)}
    f = db.smf_newforms.lucky(e, ['trace_lambda_p', 'trace_lambda_p_square'])
    aps = f['trace_lambda_p']
    aps2 = f['trace_lambda_p_square']
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
       if (j % 2 == 1) or (k == 1):
           continue
       # right now we only have implemented forms for full level
       if (N == 1):
           query = {'hecke_orbit_code' : make_orbit_code(g, F, N, k, j, 0, 0)}
           num_orbits = db.smf_newspaces.lucky(query, ['num_forms'])['num_forms']
           for i in range(num_orbits):
               entries += get_form_entries(2, 'P', N, k, j, 0, i)
    return entries

def populate_smf_hecke_traces(triple_list):
    table = db.smf_hecke_traces
    aux_fname = "smf_lmfdb/db_tables/smf_hecke_traces_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
