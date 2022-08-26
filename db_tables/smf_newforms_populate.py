from sage.all import *
from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26

import os
cwd = os.getcwd()
os.chdir("smf_lmfdb/Hecke_Eigenvalues/Siegel_Eisenstein_series")
load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')
os.chdir("../Klingen_Eisenstein_series")
load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')
os.chdir(cwd)

from lmfdb import db

def make_orbit_code(g, F, N, k, j, i, X):
    return g + (ord(F)<<8) + (N<<12) + (k<<20) + (j<<28) + ((i-1)<<36) + ((X-1)<<52)

def entry_add_columns(e, ext_data):
    e = entry_add_common_columns(e, ext_data)
    e['space_label'] = make_space_label(e)
    e['hecke_orbit'] = ext_data['num_forms']
    e['label'] = e['space_label'] + '.' + base_26(e['hecke_orbit'])
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], e['hecke_orbit'])
    # for now our forms are always over the rational field
    # TODO - have Fabien write down the dimension of each form
    e['dim'] = 1
    e['relative_dim'] = e['dim'] // e['char_degree']
    e['nf_label'] = '1.1.1.1'
    e['field_disc'] = 1
    e['field_poly_is_cyclotomic'] = False
    e['field_poly_is_real_cyclotomic'] = False
    e['field_poly'] = [0,1]
    e['field_poly_root_of_unity'] = 0
    e['field_disc_factorization'] = []
    e['is_cuspidal'] = (e['aut_rep_type'] in ['Y', 'P', 'G'])
    e['lift_type'] = e['aut_rep_type']
    # for now we don't populate these fields
    e['trace_hash'] = 'NULL'
    e['analytic_rank'] = 'NULL'
    e['analytic_rank_proved'] = False
    e['qexp_display'] = 'NULL'
    e['related_objects'] = 'NULL'
    e['embedded_related_objects'] = 'NULL'
    e['trace_display'] = 'NULL'
    e['traces'] = 'NULL'
    return e

def populate_smf_newforms(triple_list):
    table = db.smf_newforms
    aux_fname = "smf_newforms_table.dat"
    entries = []
    for triple in triple_list:
       k,j,e = triple
       entry = common_entry_values(k,j,e)
       hecke_types = {1 : ['p'],
                      2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein_Series_All,
                    'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac}
       for sub in sub_funcs.keys():
           entry_sub = {key : entry[key] for key in entry.keys()}
           entry_sub['aut_rep_type'] = sub[-1]
           for deg in hecke_types.keys():
               for hecke_type in hecke_types[deg]:
                   key = 'trace_lambda_' + hecke_type
                   entry_sub[key] = get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e)
           entries.append(entry_sub)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
