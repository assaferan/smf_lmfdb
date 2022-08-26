from sage.all import *
from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26
from lmfdb import db
import os
cwd = os.getcwd()
os.chdir('smf_lmfdb/Dimension_formulas')
load('dimformula_smf_degree2_level_1.sage')
os.chdir("../Hecke_Eigenvalues/Siegel_Eisenstein_series")
load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')
os.chdir("../Klingen_Eisenstein_series")
load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')
os.chdir(cwd)

def entry_add_columns(e, ext_data):
    e = entry_add_common_columns(e, ext_data)
    e['label'] = make_space_label(e)
    return e

# triple_list consists of triples (k,j,e) of weight and character
def populate_smf_newspaces(triple_list):
    table = db.smf_newspaces
    aux_fname = "smf_lmfdb/db_tables/smf_newspaces_table.dat"
    entries = []
    for triple in triple_list:
       k,j,e = triple
       entry = dim_splitting_smf_degree_2_level_1(j,k,e)
       entry.update(common_entry_values(k,j,e))
       hecke_types = {1 : ['p'],
       		      2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein_Series_All,
       		    'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac}
       for sub in sub_funcs.keys():
          for deg in hecke_types.keys():
       	      for hecke_type in hecke_types[deg]:
              	  entry[sub + '_lambda_' + hecke_type] = get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e)
       entries.append(entry)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
