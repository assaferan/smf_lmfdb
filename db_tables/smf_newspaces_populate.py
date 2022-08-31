# from sage.all import *
from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26
from smf_lmfdb.db_tables.commmon_create_table import SUBSPACE_TYPES, HECKE_TYPES
from lmfdb import db

from smf_lmfdb.Dimension_formulas.dimformula_smf_degree2_level_1 import smf_dims_degree_2_level_1
from smf_lmfdb.Hecke_Eigenvalues.Siegel_Eisenstein_series.Hecke_Eigenvalues_Siegel_Eisenstein import Hecke_Eigenvalues_Siegel_Eisenstein
from smf_lmfdb.Hecke_Eigenvalues.Klingen_Eisenstein_series.Hecke_Eigenvalues_Klingen_Eisenstein import Hecke_Eigenvalue_Traces_Klingen_Eisenstein

#import os
#cwd = os.getcwd()
#os.chdir('smf_lmfdb/Dimension_formulas')
#load('dimformula_smf_degree2_level_1.sage')
#os.chdir("../Hecke_Eigenvalues/Siegel_Eisenstein_series")
#load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')
#os.chdir("../Klingen_Eisenstein_series")
#load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')
#os.chdir(cwd)

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
       if (j % 2 == 1) or (k == 1):
           continue
       entry = smf_dims_degree_2_level_1(j,k,e)
       entry.update(common_entry_values(k,j,e))
       hecke_types = {1 : ['p'],
       		      2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein,
       		    'eis_Q' : Hecke_Eigenvalue_Traces_Klingen_Eisenstein}
       for sub in sub_funcs.keys():
          for deg in hecke_types.keys():
       	      for hecke_type in hecke_types[deg]:
              	  entry[sub + '_lambda_' + hecke_type] = get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e)
       non_existing = [sub for sub in SUBSPACE_TYPES.keys() if sub not in sub_funcs.keys()]
       for sub in non_existing:
           for hecke in HECKE_TYPES:
               entry[sub + '_lambda_' + hecke_type] = 'NULL'
       entries.append(entry)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
