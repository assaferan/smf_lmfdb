from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26
from smf_lmfdb.db_tables.common_create_table import SUBSPACE_TYPES, HECKE_TYPES
from smf_lmfdb.db_tables.sage_functions import smf_dims_degree_2_level_1, smf_dims_degree_2_level_2, Hecke_Eigenvalues_Traces_Siegel_Eisenstein, Hecke_Eigenvalue_Traces_Klingen_Eisenstein, Hecke_Eigenvalue_Traces_Saito_Kurokawa, Hecke_Eigenvalue_Traces_Yoshida, num_forms_Siegel_Eisenstein, num_forms_Klingen_Eisenstein, num_forms_Saito_Kurokawa, num_forms_Yoshida
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import smf_dims_paramodular

from lmfdb import db

def entry_add_columns(e, ext_data):
    e = entry_add_common_columns(e, ext_data)
    e['label'] = make_space_label(e)
    return e

def smf_level1_space(k,j,e):
    entry = smf_dims_degree_2_level_1(j,k,e)
    entry.update(common_entry_values(k,j,e))
    hecke_types = {1 : ['p'],
       		   2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
    sub_funcs = {'eis_F' : Hecke_Eigenvalues_Traces_Siegel_Eisenstein,
       		 'eis_Q' : Hecke_Eigenvalue_Traces_Klingen_Eisenstein,
                 'cusp_P' : Hecke_Eigenvalue_Traces_Saito_Kurokawa,
                 'cusp_Y' : Hecke_Eigenvalue_Traces_Yoshida}
    for sub in sub_funcs.keys():
        for deg in hecke_types.keys():
       	    for hecke_type in hecke_types[deg]:
              	entry[sub + '_lambda_' + hecke_type] = get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e)
    # in level 1 all forms are new, so we replicate the data, and fill zeros in oldforms
    dim_columns = { prefix + '_dim' for prefix in ['total', 'eis', 'cusp'] + list(SUBSPACE_TYPES.keys()) }
    for dim_col in dim_columns:
        entry['new_' + dim_col] = entry[dim_col]
        entry['old_' + dim_col] = 0
    num_form_funcs = [num_forms_Siegel_Eisenstein,
                      num_forms_Klingen_Eisenstein,
                      num_forms_Saito_Kurokawa,
                      num_forms_Yoshida]
    entry['num_forms'] = sum([func(k,j,e) for func in num_form_funcs])
    return entry

def smf_level2_space(k,j):
    entry = smf_dims_degree_2_level_2(k,j)
    entry.update(common_entry_values(k,j,0))
    return entry

# triple_list consists of triples (k,j,N) of weight and level
# at the moment we restrict to trivial character
def create_entries(triple_list):
    entries = []
    for triple in triple_list:
        k,j,N = triple
        # we dropped the condition checking since in cmf we have the zero-dimensional space in the DB
        # so we have the same behavior here
        print("creating entries for triple (%d,%d,%d)" %(k,j,N))
        if N == 1:
            # for level 1 we have complete data of Eisenstein series
            entry = smf_level1_space(k,j,0)
        elif N == 2:
            entry = smf_level2_space(k,j)
        else:
            entry = smf_dims_paramodular(k,j,N)
        entries.append(entry)
    return entries

def populate_smf_newspaces(triple_list):
    table = db.smf_newspaces
    aux_fname = "smf_lmfdb/db_tables/smf_newspaces_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return

# old code

# triple_list consists of triples (k,j,e) of weight and character
#def populate_smf_newspaces(triple_list):
#    table = db.smf_newspaces
#    aux_fname = "smf_lmfdb/db_tables/smf_newspaces_table.dat"
#    entries = []
#    for triple in triple_list:
#       k,j,e = triple
#       if (j % 2 == 1) or (k == 1):
#           continue
#       entry = smf_dims_degree_2_level_1(j,k,e)
#       entry.update(common_entry_values(k,j,e))
#       hecke_types = {1 : ['p'],
#       		      2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
#       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein,
#       		    'eis_Q' : Hecke_Eigenvalue_Traces_Klingen_Eisenstein}
#       for sub in sub_funcs.keys():
#          for deg in hecke_types.keys():
#       	      for hecke_type in hecke_types[deg]:
#              	  entry[sub + '_lambda_' + hecke_type] = get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e)
#       non_existing = [sub for sub in SUBSPACE_TYPES.keys() if sub not in sub_funcs.keys()]
#       for sub in non_existing:
#           for hecke in HECKE_TYPES:
#               entry[sub + '_lambda_' + hecke_type] = 'NULL'
#       entries.append(entry)
#    table_reload(table, entries, entry_add_columns, aux_fname)
#    return
