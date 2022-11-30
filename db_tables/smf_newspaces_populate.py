from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26
from smf_lmfdb.db_tables.common_create_table import SUBSPACE_TYPES, HECKE_TYPES
from smf_lmfdb.db_tables.sage_functions import smf_dims_degree_2_level_1, smf_dims_degree_2_level_2, Hecke_Eigenvalues_Traces_Siegel_Eisenstein, Hecke_Eigenvalue_Traces_Klingen_Eisenstein, Hecke_Eigenvalue_Traces_Saito_Kurokawa, Hecke_Eigenvalue_Traces_Yoshida, num_forms_Siegel_Eisenstein, num_forms_Klingen_Eisenstein, num_forms_Saito_Kurokawa, num_forms_Yoshida
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import smf_dims_paramodular
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code

from lmfdb import db

def entry_add_columns(e, ext_data):
    e = entry_add_common_columns(e, ext_data)
    e['label'] = make_space_label(e)
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], 1)
    return e

def smf_level1_space(k,j,e):
    entry = smf_dims_degree_2_level_1(j,k,e)
    entry.update(common_entry_values(k,j,e+1))
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
    if (k == 20) and (j == 0) and (e == 0):
        # manually fixing the (20,0) space for demonstration
        entry['num_forms'] += 1
        entry['cusp_G_lambda_p'] = [-840960,346935960,-73262366720,-5232247240500,2617414076964400,-724277370534455340,1427823701421564744,-83773835478688698980,14156088476175218899620,146957560176221097673720]
    return entry

def smf_level2_space(k,j):
    assert (j % 2 == 0) and (k >= 3)
    entry = smf_dims_degree_2_level_2(k,j)
    entry.update(common_entry_values(k,j,2))
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
        else:
            if (N == 2) and (j % 2 == 0) and (k >= 3):
                entry = smf_level2_space(k,j)
                entries.append(entry)
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
