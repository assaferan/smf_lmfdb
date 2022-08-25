import os
cwd = os.getcwd()
os.chdir("../../lmfdb")
from lmfdb import db
os.chdir(cwd)

FAMILY_DICT = {
    'paramodular' : 'K',
    'Siegel'      : 'S',
    'principal'   : 'P'
}

table_name = "smf_newspaces"
table_desc = "Spaces of Siegel modular forms" 

# Here we list the columns for the table smf_newspaces

def generate_dim_column_names():
    dim_columns = {'total_dim', 'eis_dim', 'cusp_dim', 'eis_Q_dim', 'eis_F_dim', 'cusp_P_dim', 'cusp_Y_dim', 'cusp_G_dim'}
#    old_dim_columns = {'old_' + col for col in dim_columns}
#    new_dim_columns = {'new_' + col for col in new_dim_columns}
#    dim_columns.update(old_dim_columns)
#    dim_columns.update(new_dim_columns)
    return dim_columns

def generate_column_types(dim_columns):
    col_type = {}
    col_type['degree'] = 'smallint'
    col_type['family'] = 'text'
    col_type['level'] = 'integer'
    col_type['weight'] = 'smallint[]'
    col_type['char_orbit'] = 'smallint'
    col_type['char_orbit_label'] = 'text'
    col_type['char_order'] = 'integer'
    col_type['char_degree'] = 'integer'
    col_type['conrey_indexes'] = 'integer[]'
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    for subspace in ['eis_F']:
        for hecke_type in hecke_types:
            col_type[subspace + '_lambda_' + hecke_type] = 'numeric[]'
#    col_type['num_forms'] = 'integer'
#    col_type['traces'] = 'integer[]'
    col_type['label'] = 'text'
    for col_name in dim_columns:
        col_type[col_name] = 'integer'
    return col_type

def generate_search_columns(col_type):
    search_col = {}
    types = set(col_type.values())
    for t in types:
        search_col[t] = [k for k in col_type.keys() if col_type[k] == t]
    return search_col

def generate_col_desc():
    col_desc = {}
    col_desc['degree'] = 'Degree g of this newform (automorphic with repsect to Sp(2g, Q))'
    col_desc['family'] = "Family of arithmetic subgroups ('K' = paramodular, 'S' = Siegel, 'P' =  principal)"
    col_desc['level'] = 'Level N in the family'
    col_desc['weight'] = 'Weight of this newform (highest weight of the corresponding irreducible representation of GL(g))'
    col_desc['char_orbit'] = 'ordinal i identifying the Galois orbit of the character of this newform (base26 encoded in the newform label / character orbit label)'
    col_desc['char_orbit_label'] = 'alphabetic representation of the orbit'
    col_desc['char_order'] = 'the order of the character chi'
    col_desc['char_degree'] = 'Degree of the (cyclotomic) character field'
    col_desc['conrey_indexes'] = 'Sorted list of Conrey indexes of characters in this Galois orbit'
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    for subspace in ['eis_F']:
        for hecke_type in hecke_types:
            col_desc[subspace + '_lambda_' + hecke_type] = 'List of traces of Hecke operators T_' + hecke_type + ' on the Siegel-Eisenstein (F) space corresponding to the values up to 200'
#    col_desc['num_forms'] = 'number of Hecke orbits (each corresponds to a Galois conjugacy class of modular forms)'
#    col_desc['traces'] = 'integer coefficients a_n of the trace form (sum of all newforms in the space) for n from 1 to 1000, only set when dim > 0 and not yet computed in every case.'
    col_desc['label'] = 'Label g.C.N.w.a of this newspace'
#    for prefix in ['old_', 'new_', '']:
    for prefix in ['']:
        pre = prefix[:-1]
        col_desc[prefix+'total_dim'] = 'Total dimension of the space of ' + pre +  ' modular forms'
        col_desc[prefix+'eis_dim'] = 'Dimension of the space of ' + pre + ' Eisenstein series'
        col_desc[prefix+'cusp_dim'] = 'Dimension of the space of ' + pre + ' cusp forms'
        col_desc[prefix+'eis_F_dim'] = 'Dimension of the space of ' + pre + ' Siegel-Eisenstein series'
        col_desc[prefix+'eis_Q_dim'] = 'Dimension of the space of ' + pre + ' Klingen-Eisenstein series'
        col_desc[prefix+'cusp_P_dim'] = 'Dimension of the space of ' + pre + ' Saito-Kurokawa lifts'
        col_desc[prefix+'cusp_Y_dim'] = 'Dimension of the space of ' + pre + ' Yoshida lifts'
        col_desc[prefix+'cusp_G_dim'] = 'Dimension of the space of ' + pre + ' cuspforms of general type'
    return col_desc

def generate_table():
    if table_name in db.tablenames:
        db.drop_table(table_name)
    dim_columns = generate_dim_column_names()
    col_type = generate_column_types(dim_columns)
    search_col = generate_search_columns(col_type)
    col_desc = generate_col_desc()
    db.create_table(table_name, search_col, 'label', table_desc, col_desc)
    return
