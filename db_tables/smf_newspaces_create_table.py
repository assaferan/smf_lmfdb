from smf_lmfdb.db_tables.common_create_table import generate_common_column_types, generate_common_col_desc, hecke_types, generate_table

from lmfdb import db

FAMILY_DICT = {
    'paramodular' : 'K',
    'Siegel'      : 'S',
    'principal'   : 'P'
}

# Here we list the columns for the table smf_newspaces

def generate_dim_column_names():
    dim_columns = {'total_dim', 'eis_dim', 'cusp_dim', 'eis_Q_dim', 'eis_F_dim', 'cusp_P_dim', 'cusp_Y_dim', 'cusp_G_dim'}
#    old_dim_columns = {'old_' + col for col in dim_columns}
#    new_dim_columns = {'new_' + col for col in new_dim_columns}
#    dim_columns.update(old_dim_columns)
#    dim_columns.update(new_dim_columns)
    return dim_columns

def generate_column_types():
    col_type = generate_common_column_types()
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    for subspace in ['eis_F', 'eis_Q']:
        for hecke_type in hecke_types:
            col_type[subspace + '_lambda_' + hecke_type] = 'numeric[]'
#    col_type['num_forms'] = 'integer'
#    col_type['traces'] = 'integer[]'
    dim_columns = generate_dim_column_names()
    for col_name in dim_columns:
        col_type[col_name] = 'integer'
    return col_type

def generate_column_desc():
    col_desc = generate_common_col_desc()
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    for subspace in ['eis_F', 'eis_Q']:
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

def create_table_smf_newspaces():
    table_name = "smf_newspaces"
    table_desc = "Spaces of Siegel modular forms" 
    generate_table(table_name, table_desc,
                   generate_column_types, generate_column_desc)
    return
