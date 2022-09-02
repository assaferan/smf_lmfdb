from smf_lmfdb.db_tables.common_create_table import generate_common_column_types, generate_common_col_desc, generate_table, FAMILY_DICT, HECKE_TYPES, SUBSPACE_TYPES

from lmfdb import db

# Here we list the columns for the table smf_newspaces

def generate_dim_column_names():
    dim_columns = { prefix + '_dim' for prefix in ['total', 'eis', 'cusp'] + list(SUBSPACE_TYPES.keys()) }
    old_dim_columns = {'old_' + col for col in dim_columns}
    new_dim_columns = {'new_' + col for col in dim_columns}
    dim_columns.update(old_dim_columns)
    dim_columns.update(new_dim_columns)
    return dim_columns

def generate_column_types():
    col_type = generate_common_column_types()
    for subspace in SUBSPACE_TYPES.keys():
        for hecke_type in HECKE_TYPES:
            col_type[subspace + '_lambda_' + hecke_type] = 'numeric[]'
    col_type['num_forms'] = 'integer'
    col_type['traces'] = 'integer[]'
    dim_columns = generate_dim_column_names()
    for col_name in dim_columns:
        col_type[col_name] = 'integer'
    return col_type

def generate_column_desc():
    col_desc = generate_common_col_desc()
    for subspace in SUBSPACE_TYPES.keys():
        for hecke_type in HECKE_TYPES:
            col_desc[subspace + '_lambda_' + hecke_type] = 'List of traces of Hecke operators T_' + hecke_type + ' on the Siegel-Eisenstein (F) space corresponding to the values up to 200'
    col_desc['num_forms'] = 'number of Hecke orbits (each corresponds to a Galois conjugacy class of modular forms)'
    col_desc['traces'] = 'integer coefficients a_n of the trace form (sum of all newforms in the space) for n from 1 to 1000, only set when dim > 0 and not yet computed in every case.'
    col_desc['label'] = 'Label g.C.N.w.a of this newspace'
    for prefix in ['old_', 'new_', '']:
        pre = prefix[:-1]
        col_desc[prefix+'total_dim'] = 'Total dimension of the space of ' + pre +  ' modular forms'
        col_desc[prefix+'eis_dim'] = 'Dimension of the space of ' + pre + ' Eisenstein series'
        col_desc[prefix+'cusp_dim'] = 'Dimension of the space of ' + pre + ' cusp forms'
        for subspace in SUBSPACE_TYPES.keys():
            col_desc[prefix+subspace+'_dim'] = 'Dimension of the space of ' + pre + ' ' + SUBSPACE_TYPES[subspace]
    return col_desc

def create_table_smf_newspaces():
    table_name = "smf_newspaces"
    table_desc = "Spaces of Siegel modular forms" 
    generate_table(table_name, table_desc,
                   generate_column_types, generate_column_desc)
    return
