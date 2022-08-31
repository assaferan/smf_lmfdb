from smf_lmfdb.db_tables.common_populate import make_space_label, table_reload, get_hecke, common_entry_values, base_26, dict_to_json, MAX_P, MAX_P_SQUARE
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.db_tables.sage_functions import smf_dims_degree_2_level_1, Hecke_Eigenvalues_Siegel_Eisenstein

from lmfdb import db

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    if (e['char_orbit_index'] == 2):
        e['family'] = 'P'
        e['level'] = 2
    e['char_orbit_label'] = base_26(e['char_orbit_index'])
    space_label = make_space_label(e)
    dummy = e.pop('char_orbit_label')
    hecke_orbit = ext_data['num_forms']
    e['label'] = space_label + '.' + base_26(hecke_orbit)
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], hecke_orbit)
    # for now our forms are always over the rational field
    # TODO - have Fabien write down the dimension of each form
    e['field_poly'] = [0,1]
    e['hecke_ring_rank'] = 1
    e['hecke_ring_power_basis'] = True
    e['hecke_ring_cyclotomic_generator'] = 0
    e['hecke_ring_character_values'] = 'NULL'
    e['hecke_ring_numerators'] = 'NULL'
    e['hecke_ring_denominators'] = 'NULL'
    e['hecke_ring_inverse_numerators'] = 'NULL'
    e['hecke_ring_inverse_denominators'] = 'NULL'
    e['an'] = 'NULL'
    e['maxp'] = MAX_P
    e['maxp_square'] = MAX_P_SQUARE
    return e

def populate_smf_hecke_nf(triple_list):
    table = db.smf_hecke_nf
    aux_fname = "smf_hecke_nf_table.dat"
    entries = []
    for triple in triple_list:
       k,j,e = triple
       if (j % 2 == 1) or (k == 1):
           continue
       entry = common_entry_values(k,j,e)
       hecke_types = {1 : ['p'],
                      2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein }
       dims = smf_dims_degree_2_level_1(j,k,e)
       for sub in sub_funcs.keys():
           if dims[sub + '_dim'] == 0:
               continue
           entry_sub = {key : entry[key] for key in entry.keys()}
           for deg in hecke_types.keys():
               for hecke_type in hecke_types[deg]:
                   key = 'lambda_' + hecke_type
                   entry_sub[key] = dict_to_json(get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e,ret_type='dict'))
           entries.append(entry_sub)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
