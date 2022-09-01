from smf_lmfdb.db_tables.common_populate import make_space_label, table_reload, get_hecke, common_entry_values, base_26, MAX_P, MAX_P_SQUARE
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.db_tables.sage_functions import Hecke_Eigenvalues_Siegel_Eisenstein, Hecke_Eigenvalues_Klingen_Eisenstein, Get_All_Hecke_Eigenvalues_Up_To

from lmfdb import db

# TODO - This is now assuming that everything is in Q!!
def nf_lists_to_elements(coeffs):
    return [coeff[0] for coeff in coeffs]

def nf_elts_to_lists(coeffs):
    return [[coeff] for coeff in coeffs]

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    e['char_orbit_label'] = base_26(e['char_orbit_index'])
    space_label = make_space_label(e)
    dummy = e.pop('char_orbit_label')
    hecke_orbit = ext_data['num_forms']
    e['label'] = space_label + '.' + base_26(hecke_orbit)
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], hecke_orbit)
    e['an'] = nf_elts_to_lists(Get_All_Hecke_Eigenvalues_Up_To(MAX_P+1, nf_lists_to_elements(e['lambda_p']),
                                                               nf_lists_to_elements(e['lambda_p_square']), e['weight']))
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
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein,
                    'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein}
       for sub in sub_funcs.keys():
           evs = sub_funcs[sub](k,j,e)
           for ev in evs:
               entry_sub = entry.copy()
               entry_sub.update(ev)
               entries.append(entry_sub)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
