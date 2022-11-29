from smf_lmfdb.db_tables.common_populate import make_space_label, table_reload, get_hecke, common_entry_values, base_26, MAX_P, MAX_P_SQUARE
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.db_tables.sage_functions import Hecke_Eigenvalues_Siegel_Eisenstein, Hecke_Eigenvalues_Klingen_Eisenstein, Hecke_Eigenvalues_Saito_Kurokawa, Hecke_Eigenvalues_Yoshida, Get_All_Hecke_Eigenvalues_Up_To
from smf_lmfdb.db_tables.nf_elt import get_nf_basis, nf_lists_to_elements, nf_elts_to_lists
from smf_lmfdb.qExpansions.qexp_display import get_qexp_F20G

from lmfdb import db



def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    e['char_orbit_label'] = base_26(e['char_orbit_index'])
    space_label = make_space_label(e)
    dummy = e.pop('char_orbit_label')
    hecke_orbit = ext_data['num_forms']
    e['label'] = space_label + '.' + base_26(hecke_orbit)
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], hecke_orbit)
    # sometimes we don't have the Hecke eigenvalues
    if 'lambda_p' in e:
        basis, inv_basis = get_nf_basis(e)
        e['an'] = nf_elts_to_lists(Get_All_Hecke_Eigenvalues_Up_To(e['maxp']+1, nf_lists_to_elements(e['lambda_p'], basis),
                                                                   nf_lists_to_elements(e['lambda_p_square'], basis), e['weight']),
                                   inv_basis)
    return e

def create_entries(triple_list):
    entries = []
    for triple in triple_list:
        k,j,N = triple
        if (j % 2 == 1) or (k == 1):
            continue
        # right now we only have implemented forms for full level
        if (N > 1):
            continue
        for e in [0,1]:
            entry = common_entry_values(k,j,e+1)
            sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein,
                         'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein,
                         'cusp_P' : Hecke_Eigenvalues_Saito_Kurokawa,
                         'cusp_Y' : Hecke_Eigenvalues_Yoshida}
            for sub in sub_funcs.keys():
                evs = sub_funcs[sub](k,j,e)
                for ev in evs:
                    entry_sub = entry.copy()
                    entry_sub.update(ev)
                    entries.append(entry_sub)
        if (k == 20) and (j == 0) and (N == 1):
            entry = common_entry_values(k,j,N)
            entry['hecke_ring_rank'] = 1
            entry['hecke_ring_power_basis'] = True
            entry['hecke_ring_cyclotomic_generator'] = 0
            entry['field_poly'] = [0,1]
            entry['qexp'] = get_qexp_F20G()
            entries.append(entry)
    return entries

def populate_smf_hecke_nf(triple_list):
    table = db.smf_hecke_nf
    aux_fname = "smf_lmfdb/db_tables/smf_hecke_nf_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
