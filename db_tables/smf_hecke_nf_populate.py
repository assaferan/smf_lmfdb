from sage.all import (is_square)
from smf_lmfdb.db_tables.common_populate import make_space_label, table_reload, get_hecke, common_entry_values, base_26, MAX_P, MAX_P_SQUARE
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.db_tables.sage_functions import Hecke_Eigenvalues_Siegel_Eisenstein, Hecke_Eigenvalues_Klingen_Eisenstein, Hecke_Eigenvalues_Saito_Kurokawa, Hecke_Eigenvalues_Yoshida, Get_All_Hecke_Eigenvalues_Up_To
from smf_lmfdb.db_tables.nf_elt import get_nf_basis, nf_lists_to_elements, nf_elts_to_lists
from smf_lmfdb.qExpansions.qexp_display import get_qexp_F20G, get_qexp_E4, get_qexp_E6, get_qexp_Chi10, get_qexp_Chi12
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import Hecke_Eigenvalues_paramodular

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
    if 'lambda_p_square' in e:
        basis, inv_basis = get_nf_basis(e, False)
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
        # at the moment we stop at 100 for paramodulars
        if (not is_square(N)) and (N <= 1000):
            if (k == 3) and (j == 0):
                entry = common_entry_values(k,j,N,'K')
                evs = Hecke_Eigenvalues_paramodular(k,j,N)
                for ev in evs:
                    entry_sub = entry.copy()
                    entry_sub.update(ev)
                    entries.append(entry_sub)
            continue
        # right now we only have implemented forms for full level in the principal series case
        for e in [0,1]:
            entry = common_entry_values(k,j,e+1,'P')
            sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein,
                         'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein,
                         'cusp_P' : Hecke_Eigenvalues_Saito_Kurokawa,
                         'cusp_Y' : Hecke_Eigenvalues_Yoshida}
            for sub in sub_funcs.keys():
                evs = sub_funcs[sub](k,j,e)
                for ev in evs:
                    entry_sub = entry.copy()
                    entry_sub.update(ev)
                    if (j == 0) and (N == 1) and (sub == 'eis_F') and (e == 0):
                        if (k == 4):
                            entry_sub['qexp'] = get_qexp_E4()
                        if (k == 6):
                            entry_sub['qexp'] = get_qexp_E6()
                    if (j == 0) and (N == 1) and (sub == 'cusp_P') and (e == 0):
                        if (k == 10):
                            entry_sub['qexp'] = get_qexp_Chi10()
                        if (k == 12):
                            entry_sub['qexp'] = get_qexp_Chi12()
                    entries.append(entry_sub)
        if (k == 20) and (j == 0) and (N == 1):
            entry = common_entry_values(k,j,N,'P')
            entry['hecke_ring_rank'] = 1
            entry['hecke_ring_power_basis'] = True
            entry['hecke_ring_cyclotomic_generator'] = 0
            entry['field_poly'] = [0,1]
            entry['qexp'] = get_qexp_F20G()
            entry['maxp'] = 19
            entry['lambda_p'] = [-840960,346935960,-73262366720,-5232247240500,2617414076964400,-724277370534455340,1427823701421564744,-83773835478688698980,14156088476175218899620,146957560176221097673720]
            entries.append(entry)
    return entries

def populate_smf_hecke_nf(triple_list):
    table = db.smf_hecke_nf
    aux_fname = "smf_lmfdb/db_tables/smf_hecke_nf_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
