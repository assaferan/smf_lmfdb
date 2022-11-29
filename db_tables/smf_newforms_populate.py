from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26, MAX_P
from smf_lmfdb.db_tables.sage_functions import Hecke_Eigenforms_Siegel_Eisenstein, Hecke_Eigenforms_Klingen_Eisenstein, Hecke_Eigenforms_Saito_Kurokawa, Hecke_Eigenforms_Yoshida, Get_All_Hecke_Eigenvalues_Up_To
from smf_lmfdb.qExpansions.qexp_display import get_qexp_display_F20G, get_qexp_display_E4, get_qexp_display_E6, get_qexp_display_Chi10, get_qexp_display_Chi12,

from lmfdb import db

def make_orbit_code(g, F, N, k, j, i, X):
    return g + (ord(F)<<8) + (N<<12) + (k<<20) + (j<<28) + ((i-1)<<36) + ((X-1)<<52)

def entry_add_columns(e, ext_data):
    e = entry_add_common_columns(e, ext_data)
    e['space_label'] = make_space_label(e)
    e['hecke_orbit'] = ext_data['num_forms']
    e['label'] = e['space_label'] + '.' + base_26(e['hecke_orbit'])
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], e['hecke_orbit'])
    # for now we don't populate these fields
    e['trace_hash'] = 'NULL'
    e['analytic_rank'] = 'NULL'
    e['analytic_rank_proved'] = False
    if 'qexp_display' not in e:
        e['qexp_display'] = 'NULL'
    e['embedded_related_objects'] = []
    # sometimes we just have q-expansions and no hecke eigenvalues
    if 'trace_lambda_p' in e:
        e['trace_display'] = e['trace_lambda_p'][:4]
        e['traces'] = Get_All_Hecke_Eigenvalues_Up_To(MAX_P+1, e['trace_lambda_p'], e['trace_lambda_p_square'], e['weight'])
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
            sub_funcs = {'eis_F' : Hecke_Eigenforms_Siegel_Eisenstein,
                         'eis_Q' : Hecke_Eigenforms_Klingen_Eisenstein,
                         'cusp_P': Hecke_Eigenforms_Saito_Kurokawa,
                         'cusp_Y': Hecke_Eigenforms_Yoshida}
            for sub in sub_funcs.keys():
                forms = sub_funcs[sub](k,j,e)
                for f in forms:
                    entry_sub = entry.copy()
                    entry_sub.update(f)
                    if (j == 0) and (N == 1) and (sub == 'eis_F'):
                        if (k == 4):
                            entry['qexp_display'] = get_qexp_display_E4()
                        if (k == 6):
                            entry['qexp_display'] = get_qexp_display_E6()
                    if (j == 0) and (N == 1) and (sub == 'cusp_P'):
                        if (k == 10):
                            entry['qexp_display'] = get_qexp_display_Chi10()
                        if (k == 12):
                            entry['qexp_display'] = get_qexp_display_Chi12()
                    entries.append(entry_sub)
        # adding for demonstration a single function
        if (k == 20) and (j == 0) and (N == 1):
            entry = common_entry_values(k,j,N)
            entry['is_cuspidal'] = True
            entry['aut_rep_type'] = 'G'
            entry['qexp_display'] = get_qexp_display_F20G()
            entry['dim'] = 1
            entry['nf_label'] = '1.1.1.1'
            entry['hecke_ring_index'] = 1
            entry['field_poly_is_cyclotomic'] = False
            entry['field_poly_root_of_unity'] = 0
            entry['related_objects'] = []
            entry['field_poly_is_real_cyclotomic'] = False
            entry['hecke_ring_index_proved'] = True
            entry['hecke_ring_generator_nbound'] = 1
            entry['field_disc'] = 1
            entry['field_disc_factorization'] = []
            entry['field_poly'] = [0,1]
            entry['hecke_ring_index_factorization'] = []
            entry['relative_dim'] = 1
            entry['is_polredabs'] = True
            entries.append(entry)
    return entries

def populate_smf_newforms(triple_list):
    table = db.smf_newforms
    aux_fname = "smf_lmfdb/db_tables/smf_newforms_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
