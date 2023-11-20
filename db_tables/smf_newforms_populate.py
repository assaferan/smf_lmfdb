from sage.all import (nth_prime, is_square)

from smf_lmfdb.db_tables.common_create_table import FAMILY_DICT
from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26, MAX_P, write_data_from_files, fill_nulls
from smf_lmfdb.db_tables.sage_functions import Hecke_Eigenforms_Siegel_Eisenstein, Hecke_Eigenforms_Klingen_Eisenstein, Hecke_Eigenforms_Saito_Kurokawa, Hecke_Eigenforms_Yoshida, Get_All_Hecke_Eigenvalues_Up_To, Get_All_Dirichlet_Coeffs_Up_To
from smf_lmfdb.qExpansions.qexp_display import get_qexp_display_F20G, get_qexp_display_E4, get_qexp_display_E6, get_qexp_display_Chi10, get_qexp_display_Chi12
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import Hecke_Eigenforms_paramodular

from lmfdb import db

def make_orbit_code(g, F, N, k, j, i, X):
    return g + (ord(F)<<8) + (N<<16) + (k<<32) + (j<<40) + ((i-1)<<48) + ((X-1)<<56)

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
    if ('trace_lambda_p_square' in e) and (e['trace_lambda_p_square'] != 'NULL'):
        max_p = min(MAX_P, nth_prime(len(e['trace_lambda_p'])))
        bad_ps = prime_divisors(e['level'])
        eps = { p : 1 for p in bad_ps }
        if len(bad_ps) > 0: 
            eps[max(bad_ps)] = -1
            assert e['atkin_lehner_eigenvals'][-1][0] == max(bad_ps)
            e['atkin_lehner_eigenvals'][-1][1] *= (-1)
            if (e['atkin_lehner_string'][-1] == '+'):
                e['atkin_lehner_string'][-1] = '-'
            else:
                e['atkin_lehner_string'][-1] = '+'
            e['fricke_eigenval'] = reduce(lambda x,y:x*y, [x[1] for x in e['atkin_lehner_eigenvals']], 1)
        e['trace_display'] = e['trace_lambda_p'][:4]
        # e['traces'] = Get_All_Hecke_Eigenvalues_Up_To(max_p+1, e['trace_lambda_p'], e['trace_lambda_p_square'], e['weight'])
        e['traces'] = Get_All_Dirichlet_Coeffs_Up_To(max_p+1, e['trace_lambda_p'], e['trace_lambda_p_square'], e['weight'], e['level'], eps)
    return e

def write_temp_entry(entry, folder, ext_data, table):
    if entry['level'] == 1:
        families = FAMILY_DICT.values()
    else:
        families = [entry['family']]
    num_entries = 0
    for family in families:
        entry['family'] = family
        e = entry.copy()
        e = entry_add_columns(e, ext_data)
        e['id'] = entry['id'] + num_entries
        e = fill_nulls(e, table)
        fname = "smf_lmfdb/db_tables/data/" + folder + "/" + str(e['id'])
        f = open(fname, 'w')
        f.write(str(e))
        f.close()
        num_entries += 1
    return num_entries

def create_entries(triple_list, folder, table):
    space_num_forms = {}
    entries = []
    idx = 0
    for triple in triple_list:
        print("creating newform entry for triple", triple)
        k,j,N = triple
        if (j % 2 == 1) or (k == 1):
            continue
        # right now we only have implemented forms for full level
        # for paramodular we stop at 1000 at the moment
        if (not is_square(N)) and (N < 1000):
            if (((k == 3) and (j == 0)) or
                ((k == 3) and (j == 2) and (N == 19)) or
                ((k == 4) and (j == 0) and (N == 31))):
                entry = common_entry_values(k,j,N, 'K')
                forms = Hecke_Eigenforms_paramodular(k,j,N)
                forms = sorted(forms, key=lambda f : [f['dim']] + f['trace_lambda_p'])
                for f in forms:
                    entry_sub = entry.copy()
                    entry_sub.update(f)
                    entries.append(entry_sub)
                    space_label = make_space_label(entry_sub, False)
                    space_num_forms[space_label] = space_num_forms.get(space_label,0) + 1
                    idx += write_temp_entry(entry_sub, folder,  {'id' : idx, 'num_forms' : space_num_forms[space_label]}, table)
            continue
                
        for e in [0,1]:
            entry = common_entry_values(k,j,e+1, 'P')
            sub_funcs = {'eis_F' : Hecke_Eigenforms_Siegel_Eisenstein,
                         'eis_Q' : Hecke_Eigenforms_Klingen_Eisenstein,
                         'cusp_P': Hecke_Eigenforms_Saito_Kurokawa,
                         'cusp_Y': Hecke_Eigenforms_Yoshida}
            for sub in sub_funcs.keys():
                forms = sub_funcs[sub](k,j,e)
                for f in forms:
                    entry_sub = entry.copy()
                    entry_sub.update(f)
                    # adding several qexpansions for demonstration purposes
                    if (j == 0) and (N == 1) and (sub == 'eis_F'):
                        if (k == 4):
                            entry_sub['qexp_display'] = get_qexp_display_E4()
                        if (k == 6):
                            entry_sub['qexp_display'] = get_qexp_display_E6()
                    if (j == 0) and (N == 1) and (sub == 'cusp_P'):
                        if (k == 10):
                            entry_sub['qexp_display'] = get_qexp_display_Chi10()
                        if (k == 12):
                            entry_sub['qexp_display'] = get_qexp_display_Chi12()
                    entries.append(entry_sub)
                    space_label = make_space_label(entry_sub, False)
                    space_num_forms[space_label] = space_num_forms.get(space_label,0) + 1
                    write_temp_entry(entry_sub, folder,  {'id' : idx, 'num_forms' : space_num_forms[space_label]}, table)
                    idx += 1
        # adding for demonstration a single function
        if (k == 20) and (j == 0) and (N == 1):
            entry = common_entry_values(k,j,N, 'P')
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
            entry['trace_lambda_p'] = [-840960,346935960,-73262366720,-5232247240500,2617414076964400,-724277370534455340,1427823701421564744,-83773835478688698980,14156088476175218899620,146957560176221097673720]
            entries.append(entry)
            space_label = make_space_label(entry, False)
            space_num_forms[space_label] = space_num_forms.get(space_label,0) + 1
            write_temp_entry(entry, folder, {'id' : idx, 'num_forms' : space_num_forms[space_label]}, table)
            idx += 1
    return entries

def populate_smf_newforms(triple_list):
    table = db.smf_newforms
    aux_fname = "smf_lmfdb/db_tables/smf_newforms_table.dat"
    entries = create_entries(triple_list, "newforms", table)
    table_reload(table, entries, entry_add_columns, aux_fname, "newforms")
    return

def update_yoshida(idx, forms_folder, space_num_forms):
    form_file = open(forms_folder + str(idx))
    form_data = eval(form_file.read())
    form_file.close()
    # Changing Yoshida to be Siegel
    if (form_data['family'] == 'K') and (form_data['aut_rep_type'] == 'Y'):
        form_data['family'] = 'S'
        form_data['space_label'] = make_space_label(form_data)
        
    space_num_forms[form_data['space_label']] = space_num_forms.get(form_data['space_label'],0)+1
    form_data['hecke_orbit'] = space_num_forms[form_data['space_label']]
    ret = {'old_label' : form_data['label'], 'old_orbit_code' : form_data['hecke_orbit_code']}
    form_data['label'] = form_data['space_label'] + '.' + base_26(form_data['hecke_orbit'])
    form_data['hecke_orbit_code'] = make_orbit_code(form_data['degree'], form_data['family'], form_data['level'], form_data['weight'][0], form_data['weight'][1], form_data['char_orbit_index'], form_data['hecke_orbit'])
    ret.update({'new_label' : form_data['label'], 'new_orbit_code' : form_data['hecke_orbit_code']})
    form_file = open(forms_folder + str(idx), "w")
    form_file.write(str(form_data))
    form_file.close()
    return ret, space_num_forms

def update_all_yoshida(forms_folder):
    table = db.smf_newforms
    aux_fname = "smf_lmfdb/db_tables/smf_newforms_table.dat"
    space_num_forms = {}
    track_changes = []
    for idx in range(table.count()):
        print("updating Yoshida lifts for idx =  ", idx, "out of ", table.count())
        yosh, space_num_forms = update_yoshida(idx, forms_folder, space_num_forms)
        track_changes.append(yosh)
    label_dict = {t['old_label'] : t['new_label'] for t in track_changes}
    orbit_code_dict = {t['old_orbit_code'] : t['new_orbit_code'] for t in track_changes}
    write_data_from_files(table, aux_fname, forms_folder)
    table.reload(aux_fname, null="NULL")
    return label_dict, orbit_code_dict
