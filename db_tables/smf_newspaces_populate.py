from sage.all import (is_square, is_squarefree, prime_range, factor, divisors, WeightedIntegerVectors)
from smf_lmfdb.db_tables.common_populate import make_space_label, entry_add_common_columns, table_reload, get_hecke, common_entry_values, base_26, write_data_from_files
from smf_lmfdb.db_tables.common_create_table import SUBSPACE_TYPES, HECKE_TYPES
from smf_lmfdb.db_tables.sage_functions import smf_dims_degree_2_level_1, smf_dims_degree_2_level_2, Hecke_Eigenvalues_Traces_Siegel_Eisenstein, Hecke_Eigenvalue_Traces_Klingen_Eisenstein, Hecke_Eigenvalue_Traces_Saito_Kurokawa, Hecke_Eigenvalue_Traces_Yoshida, num_forms_Siegel_Eisenstein, num_forms_Klingen_Eisenstein, num_forms_Saito_Kurokawa, num_forms_Yoshida
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import smf_dims_paramodular, Yoshida_lift_dim_orth, Yoshida_new_lift_dim_orth, Saito_Kurokawa_lift_dim
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_paramodular import Hecke_Eigenvalues_Traces_paramodular, num_forms_paramodular, parse_omf5
from lmfdb import db

def entry_add_columns(e, ext_data):
    e = entry_add_common_columns(e, ext_data)
    e['label'] = make_space_label(e)
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], 1)
    return e

def smf_level1_space(k,j,e):
    entry = smf_dims_degree_2_level_1(j,k,e)
    entry.update(common_entry_values(k,j,e+1,'P'))
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
    # manually handling (20,0,1) for demo
    if (k == 20) and (j == 0) and (e == 0):
        # manually fixing the (20,0) space for demonstration
        entry['num_forms'] += 1
        entry['cusp_G_lambda_p'] = [-840960,346935960,-73262366720,-5232247240500,2617414076964400,-724277370534455340,1427823701421564744,-83773835478688698980,14156088476175218899620,146957560176221097673720]
    return entry

def smf_level2_space(k,j):
    assert (j % 2 == 0) and (k >= 3)
    entry = smf_dims_degree_2_level_2(k,j)
    entry.update(common_entry_values(k,j,2,'P'))
    return entry

def num_level_raise(M,N):
    ret = 1
    for p, e in factor(N // M):
        ret *= WeightedIntegerVectors(e, [1,1,2]).cardinality()
    return ret

def count_old_G_forms(k,j,N):
    c = 0
    divs = [d for d in divisors(N) if d != N]
    for M in divs:
        _, dim_G_new = num_forms_paramodular(k,j,M)
        c += num_level_raise(M,N)*dim_G_new
    return c

def num_classical_minus_cusp_dim(k,N):
    '''
    Returns the dimension of S_k^{new, -}(N).
    We do it by querying the LMFDB for the total dimension and the
    dimension of the plus subspace

    Example #1 (p. 618):
    >>> [[classical_minus_cusp_dim(2*k-2,N) for k in range(3,12)]
    ... for N in [6,3,2,1]] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 1, 0, 1, 1, 1],
     [0, 0, 0, 1, 0, 1, 1, 1, 1],
     [0, 0, 0, 0, 0, 1, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 1, 0]]
    '''
    sgn = (-1)**(k//2-1)
    fs = db.mf_newforms.search({'weight' : k, 'level' : N, 'fricke_eigenval' : sgn,
                                'char_order' : 1},['label'])
    num_forms = len([f for f in fs])
    return num_forms

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
            # we temporarily go only up to a 1000 in paramodular
            if (k == 3) and (j == 0) and (not is_square(N)) and (N < 1000):
                traces, dim_G_new, al_dims_G = Hecke_Eigenvalues_Traces_paramodular(k,j,N)
                entry.update(traces)
                entry['num_forms'], dim_G_new = num_forms_paramodular(k,j,N)
                entry['num_forms'] += num_classical_minus_cusp_dim(2*k-2,N)
                entry['ALdims_G'] = al_dims_G
                entry['ALdims_P'] = [0 for al in al_dims_G]
                entry['ALdims'] = [0 for al in al_dims_G]
                divs = [d for d in divisors(N) if is_squarefree(d)]
                for i in range(len(divs)):
                    entry['ALdims_P'][i] = Saito_Kurokawa_lift_dim(k,N,al=divs[i])
                    entry['ALdims'][i] = entry['ALdims_G'][i] + entry['ALdims_P'][i]
                if not is_squarefree(N):
                    entry['new_cusp_G_dim'] = dim_G_new
                    entry['cusp_dim'] = cusp_dim
                    entry['old_cusp_G_dim'] = count_old_G_forms(k,j,N)
                    entry['cusp_G_dim'] = entry['old_cusp_G_dim'] + entry['new_cusp_G_dim']
                    entry['cusp_dim'] = entry['cusp_G_dim'] + entry['cusp_P_dim']
                    entry['new_cusp_dim'] = dim_G_new + entry['new_cusp_P_dim']
                    entry['old_cusp_dim'] = entry['cusp_dim'] - entry['new_cusp_dim']
                    assert entry['old_cusp_G_dim'] >= 0
                    assert entry['old_cusp_dim'] >= 0
                    assert entry['cusp_G_dim'] >= 0
                else:
                    assert entry['new_cusp_G_dim'] == dim_G_new
                    assert entry['old_cusp_G_dim'] == count_old_G_forms(k,j,N)
        entries.append(entry)
    return entries

def populate_smf_newspaces(triple_list):
    table = db.smf_newspaces
    aux_fname = "smf_lmfdb/db_tables/smf_newspaces_table.dat"
    entries = create_entries(triple_list)
    table_reload(table, entries, entry_add_columns, aux_fname, "newspaces")
    return

def update_cusp_dim(idx, space_folder, hecke_row_folder):
    space_file = open(space_folder + str(idx))
    space_data = eval(space_file.read())
    space_file.close()
    # if (space_data['cusp_dim'] != 'NULL'):
    #    return
    # assert (space_data['family'] == 'K')
    if (space_data['family'] != 'K'):
        return
    N = space_data['level']
    # Why do we have these when they are not squarefree?
    if (N >= 1000) or (N == 1):
        return
    # We add +2 for the case N = 2
    p = [x for x in prime_range(N+2) if N % x != 0][0]
    hecke_row_file = open(hecke_row_folder + "hecke_row_" + str(N) + "_" + str(p) + ".dat")
    hecke_row_data = eval(hecke_row_file.read())
    hecke_row_file.close()
    # our choice of p0
    exactly_divides = [p for p,e in factor(N) if e == 1]
    if (len(exactly_divides) == 0):
        return
    p0 = max(exactly_divides)
    yosh_all = Yoshida_lift_dim_orth(3,0,N,p0)
    yosh_new = Yoshida_new_lift_dim_orth(3,0,N,p0)
    # Subtract Eisenstein series and Yoshida lifts
    space_data['cusp_dim'] = sum([len(hecke_row_data[k]) for k in hecke_row_data.keys()])-1-yosh_all
    space_data['new_cusp_dim'] = space_data['new_cusp_G_dim'] + space_data['new_cusp_P_dim']
    space_data['old_cusp_dim'] = space_data['cusp_dim'] - space_data['new_cusp_dim']
    space_data['cusp_G_dim'] = space_data['cusp_dim'] - space_data['cusp_P_dim']
    space_data['old_cusp_G_dim'] = space_data['cusp_G_dim'] - space_data['new_cusp_G_dim']
    forms = parse_omf5(3,0,N)
    # Also updating the forms to include only new non-Yoshida
    space_data['num_forms'] = len([f for f in forms if f['aut_rep_type'] not in ['O', 'Y']])
    space_file = open(space_folder + str(idx), "w")
    space_file.write(str(space_data))
    space_file.close()
    return

def update_all_cusp_dim(space_folder, hecke_row_folder):
    table = db.smf_newspaces
    aux_fname = "smf_lmfdb/db_tables/smf_newspaces_table.dat"
    for idx in range(table.count()):
        print("updating cusp dim for idx = ", idx, " out of ", table.count())
        update_cusp_dim(idx, space_folder, hecke_row_folder)
    write_data_from_files(table, aux_fname, space_folder)
    table.reload(aux_fname, null="NULL")
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
