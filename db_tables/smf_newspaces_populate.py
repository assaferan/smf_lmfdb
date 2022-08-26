import os
cwd = os.getcwd()
os.chdir("../Dimension_formulas")
load('dimformula_smf_degree2_level_1.sage')
os.chdir("../Hecke_Eigenvalues/Siegel_Eisenstein_series/")
load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')
os.chdir("../Klingen_Eisenstein_series/")
load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')
os.chdir(cwd)
os.chdir("../../lmfdb")
from lmfdb import db
os.chdir(cwd)

aux_fname = "smf_newspaces_table.dat"

def make_label(e):
    return '.'.join([str(x) for x in [e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_label']]])

def write_data(entries):
    keys = [k for k in db.smf_newspaces.col_type.keys()]
    types = [db.smf_newspaces.col_type[k] for k in keys]
    column_names = '|'.join(keys)
    column_types = '|'.join(types)
    e_data = []
    for id,e in enumerate(entries):
        e['char_orbit_label'] = chr(ord('a') + e['char_orbit'])
        e['char_order'] = 1 + e['char_orbit']
        if (e['char_order'] == 2):
           e['level'] = 2
        e['char_degree'] = 1
        e['conrey_indexes'] = [1+2*e['char_orbit']]
        e['label'] = make_label(e)
        e['id'] = id
        e_datum = '|'.join([str(e[k]) for k in keys])
        e_datum = e_datum.replace('[', '{').replace(']','}')
        e_datum = e_datum.replace('(', '{').replace(')','}')
        e_data.append(e_datum)
    write_data = "\n".join([column_names, column_types, ""] + e_data)
    f = open(aux_fname, "w")
    f.write(write_data)
    f.close()
    return

def table_reload(entries):
    write_data(entries)
    db.smf_newspaces.reload(aux_fname, null="")
    db.smf_newspaces.cleanup_from_reload()
    return

def get_hecke(func,deg,hecke_type,j,k,e,prime_bound=200):
    prime_bound = prime_bound^(1/deg)
    vals = func(k,j,e)['lambda_' + hecke_type]
    return [vals[p] for p in prime_range(prime_bound)]

# triple_list consists of triples (k,j,e) of weight and character
def generate_dimensions(triple_list):
    entries = []
    for triple in triple_list:
       k = triple[0]
       j = triple[1]
       e = triple[2]
       entry = dim_splitting_smf_degree_2_level_1(j,k,e)
       hecke_types = {1 : ['p'],
       		      2 : ['p_square', 'p_square_0', 'p_square_1', 'p_square_2']}
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein_Series_All,
       		    'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac}
       for sub in sub_funcs.keys():
          for deg in hecke_types.keys():
       	      for hecke_type in hecke_types[deg]:
              	  entry[sub + '_lambda_' + hecke_type] = get_hecke(sub_funcs[sub],deg,hecke_type,j,k,e)
       entries.append(entry)
    table_reload(entries)
    return