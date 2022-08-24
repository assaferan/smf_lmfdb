import os
cwd = os.getcwd()
os.chdir("../Dimension_formulas")
load('dimformula_smf_degree2_level_1.sage')
os.chdir("../../lmfdb")
from lmfdb import db
os.chdir(cwd)

aux_fname = "smf_newspaces_table.dat"

def make_label(e):
    return '.'.join([str(x) for x in [e['degree'], e['type'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_label']]])

def make_id(e):
    return hash(e['label'])

def write_data(entries):
    keys = [k for k in db.smf_newspaces.col_type.keys()]
    types = [db.smf_newspaces.col_type[k] for k in keys]
    column_names = '|'.join(keys)
    column_types = '|'.join(types)
    e_data = []
    for e in entries:
        e['char_orbit_label'] = chr(ord('a') + e['char_orbit'])
        e['char_order'] = 1 + e['char_orbit']
        e['char_degree'] = 1
        e['conrey_indexes'] = [1+2*e['char_orbit']]
        e['label'] = make_label(e)
        e['id'] = make_id(e)
        e_datum = '|'.join([str(e[k]) for k in keys])
        e_datum = e_datum.replace('[', '{').replace(']','}')
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

# triple_list consists of triples (k,j,e) of weight and character
def generate_dimensions(triple_list):
    entries = []
    for triple in triple_list:
       k = triple[0]
       j = triple[1]
       e = triple[2]
       entries.append(dim_splitting_smf_degree_2_level_1(j,k,e))
    table_reload(entries)
    return