load('../Dimension_formulas/DimFormulaSMFVectorValuedLevel1WithoutChara\
cter.sage')
import os
cwd = os.getcwd()
os.chdir("../../lmfdb")
from lmfdb import db
os.chdir(cwd)

aux_fname = "smf_newspaces_table.dat"

def make_label(e):
    return '.'.join([str(x) for x in [e['degree'], e['type'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit']]])

def make_id(e):
    return hash(e['label'])

def write_data(entries):
    keys = [k for k in db.smf_newspaces.col_type.keys()]
    types = [db.smf_newspaces.col_type[k] for k in keys]
    column_names = '|'.join(keys)
    column_types = '|'.join(types)
    e_data = []
    for e in entries:
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

def generate_dimensions(weight_list):
    entries = []
    for weight in weight_list:
       k = weight[0]
       j = weight[1]
       entries.append(dim_splitting_VV_All_weight(k,j))
    table_reload(entries)
    return