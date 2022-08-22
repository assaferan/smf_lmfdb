from lmfdb import db

aux_fname = "smf_newspaces_table.dat"

def make_label(e):
    return '.'.join([e['degree'], e['type'], e['level'], e['weight'][0], e['weight'][1]])

def write_data(entries):
    keys = [k for k in db.smf_spaces.col_type.keys()]
    types = [db.smf_spaces.col_type[k] for k in keys]
    column_names = '|'.join(keys)
    column_types = '|'.join(types)
    e_data = []
    for e in entries:
        e['label'] = make_label(e)
        e_datum = '|'.join([str(e[k]) for k in keys])
        e_datum = e_datum.replace('[', '{').replace(']','}')
        e_data.append(e_datum)
    write_data = "\n".join([column_names, column_types, ""] + e_data)
    f = open(aux_fname, "w")
    f.write(write_data)
    f.close()
    return

def table_reload(entries):    
    db.smf_newforms.reload(aux_fname, null="")
    db.smf_newforms.cleanup_from_reload()
    return
