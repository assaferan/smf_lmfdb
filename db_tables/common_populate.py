from sage.all import prime_range, is_prime, is_prime_power, is_square, is_squarefree, prime_divisors, radical

from lmfdb import db

def make_space_label(e, label=True):
    last_key = 'char_orbit'
    if label:
        last_key += '_label'
    return '.'.join([str(x) for x in [e['degree'], e['family'], e['level']] + list(e['weight']) +[e[last_key]]])

def base_26(index):
    res = index - 1
    letter = chr(ord('a') + (res % 26))
    ret = letter
    res //= 26
    while (res > 0):
        letter = chr(ord('a') + (res % 26))
        ret = letter + ret
        res //= 26
    return ret

# TODO - the generating functions should return that to us
def char_order(char_orbit_index):
    # right now just working for Fabien's code (index = 0 or 1)
    return char_orbit_index

def conrey_indexes(char_orbit_index):
    # right now just working for Fabien's code (index = 1 or 2)
    return [2*char_orbit_index-1]

def common_entry_values(k,j,e):
    entry = {}
    entry['degree'] = 2
    entry['family'] = 'S'
    entry['weight'] = [k,j]
    entry['char_orbit'] = e
    entry['level'] = 1
    return entry

def entry_add_common_columns(e, ext_data):
    # we do that because the non-trivial level have level 2
    if (e['char_orbit'] == 1):
        e['family'] = 'P'
        e['level'] = 2
#        e['char_conductor'] = 4
#    else:
    e['char_orbit'] = 0
    e['char_conductor'] = 1
    e['char_orbit_index'] = e['char_orbit'] + 1
    e['prim_orbit_index'] = e['char_orbit_index']
    e['char_orbit_label'] = base_26(e['char_orbit_index'])
    e['char_order'] = char_order(e['char_orbit_index'])
    # In our cases the degree is always 1
    e['char_degree'] = 1
    e['conrey_indexes'] = conrey_indexes(e['char_orbit_index'])
    e['id'] = ext_data['id']
    e['level_radical'] = radical(e['level'])
    e['char_parity'] = 3-2*e['char_orbit_index']
    e['char_is_minimal'] = True
    e['char_is_real'] = True
    e['char_values'] = "NULL"
    e['level_is_prime'] = is_prime(e['level'])
    e['level_is_prime_power'] = is_prime_power(e['level'])
    e['level_is_square'] = is_square(e['level'])
    e['level_is_squarefree'] = is_squarefree(e['level'])
    e['level_primes'] = prime_divisors(e['level'])
    return e

def write_data(table, entries, entry_postprocess, aux_fname):
    keys = [k for k in table.col_type.keys()]
    types = [table.col_type[k] for k in keys]
    column_names = '|'.join(keys)
    column_types = '|'.join(types)
    e_data = []
    space_num_forms = {}
    for i,e in enumerate(entries):
        space_label = make_space_label(e, False)
        space_num_forms[space_label] = space_num_forms.get(space_label,0) + 1
        e = entry_postprocess(e, {'id' : i,
                                  'num_forms' : space_num_forms[space_label]})
        e_datum = '|'.join([str(e[k]) for k in keys])
        e_datum = e_datum.replace('[', '{').replace(']','}')
        e_datum = e_datum.replace('(', '{').replace(')','}')
        e_data.append(e_datum)
    write_data = "\n".join([column_names, column_types, ""] + e_data)
    f = open(aux_fname, "w")
    f.write(write_data)
    f.close()
    return

def table_reload(table, entries, entry_postprocess, aux_fname):
    write_data(table, entries, entry_postprocess, aux_fname)
    table.reload(aux_fname, null="NULL")
    table.cleanup_from_reload()
    return

def get_hecke(func,deg,hecke_type,j,k,e,prime_bound=200):
    prime_bound = prime_bound**(1/deg)
    vals = func(k,j,e)['lambda_' + hecke_type]
    return [vals[p] for p in prime_range(prime_bound)]
