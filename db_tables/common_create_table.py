import os
cwd = os.getcwd()
os.chdir("../../lmfdb")
from lmfdb import db
os.chdir(cwd)

FAMILY_DICT = {
    'paramodular' : 'K',
    'Siegel'      : 'S',
    'principal'   : 'P'
}

hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']

subspace_types = ['eis_F', 'eis_Q']

def generate_common_column_types():
    col_type = {}
    col_type['degree'] = 'smallint'
    col_type['family'] = 'text'
    col_type['level'] = 'integer'
    col_type['level_radical'] = 'integer'
    col_type['weight'] = 'smallint[]'
    col_type['char_orbit'] = 'smallint'
    col_type['char_orbit_index'] = 'smallint'
    col_type['char_orbit_label'] = 'text'
    col_type['char_order'] = 'integer'
    col_type['char_degree'] = 'integer'
    col_type['char_parity'] = 'smallint'
    col_type['char_is_minimal'] = 'boolean'
    col_type['char_conductor'] = 'integer'
    col_type['prim_orbit_index'] = 'smallint'
    col_type['char_is_real'] = 'boolean'
    col_type['char_values'] = 'jsonb'
    col_type['conrey_indexes'] = 'integer[]'
    col_type['level_is_prime'] = 'boolean'
    col_type['level_is_prime_power'] = 'boolean'
    col_type['level_is_square'] = 'boolean'
    col_type['level_is_squarefree'] = 'boolean'
    col_type['level_primes'] = 'integer[]'
    col_type['label'] = 'text'
    return col_type

def generate_search_columns(col_type):
    search_col = {}
    types = set(col_type.values())
    for t in types:
        search_col[t] = [k for k in col_type.keys() if col_type[k] == t]
    return search_col

def generate_common_col_desc():
    col_desc = {}
    col_desc['degree'] = 'Degree g (automorphic with repsect to Sp(2g, Q))'
    col_desc['family'] = "Family of arithmetic subgroups ('K' = paramodular, 'S' = Siegel, 'P' =  principal)"
    col_desc['level'] = 'Level N in the family'
    col_desc['level_is_prime'] = 'true if N is prime (1 is not prime)'
    col_desc['level_is_prime_power'] = 'true if N is a prime power (1 is not a prime power)'
    col_desc['level_is_square'] = 'true if N is square'
    col_desc['level_is_squarefree'] = 'true if N is squarefree'
    col_desc['level_primes'] = 'sorted list of prime divisors of N'
    col_desc['level_radical'] = 'product of prime divisors of N'
    col_desc['weight'] = 'Weight (k,j)'
    col_desc['char_parity'] = '1 for even, -1 for odd'
    col_desc['char_is_minimal'] =  "true if the character chi is {{KNOWL('character.dirichlet.minimal','minimal')}}"
    col_desc['char_orbit'] = 'ordinal i identifying the Galois orbit of the characterof this newform (base26 encoded in the newform label / character orbit label)'
    col_desc['char_orbit_index'] = 'ordinal i identifying the Galois orbit of the characterof this newform (base26 encoded in the newform label / character orbit label)'
    col_desc['char_orbit_label'] = 'base26-encoding of char_orbit-1'
    col_desc['char_order'] = 'the order of the character chi'
    col_desc['char_degree'] = 'Degree of the (cyclotomic) character field'
    col_desc['char_conductor'] = 'Conductor of the Dirichlet character chi of this newform'
    col_desc['char_is_real'] = 'true if the character takes only real values (trivial or quadratic)'
    col_desc['char_values'] = 'quadruple <N,n,u,v> where N is the level, n is the order of the character, u is a list of generators for the unit group of Z/NZ, and v is a corresponding list of integers for which chi(u[i]) = zeta_n^v[i]'
    col_desc['prim_orbit_index'] = 'char_orbit for the primitive version of this character'
    col_desc['conrey_indexes'] = 'Sorted list of Conrey indexes of characters in this Galois orbit'
    return col_desc

def generate_table(table_name, table_desc, col_type_func, col_desc_func):
    if table_name in db.tablenames:
       db.drop_table(table_name)
    col_type = col_type_func()
    search_col = generate_search_columns(col_type)
    col_desc = col_desc_func()
    db.create_table(table_name, search_col, 'label', table_desc, col_desc)
    return
