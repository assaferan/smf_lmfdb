from smf_lmfdb.db_tables.common_create_table import generate_table

def generate_column_types():
    col_type = {}
    col_type['label'] = 'text'
    col_type['hecke_orbit_code'] = 'bigint'
    col_type['hecke_ring_rank'] = 'integer'
    col_type['hecke_ring_power_basis'] = 'boolean'
    col_type['hecke_ring_cyclotomic_generator'] = 'integer'
    col_type['hecke_ring_character_values'] = 'jsonb'
    col_type['maxp'] = 'integer'
    col_type['maxp_square'] = 'integer'
    col_type['an'] = 'jsonb'
#    col_type['ap'] = 'jsonb'
    col_type['lambda_p'] = 'jsonb'
    col_type['lambda_p_square'] = 'jsonb'
    col_type['lambda_p_square_0'] = 'jsonb'
    col_type['lambda_p_square_1'] = 'jsonb'
    col_type['lambda_p_square_2'] = 'jsonb'
    col_type['field_poly'] = 'numeric[]'
    col_type['hecke_ring_numerators'] = 'numeric[]'
    col_type['hecke_ring_denominators'] = 'numeric[]'
    col_type['hecke_ring_inverse_numerators'] = 'numeric[]'
    col_type['hecke_ring_inverse_denominators'] = 'numeric[]'
    col_type['degree'] = 'smallint'
    col_type['family'] = 'text'
    col_type['level'] = 'integer'
    col_type['weight'] = 'smallint[]'
    col_type['char_orbit_index'] = 'smallint[]'
    return col_type

def generate_column_desc():
    col_desc = {}
    col_desc['label'] = 'Label g.C.N.w.a.x of this newform'
    col_desc['hecke_orbit_code'] = 'encoding of the tuple (g.C.N.w.i.x) into 64 bits, used in eigenvalue tables.  g + (ord(C)<<8) + (N<<12) + (k<<20) + (j<<28) + ((i-1)<<36) + ((X-1)<<52).'
    col_desc['hecke_ring_rank'] = 'rank of Hecke ring as a free Z-module = dimension of newform = degree of field_lpoly'
    col_desc['hecke_ring_power_basis'] = 'True if we are using the power basis specified by field_poly as the basis for the hecke ring'
    col_desc['hecke_ring_cyclotomic_generator'] = 'zero or an integer m suth that an and ap are encoded as sparse integer polynomials in zeta_m (typically zeta_m is a root of the field poly but not necessarily)'
    col_desc['hecke_ring_character_values'] = 'list of pairs [[m1,[a11,...,a1n]],...[mr,[a1r,...,arn]]] where mi are generators of (Z/NZ)* and [ai1,...,ain] is the value of chi(mi) expressed in terms of the Hecke ring basis or in cyclotomic representation [[c,e]] encoding c x zeta_m^e where m is hecke_ring_cyclotomic_generator'
    col_desc['maxp'] = 'largest prime for which a_p appears in the list ap'
    col_desc['maxp_square'] = 'largest prime p for which lambda_p_square appears in the list lambdap_square'
    col_desc['an'] = 'list of lists encoding the coefficients of the L-series a_1,a_2,...,a_100 either as a linear combination of the basis specified in the orbit table or as a list of pairs [ci,ei] encoding prod_i c_i*zeta_m^e_i, where m it the value of hecke_ring_cyclotomic_generator'
    col_desc['lambda_p'] = 'list of lists encoding Hecke eigenvalues a_p (same format as a_n) for primes p up to maxp'
    col_desc['lambda_p_square'] = 'list of lists encoding Hecke eigenvalues of a_{p^2} (same format as a_n) for primes p up to maxp_square'
    col_desc['lambda_p_square_0'] = 'list of lists encoding Hecke eigenvalues of a_{p^2,0} (same format as a_n) for primes p up to maxp_square'
    col_desc['lambda_p_square_1'] = 'list of lists encoding Hecke eigenvalues of a_{p^2,1} (same format as a_n) for primes p up to maxp_square'
    col_desc['lambda_p_square_2'] = 'list of lists encoding Hecke eigenvalues of a_{p^2,2} (same format as a_n) for primes p up to maxp_square'
    col_desc['field_poly'] = 'list of integer coefficients of a defining polynomial for the Hecke field'
    col_desc['hecke_ring_numerators'] = 'List of lists of integers giving numerators of the basis for the hecke ring in terms of the power basis (if hecke_ring_power_basis is false)'
    col_desc['hecke_ring_denominators'] = 'List of integers giving denominators of the basis for the hecke ring in terms of the power basis (if hecke_ring_power_basis is false)'
    col_desc['hecke_ring_inverse_numerators'] = 'List of lists of integers giving numerators of the inverse basis that represents power basis in terms of the basis for the hecke ring'
    col_desc['hecke_ring_inverse_denominators'] = 'List of integers giving denominators of the inverse basis that represents power basis in terms of the basis for the hecke ring'
    col_desc['degree'] = 'Degree g (automorphic with repsect to Sp(2g, Q))'
    col_desc['family'] = "Family of arithmetic subgroups ('K' = paramodular, 'S' = Siegel, 'P' =  principal)"
    col_desc['level'] = 'Level N in the family'
    col_desc['weight'] = 'Weight (k,j)'
    col_desc['char_orbit_index'] = 'ordinal i identifying the Galois orbit of the characterof this newform (base26 encoded in the newform label / character orbit label)'
    return col_desc
    
def create_table_smf_hecke_nf():
    table_name = "smf_hecke_nf"
    table_desc = "Hecke eigenvalues of Siegel modular forms over number fields"
    generate_table(table_name, table_desc,
                   generate_column_types, generate_column_desc)
    return
