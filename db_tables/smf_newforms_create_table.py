from smf_lmfdb.db_tables.common_create_table import generate_common_column_types, generate_common_col_desc, generate_table, HECKE_TYPES

def generate_column_types():
    col_type = generate_common_column_types()
    col_type['aut_rep_type'] = 'text'
    col_type['space_label'] = 'text'
    col_type['hecke_orbit'] = 'integer'
    col_type['hecke_orbit_code'] = 'bigint'
    col_type['dim'] = 'integer'
    col_type['relative_dim'] = 'integer'
    col_type['nf_label'] = 'text'
    col_type['field_disc'] = 'numeric'
    col_type['field_poly_is_cyclotomic'] = 'boolean'
    col_type['field_poly'] = 'numeric[]'
    col_type['field_poly_is_real_cyclotomic'] = 'boolean'
    col_type['field_poly_root_of_unity'] = 'integer'
    col_type['field_disc_factorization'] = 'numeric[]'
    col_type['trace_hash'] = 'bigint'
    col_type['analytic_rank'] = 'smallint'
    col_type['analytic_rank_proved'] = 'boolean'
    col_type['qexp_display'] = 'text'
    col_type['related_objects'] = 'text[]'
    col_type['embedded_related_objects'] = 'text[]'
    col_type['trace_display'] = 'numeric[]'
    col_type['traces'] = 'numeric[]'
    col_type['is_cuspidal'] = 'boolean'
    col_type['lift_type'] = 'text'
    for hecke_type in HECKE_TYPES:
    	col_type['trace_lambda_'+hecke_type] = 'numeric[]'
    return col_type

def generate_column_desc():
    col_desc = generate_common_col_desc()
    col_desc['label'] = 'Label g.C.N.w.a.x of this newform'
    col_desc['aut_rep_type'] = 'Type of the automorphic representation corresponding to this form - one of (F, B, P, Q, Y, G) according to the classification of Arthur parameters by Schmidt'
    col_desc['space_label'] = 'label g.C.N.w.a of the newspace containing this newform'
    col_desc['hecke_orbit'] = '(X) An integer that is encoded into x in the label via 1=a, 2=b, 26=z, 27=ba, 28=bb.  Note the shift: the letter is the Cremona code for X-1.'
    col_desc['hecke_orbit_code'] = 'encoding of the tuple (g.C.N.w.i.x) into 64 bits, used in eigenvalue tables.  g + (ord(C)<<8) + (N<<12) + (k<<20) + (j<<28) + ((i-1)<<36) + ((X-1)<<52).'
    col_desc['dim'] = 'Degree of the number field over which the form is defined'
    col_desc['relative_dim'] = 'the Q(chi)-dimension of this Hecke orbit (=dim/char_degree)'
    col_desc['nf_label'] = 'LMFDB label for the corresponding number field (can be NULL)'
    col_desc['trace_hash'] = 'appropriate linear combination of the a_{p,i} between 2^12 and 2^13'
    col_desc['analytic_rank'] = 'the analytic rank of the L-function of the embedded newforms in this newform orbit'
    col_desc['analytic_rank_proved'] = 'true if the analytic rank has been proved'
    col_desc['qexp_display'] = 'latexed string for display on search page results'
    col_desc['related_objects'] = 'list of text URLs of related objects (e.g. elliptic curve isogeny class, Artin rep, ...), e.g. ["EllipticCurve/Q/11/a"]'
    col_desc['embedded_related_objects'] = 'list of lists of text URLs of related objects (e.g. Artin reps), indexed by embedding_m (so first entry is a list of friends for the first embeddeded newform)'
    col_desc['field_disc'] = 'discriminant of the coefficient field, if known'
    col_desc['field_poly_is_cyclotomic'] = "true if field_poly is a cylcotomic polynomial (the field might be Q(zeta_n) even when this flage is not set if we haven't chosen a cyclotomic polynomial to define it)"
    col_desc['field_poly'] = 'list of integers giving defining polynomial for the Hecke field (standard Sage order of coefficients)'
    col_desc['field_poly_is_real_cyclotomic'] = "true if field_poly is the minimal polynomial of zeta_n + zeta_n^-1 for some n (the field might be Q(zeta_n)^+ even when this flage is not set if we haven't chosen a cyclotomic polynomial to define it)"
    col_desc['field_poly_root_of_unity'] = 'the value of n if either field_poly_is_cylotomic of field_poly_is_real_cyclotomic is set'
    col_desc['field_disc_factorization'] = 'factorization of field discriminant stored as ordered list of pairs [p,e]'
    col_desc['trace_display'] = 'list of the first four a_{p,1} traces for display on search page results'
    col_desc['traces'] = 'full list of integer traces tr(a_{n,1}) for n from 1 to 1000 (or more)'
    col_desc['is_cuspidal'] = 'true if this is a cusp form'
    col_desc['lift_type'] = 'If not a lift, states None'
    for hecke_type in HECKE_TYPES:
        col_desc['trace_lambda_'+hecke_type] = 'List of traces of Hecke operators T_' + hecke_type + ' on the Galois orbit corresponding to the values up to 200'
    return col_desc

def create_table_smf_newforms():
    table_name = "smf_newforms"
    table_desc = "Siegel modular forms"
    generate_table(table_name, table_desc,
                   generate_column_types, generate_column_desc)
    return
