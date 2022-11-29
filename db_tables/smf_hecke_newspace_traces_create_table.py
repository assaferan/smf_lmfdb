from smf_lmfdb.db_tables.common_create_table import generate_table

def generate_column_types():
    col_type = {}
    col_type['n'] = 'integer'
    col_type['hecke_orbit_code'] = 'bigint'
    col_type['trace_an'] = 'numeric'
    return col_type

def generate_column_desc():
    col_desc = {}
    col_desc['n'] = 'index n of a_n'
    col_desc['hecke_orbit_code'] = 'encoding of the tuple (g.C.N.w.i) into 64 bits, used in eigenvalue tables.  g + (ord(C)<<8) + (N<<12) + (k<<20) + (j<<28) + ((i-1)<<36).'
    col_desc['trace_an'] = 'integer containing the nth coefficient of the trace form for the entire newspace (sum of trace forms)'
    return col_desc

def create_table_smf_hecke_newspace_traces():
    table_name = "smf_hecke_newspace_traces"
    table_desc = "Traces of Hecke operators on spaces of Siegel modular forms over number fields"
    generate_table(table_name, table_desc,
                   generate_column_types, generate_column_desc, label_col=None)
    return
    
