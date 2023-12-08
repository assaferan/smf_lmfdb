from sage.arith.misc import previous_prime
from sage.functions.other import floor
from sage.rings.fast_arith import prime_range

from lmfdb import db
from smf_lmfdb.db_tables.common_populate import MAX_P
from smf_lmfdb.db_tables.nf_elt import apply_to_nf_elt

def ms_label(k,N):
    return '.'.join([str(N), str(k), "a"])

def SK_lambda_p(k,p,a_p):
    return p**(k-2)+ p**(k-1)+a_p

def SK_lambda_p_square(k,p,a_p):
    return a_p**2+(p+1)*p**(k-2)*a_p+p**(2*k-2)

def SK_lambda_p_square_0(k,p,a_p):
    return a_p**2+p**(k-3)*(p^2-1)*(a_p+p**(k-1))

def SK_lambda_p_square_1(k,p,a_p):
    return p**(k-3)*(p+1)*a_p+p**(2*k-6)*(p**2-1)

def SK_lambda_p_square_2(k,p,a_p):
    return p**(2*k-6)

def SK_Hecke_p(k):
    w = 2*k-2 
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(MAX_P+1): 
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(MAX_P+1):
            L[p] = p**(k-2)+ p**(k-1)+Tr[p-1]
    return L 

def SK_Hecke_p_square(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(MAX_P+1)):
        L[p**2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(MAX_P+1)):
            L[p**2] = Tr[p-1]**2+(p+1)*p**(k-2)*Tr[p-1]+p**(2*k-2)
    return L
    
def SK_Hecke_p_square_0(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(MAX_P+1)):
        L[p**2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(MAX_P+1)):
            L[p**2] = Tr[p-1]**2+p**(k-3)*(p**2-1)*(Tr[p-1]+p**(k-1))
    return L

def SK_Hecke_p_square_1(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(MAX_P+1)):
        L[p**2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(MAX_P+1)):
            L[p**2] = p**(k-3)*(p+1)*Tr[p-1]+p**(2*k-6)*(p**2-1)
    return L

def SK_Hecke_p_square_2(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(MAX_P+1)):
        L[p**2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(MAX_P+1)):
            L[p**2] = p**(2*k-6)
    return L


def Hecke_Eigenvalues_SK(k):
    """
    Compute                     
    """
    L={}
    L['lambda_p'] = SK_Hecke_p(k)
    L['lambda_p_square'] = SK_Hecke_p_square(k)
    L['lambda_p_square_0'] = SK_Hecke_p_square_0(k)
    L['lambda_p_square_1'] = SK_Hecke_p_square_1(k)
    L['lambda_p_square_2'] = SK_Hecke_p_square_2(k)
    return L   

def Make_Zero_p():
    L={}
    for p in prime_range(MAX_P+1):
        L[p]=0
    return L

def Make_Zero_p_square():
    L={}
    for p in prime_range(sqrt(MAX_P+1)):
        L[p**2]=0
    return L

def Make_Zero():
    """
    Compute
    """
    L={}
    L['lambda_p'] = Make_Zero_p()
    L['lambda_p_square'] = Make_Zero_p_square()
    L['lambda_p_square_0'] = Make_Zero_p_square()
    L['lambda_p_square_1'] = Make_Zero_p_square()
    L['lambda_p_square_2'] = Make_Zero_p_square()
    return L


def Hecke_Eigenvalues_SK_All(k,j):
    if j > 0 : 
       return Make_Zero()   
    elif j == 0: 
       return Hecke_Eigenvalues_SK(k)

# Unfortunately mf_newspaces does not have the data of traces on +/- spaces,
# so we resort to actually getting the forms
def Hecke_Traces_Eigenvalues_Saito_Kurokawa(k,j,N,prime_bound=200):
    """
    Compute                     
    """
    w = 2*k-2
    query = {'level' : N, 'weight': str(w), 'fricke_eigenval' : -1}
    res = db.mf_newforms.search(query, ['traces'])
    res = {'traces' : sum([vector(f['traces']) for f in res])}
    SK_func = {ht : eval('SK_' + ht) for ht in hecke_types}
    exp = {ht : 1 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound**(1/exp[ht])))+1 for ht in hecke_types}
    ranges = {ht : prime_range(bound[ht]) for ht in hecke_types}  
    L = { ht : {p : 0 for p in ranges[ht]}  for ht in hecke_types }
    if (res and (j == 0)):
        Tr = res['traces']
        for ht in hecke_types:
            for p in ranges[ht]:
                L[ht][p] = SK_func[ht](k,p,Tr[p**exp[ht]-1]) 
    return L   

def Hecke_Eigenvalues_PM_SK_new_num_forms(k,j,N):
    if (j > 0):
        return 0
    w = 2*k-2
    query = {'level' : N, 'weight': str(w), 'fricke_eigenval' : -1}
    num_orbits = db.mf_newforms.count(query)
    return num_orbits

# We probably don't need this one, but just in case
def Hecke_Eigenvalues_PM_SK_num_forms(k,j,N):
    return sum([Hecke_Eigenvalues_PM_SK_new_num_forms(k,j,M) for M in divisors(N)])

def Hecke_Eigenvalues_PM_SK_all_forms(k,j,N,prime_bound=MAX_P+1):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = []
    if (j > 0):
        return forms
    w = 2*k-2
    query = {'level' : N, 'weight': str(w), 'fricke_eigenval' : -1}
        
    orbits = db.mf_newforms.search(query, ['label', 'dim', 'nf_label', 'hecke_ring_index', 'field_poly_is_cyclotomic', 'field_poly_root_of_unity', 'field_poly_is_real_cyclotomic', 'hecke_ring_index_proved', 'hecke_ring_generator_nbound', 'field_disc', 'field_disc_factorization', 'field_poly', 'hecke_ring_index_factorization', 'relative_dim', 'traces', 'is_polredabs'])
    hecke_types = ['lambda_p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    SK_func = {ht : eval('SK_' + ht) for ht in hecke_types}
    exp = {ht : 1 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound**(1/exp[ht])))+1 for ht in hecke_types}
    for orbit in orbits:
        orbit['related_objects'] = [ 'ModularForm/GL2/Q/holomorphic/' + '/'.join(orbit['label'].split('.'))]
        orbit['is_cuspidal'] = True
        orbit['aut_rep_type'] = 'P'
        for ht in hecke_types:
            orbit['trace_' + ht] = [SK_func[ht](k,p,orbit['traces'][p**exp[ht]-1]) for p in prime_range(bound[ht])]
        # We don't want to accidentally use the same label or traces for the Saito-Kurokawa lift
        for field_name in ['label', 'traces']:
            dummy = orbit.pop(field_name)
        forms.append(orbit)
    return forms

# For now we restrict to 100 since column 'an' of mf_hecke_nf only stores up to 100
# !! TODO - the ap 'column' stores further (up to 997). W can use it to get to the a_{p^2}
# However, for that we will need to perform actual field arithmetic so we save it for later

def Hecke_Eigenvalues_PM_SK_all_evs(k,j,N,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    evs = []
    if (j > 0):
        return evs
    w = 2*k-2
    query = {'level' : N, 'weight': str(w), 'fricke_eigenval' : -1}
    orbits = db.mf_newforms.search(query, ['hecke_orbit_code'])
    hecke_types = ['lambda_p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    SK_func = {ht : apply_to_nf_elt(eval('SK_' + ht)) for ht in hecke_types}
    exp = {ht : 1 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound**(1/exp[ht])))+1 for ht in hecke_types}
    for orbit in orbits:
        ev = db.mf_hecke_nf.lucky({'hecke_orbit_code' : orbit['hecke_orbit_code']})
        if ev:
            ev['maxp'] = bound['lambda_p'] -1 
            ev['maxp_square'] = bound['lambda_p_square'] - 1
            for ht in hecke_types:
                ev[ht] = [SK_func[ht](k,p,ev['an'][p**exp[ht]-1],ev) for p in prime_range(bound[ht])]
            # We don't want to accidentally use the same fields for the Saito-Kurokawa lift
            for field_name in ['label', 'hecke_orbit_code', 'an', 'ap', 'weight']:
                dummy = ev.pop(field_name)
            evs.append((ev))
    return evs
