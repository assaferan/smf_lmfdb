from lmfdb import db

def KE_lambda_p(k,j,p,a_p):
    return (1+p^(k-2))*a_p

def KE_lambda_p_square(k,j,p,a_p2):
    return (1+p^(k-2)+p^(2*k-4))*a_p2+p^(2*k+j-4)*(p-1)

def KE_lambda_p_square_0(k,j,p,a_p2):
    return (1+(p-1)*p^(k-3)+p^(2*k-4))*a_p2+p^(k+j-2)*(p^(k-2)*(p-1)-p^(2*k-4)-1)

def KE_lambda_p_square_1(k,j,p,a_p2):
    return p^(k-3)*a_p2+p^(k+j-2)*(1-p^(k-4)+p^(2*k-4))

def KE_lambda_p_square_2(k,j,p,a_p):
    return p^(2*k+j-6)

def apply_to_nf_elt(func):
    # Here we use the fact that all the functions are linear in the number field elements,
    # so we don't need to convert back and forth, but simply scale all the coefficients
    def wrapper(*args):
        return [func(*(list(args[:-1]) + [x])) for x in args[-1]]
    return wrapper

def klingen_eis_Hecke_p(k,j):
    w = j+k 
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['traces'])
    L={}
    for p in prime_range(200): 
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(200):
            L[p] += (1+p^(k-2))*Tr[p-1]
    return L 

def klingen_eis_Hecke_p_square(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['traces'])
    L={}
    for p in prime_range(sqrt(200)):
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p] += (1+p^(k-2)+p^(2*k-4))*Tr[p^2-1]+p^(2*k+j-4)*(p-1)
    return L
    
def klingen_eis_Hecke_p_square_0(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['traces'])
    L={}
    for p in prime_range(sqrt(200)):
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p] += (1+(p-1)*p^(k-3)+p^(2*k-4))*Tr[p^2-1]+p^(k+j-2)*(p^(k-2)*(p-1)-p^(2*k-4)-1)
    return L

def klingen_eis_Hecke_p_square_1(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['traces'])
    L={}
    for p in prime_range(sqrt(200)):
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p] += p^(k-3)*Tr[p^2-1]+p^(k+j-2)*(1-p^(k-4)+p^(2*k-4))
    return L

def klingen_eis_Hecke_p_square_2(k,j):
    # do we even need to read from the database here?
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['traces'])
    L={}
    for p in prime_range(sqrt(200)):
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p] += p^(2*k+j-6)
    return L


def Hecke_Eigenvalues_Klingen_Eisenstein_Series(k,j):
    """
    Compute                     
    """
    L={}
    L['lambda_p'] = klingen_eis_Hecke_p(k,j)
    L['lambda_p_square'] = klingen_eis_Hecke_p_square(k,j)
    L['lambda_p_square_0'] = klingen_eis_Hecke_p_square_0(k,j)
    L['lambda_p_square_1'] = klingen_eis_Hecke_p_square_1(k,j)
    L['lambda_p_square_2'] = klingen_eis_Hecke_p_square_2(k,j)
    return L   
   
def Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p():
    L={}
    for p in prime_range(200):
        L[p]=0
    return L

def Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p_square():
    L={}
    for p in prime_range(sqrt(200)):
        L[p]=0
    return L

def Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac():
    """
    Compute
    """
    L={}
    L['lambda_p'] = Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p()
    L['lambda_p_square'] = Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p_square()
    L['lambda_p_square_0'] = Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p_square()
    L['lambda_p_square_1'] = Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p_square()
    L['lambda_p_square_2'] = Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac_p_square()
    return L

def Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac(k,j,e,prime_bound=200):
    empty_dic = {p : 0 for p in prime_range(prime_bound)}
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    if (k < 4):
       return {'lambda_' + hecke : empty_dic for hecke in hecke_types}
    if e == 1 : 
       return Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac()    
    elif e == 0 : 
       return Hecke_Eigenvalues_Klingen_Eisenstein_Series(k,j)

def Hecke_Eigenvalues_Klingen_Eisenstein_all_forms(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = []
    if (e == 1) or (k < 4):
        return forms
    w = j + k
    orbits = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['label', 'dim', 'nf_label', 'hecke_ring_index', 'field_poly_is_cyclotomic', 'field_poly_root_of_unity', 'field_poly_is_real_cyclotomic', 'hecke_ring_index_proved', 'hecke_ring_generator_nbound', 'field_disc', 'field_disc_factorization', 'field_poly', 'hecke_ring_index_factorization', 'relative_dim', 'traces', 'is_polredabs'])
    hecke_types = ['lambda_p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    KE_func = {ht : apply_to_nf_elt(eval('KE_' + ht)) for ht in hecke_types}
    exp = {ht : 1 if ht == 'lambda_p' else 2 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}
    for orbit in orbits:
        orbit['related_objects'] = [ 'ModularForm/GL2/Q/holomorphic/' + '/'.join(orbit['label'].split('.'))]
        orbit['is_cuspidal'] = False
        orbit['aut_rep_type'] = 'Q'
        for ht in hecke_types:
            orbit['trace_' + ht] = [KE_func[ht](k,j,p,orbit['traces'][p^exp[ht]-1]) for p in prime_range(bound[ht])]
        # We don't want to accidentally use the same label or traces for the Klingen-Eisenstein form
        for field_name in ['label', 'traces']:
            dummy = orbit.pop(field_name)
        forms.append(orbit)
    return forms

# For now we restrict to 100 since column 'an' of mf_hecke_nf only stores up to 100
# !! TODO - the ap 'column' stores further (up to 997). W can use it to get to the a_{p^2}
# However, for that we will need to perform actual field arithmetic so we save it for later

def Hecke_Eigenvalues_Klingen_Eisenstein_all_evs(k,j,e,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    evs = []
    if (e == 1) or (k < 4):
        return evs
    w = j + k
    orbits = db.mf_newforms.search({'level' : '1', 'weight': str(w)}, ['hecke_orbit_code'])
    hecke_types = ['lambda_p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    KE_func = {ht : apply_to_nf_elt(eval('KE_' + ht)) for ht in hecke_types}
    exp = {ht : 1 if ht == 'lambda_p' else 2 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}
    for orbit in orbits:
        ev = db.mf_hecke_nf.lucky({'hecke_orbit_code' : orbit['hecke_orbit_code']})
        if ev:
            ev['maxp'] = bound['lambda_p'] -1 
            ev['maxp_square'] = bound['lambda_p_square'] - 1
            for ht in hecke_types:
                ev[ht] = [KE_func[ht](k,j,p,ev['an'][p^exp[ht]-1]) for p in prime_range(bound[ht])]
            # We don't want to accidentally use the same fields for the Klingen-Eisenstein form
            for field_name in ['label', 'hecke_orbit_code', 'an', 'ap', 'weight']:
                dummy = ev.pop(field_name)
            evs.append((ev))
    return evs
