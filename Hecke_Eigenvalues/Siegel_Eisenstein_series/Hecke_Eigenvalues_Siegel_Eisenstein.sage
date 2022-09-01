def siegel_eis_Hecke_p(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series at p
    """
    k = ZZ(k)
    p = ZZ(p)
    l=(1+p^(k-1))*(1+p^(k-2))
    return l

def siegel_eis_Hecke_p_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(200):
        evs[p] = siegel_eis_Hecke_p(p,k)
    return evs

def siegel_eis_Hecke_p_square(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series at p^2
    """
    k = ZZ(k)
    p = ZZ(p)
    l=1+p^(4*k-6)+p^(k-1)+p^(k-2)+2*p^(2*k-3)+p^(3*k -5)+p^(3*k-4)+p^(2*k-2)
    return l

def siegel_eis_Hecke_p_square_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square(p,k)
    return evs

def siegel_eis_Hecke_p_square_0(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series for the Hecke operator 
    T_{2,0}(p^2)
    """
    k = ZZ(k)
    p = ZZ(p)
    l=1+p^(4*k-6)+2*p^(2*k-3)+p^(3*k-4)+p^(k-1)+p^(2*k-2)-p^(2*k-4)-p^(3*k-6)-p^(k-3)
    return l

def siegel_eis_Hecke_p_square_0_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square_0(p,k)
    return evs

def siegel_eis_Hecke_p_square_1(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series for the Hecke operator 
    T_{2,1}(p^2)
    """
    k = ZZ(k)
    p = ZZ(p)
    l=p^(2*k-4)-p^(2*k-6)+p^(3*k-6)+p^(3*k-5)+p^(k-3)+p^(k-2)
    return l

def siegel_eis_Hecke_p_square_1_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square_1(p,k)
    return evs

def siegel_eis_Hecke_p_square_2(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series for the Hecke operator 
    T_{2,2}(p^2)
    """
    k = ZZ(k)
    p = ZZ(p)
    l=p^(2*k-6)
    return l

def siegel_eis_Hecke_p_square_2_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square_2(p,k)
    return evs


def Hecke_Eigenvalues_Siegel_Eisenstein_Series(k):
    """
    Compute                     
    """
    k=ZZ(k)
    L={}
    L['lambda_p'] = siegel_eis_Hecke_p_list(k)
    L['lambda_p_square'] = siegel_eis_Hecke_p_square_list(k)
    L['lambda_p_square_0'] = siegel_eis_Hecke_p_square_0_list(k)
    L['lambda_p_square_1'] = siegel_eis_Hecke_p_square_1_list(k)
    L['lambda_p_square_2'] = siegel_eis_Hecke_p_square_2_list(k)
    return L 

def Hecke_Eigenvalues_Siegel_Eisenstein_Series_All(k,j,e,prime_bound=200):
    empty_dic = {p : 0 for p in prime_range(prime_bound)}
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    if (j != 0) or (e != 0) or (k < 4):
       return {'lambda_' + hecke : empty_dic for hecke in hecke_types}
    return Hecke_Eigenvalues_Siegel_Eisenstein_Series(k)

def Hecke_Eigenvalues_Siegel_Eisenstein_all_forms(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = []
    if (e != 0) or (j != 0) or (k < 4):
        return forms

    hecke_types = ['p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    SE_func = {ht : eval('siegel_eis_Hecke_' + ht) for ht in hecke_types}
    exp = {ht : 1 if ht == 'p' else 2 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}
    
    orbit = {}

    orbit['related_objects'] = []
    orbit['dim'] = 1
    orbit['nf_label'] = '1.1.1.1'
    orbit['hecke_ring_index'] = 1
    orbit['field_poly_is_cyclotomic'] = False
    orbit['field_poly_root_of_unity'] = 0
    orbit['field_poly_is_real_cyclotomic'] = False
    orbit['hecke_ring_index_proved'] = True
    orbit['hecke_ring_generator_nbound'] = 1
    orbit['field_disc'] = 1
    orbit['field_disc_factorization'] = []
    orbit['field_poly'] = [0,1]
    orbit['hecke_ring_index_factorization'] = []
    orbit['relative_dim'] = 1
    orbit['is_polredabs'] = True
    orbit['is_cuspidal'] = False
    orbit['aut_rep_type'] = 'F'
        
    for ht in hecke_types:
        orbit['trace_lambda_'+ ht] = [SE_func[ht](p,k) for p in prime_range(bound[ht])]
        
    forms.append(orbit)
    return forms

def Hecke_Eigenvalues_Siegel_Eisenstein_all_evs(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    evs = []
    if (e != 0) or (j != 0) or (k < 4):
        return evs

    hecke_types = ['p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    SE_func = {ht : eval('siegel_eis_Hecke_' + ht) for ht in hecke_types}
    exp = {ht : 1 if ht == 'p' else 2 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}
    
    ev = {}
    ev['maxp'] = previous_prime(prime_bound)
    ev['maxp_square'] = previous_prime(floor(sqrt(prime_bound)))
    ev['hecke_ring_character_values'] = 'NULL'
    ev['hecke_ring_numerators'] = 'NULL'
    ev['hecke_ring_denominators'] = 'NULL'
    ev['hecke_ring_inverse_numerators'] = 'NULL'
    ev['hecke_ring_inverse_denominators'] = 'NULL'
    ev['hecke_ring_rank'] = 1
    ev['hecke_ring_power_basis'] = True
    ev['hecke_ring_cyclotomic_generator'] = 0
    ev['field_poly'] = [0,1]
    
    for ht in hecke_types:
        ev['lambda_' + ht]  = [SE_func[ht](p,k) for p in prime_range(bound[ht])]
        
    evs.append(ev)
    return evs



