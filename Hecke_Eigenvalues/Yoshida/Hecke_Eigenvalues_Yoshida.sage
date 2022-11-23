import os
# os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db
from smf_lmfdb.db_tables.nf_elt import get_nf_basis, nf_lists_to_elements, nf_elts_to_lists

def Y_lambda_p(k,j,p,a_p,b_p):
    assert k >= 2
    return a_p+p^(k-2)*b_p

def Y_lambda_p_square(k,j,p,a_p,b_p):
    assert k >= 2
    return a_p^2+p^(k-2)*a_p*b_p+p^(2*k-4)*b_p^2-p^(2*k+j-4)*(1+2*p)

def Y_lambda_p_square_0(k,j,p,a_p,b_p):
    assert k >= 2
    # right now we return a p-multiple when k == 2,
    # needs to remember to display it divided by p
    if (k == 2):
        return p*(a_p^2+p^(2*k-4)*b_p^2+a_p*b_p*p^(k-3)*(p-1)-2*p^(2*k+j-4)*(1+p))
    return a_p^2+p^(2*k-4)*b_p^2+a_p*b_p*p^(k-3)*(p-1)-2*p^(2*k+j-4)*(1+p)

def Y_lambda_p_square_1(k,j,p,a_p,b_p):
    assert k >= 2
    # right now we return a p-multiple when k == 2,
    # needs to remember to display it divided by p
    if (k == 2):
        return p*(a_p*b_p*p^(k-3)+(p^2-1)*p^(2*k+j-6))
    return a_p*b_p*p^(k-3)+(p^2-1)*p^(2*k+j-6)

def Y_lambda_p_square_2(k,j,p,a_p,b_p):
    assert k >= 2
    # right now we return a p^2-multiple when k == 2 and j == 0,
    # needs to remember to display it divided by p^2
    if (k == 2) and (j == 0):
        return p^2 * p^(2*k+j-6) 
    return p^(2*k+j-6)

def ms_label(k,N):
    return '.'.join([str(N), str(k), 'a'])

## !! TODO - can't do it anymore, need to get the data about the number field from hecke_nf table
def apply_to_nf_elt(func):
    # Here we use the fact that all the functions are linear in the number field elements,
    # so we don't need to convert back and forth, but simply scale all the coefficients
    def wrapper(*args):
        return [func(*(list(args[:-1]) + [x])) for x in args[-1]]
    return wrapper

def Yoshida1_Hecke_p(k,j):
    signs = ['+','-']
    ws = [j+2, j+2*k-2]
    result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}) for w in ws] for sign in signs]
    L = { p : 0 for p in prime_range(200)}
    for sign_idx in [0,1]:
        for g in result[sign_idx][0]:
              tr_g = g['traces']
              for f in result[1-sign_idx][1]:
                  tr_f = f['traces']
                  for p in prime_range(200):
                       L[p] += tr_f[p-1]+p^(k-2)*tr_g[p-1]
    return L

def Yoshida_Hecke_p(k,j):
    if k == 2: 
       return {key: value / 2 for key, value in Yoshida1_Hecke_p(k,j).items()}
    elif k > 2:
       return Yoshida1_Hecke_p(k,j) 

def Yoshida1_Hecke_p_square(k,j):
    signs = ['+','-']
    ws = [j+2, j+2*k-2]
    result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}) for w in ws] for sign in signs]
    L = { p^2 : 0 for p in prime_range(sqrt(200))}
    for sign_idx in [0,1]:
        for g in result[sign_idx][0]:
              tr_g = g['traces']
              for f in result[1-sign_idx][1]:
                  tr_f = f['traces']
                  for p in prime_range(sqrt(200)):
                       L[p^2] += tr_f[p-1]^2+p^(k-2)*tr_f[p-1]*tr_g[p-1]+p^(2*k-4)*tr_g[p-1]^2-p^(2*k+j-4)*(1+2*p)
    return L

def Yoshida_Hecke_p_square(k,j):
    if k == 2:
       return {key: value / 2 for key, value in Yoshida1_Hecke_p_square(k,j).items()}
    elif k > 2:
       return Yoshida1_Hecke_p_square(k,j)

def Yoshida1_Hecke_p_square_0(k,j):
    signs = ['+','-']
    ws = [j+2, j+2*k-2]
    result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}) for w in ws] for sign in signs]
    L = { p^2 : 0 for p in prime_range(sqrt(200))}
    for sign_idx in [0,1]:
        for g in result[sign_idx][0]:
              tr_g = g['traces']
              for f in result[1-sign_idx][1]:
                  tr_f = f['traces']
                  for p in prime_range(sqrt(200)):
                       L[p^2] += tr_f[p-1]^2+p^(2*k-4)*tr_g[p-1]^2+tr_f[p-1]*tr_g[p-1]*p^(k-3)*(p-1)-2*p^(2*k+j-4)*(1+p)
    return L

def Yoshida_Hecke_p_square_0(k,j):
    if k == 2:
       return {key: value / 2 for key, value in Yoshida1_Hecke_p_square_0(k,j).items()}
    elif k > 2:
       return Yoshida1_Hecke_p_square_0(k,j)


def Yoshida1_Hecke_p_square_1(k,j):
    signs = ['+','-']
    ws = [j+2, j+2*k-2]
    result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}) for w in ws] for sign in signs]
    L = { p^2 : 0 for p in prime_range(sqrt(200))}
    for sign_idx in [0,1]:
        for g in result[sign_idx][0]:
              tr_g = g['traces']
              for f in result[1-sign_idx][1]:
                  tr_f = f['traces']
                  for p in prime_range(sqrt(200)):
                       L[p^2] += tr_f[p-1]*tr_g[p-1]*p^(k-3)+(p^2-1)*p^(2*k+j-6)
    return L

def Yoshida_Hecke_p_square_1(k,j):
    if k == 2:
       return {key: value / 2 for key, value in Yoshida1_Hecke_p_square_1(k,j).items()}
    elif k > 2:
       return Yoshida1_Hecke_p_square_1(k,j)



def Yoshida1_Hecke_p_square_2(k,j):
    signs = ['+','-']
    ws = [j+2, j+2*k-2]
    result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}) for w in ws] for sign in signs]
    L = { p^2 : 0 for p in prime_range(sqrt(200))}
    for sign_idx in [0,1]:
        for g in result[sign_idx][0]:
              tr_g = g['traces']
              for f in result[1-sign_idx][1]:
                  tr_f = f['traces']
                  for p in prime_range(sqrt(200)):
                       L[p^2] += p^(2*k+j-6)
    return L

def Yoshida_Hecke_p_square_2(k,j):
    if k == 2:
       return {key: value / 2 for key, value in Yoshida1_Hecke_p_square_2(k,j).items()}
    elif k > 2:
       return Yoshida1_Hecke_p_square_2(k,j)

def Hecke_Eigenvalues_Yoshida(k,j):
    """
    Compute
    """
    L={}
    L['lambda_p'] = Yoshida_Hecke_p(k,j)
    L['lambda_p_square'] = Yoshida_Hecke_p_square(k,j)
    L['lambda_p_square_0'] = Yoshida_Hecke_p_square_0(k,j)
    L['lambda_p_square_1'] = Yoshida_Hecke_p_square_1(k,j)
    L['lambda_p_square_2'] = Yoshida_Hecke_p_square_2(k,j)
    return L


def Make_Zero_p():
    L={}
    for p in prime_range(200):
        L[p]=0
    return L

def Make_Zero_p_square():
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    return L

def Make_Unknown_p():
    L={}
    for p in prime_range(200):
        L[p]='?'
    return L

def Make_Unknown_p_square():
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]='?'
    return L


def Make_Zero():
    L={}
    L['lambda_p'] = Make_Zero_p()
    L['lambda_p_square'] = Make_Zero_p_square()
    L['lambda_p_square_0'] = Make_Zero_p_square()
    L['lambda_p_square_1'] = Make_Zero_p_square()
    L['lambda_p_square_2'] = Make_Zero_p_square()
    return L

def Make_Unknown():
    L={}
    L['lambda_p'] = Make_Unknown_p()
    L['lambda_p_square'] = Make_Unknown_p_square()
    L['lambda_p_square_0'] = Make_Unknown_p_square()
    L['lambda_p_square_1'] = Make_Unknown_p_square()
    L['lambda_p_square_2'] = Make_Unknown_p_square()
    return L

def Hecke_Eigenvalues_Yoshida_All(k,j,e):
    if e == 0 or (j % 2) == 1 : 
       return Make_Zero()
    elif e ==1 :
       return Make_Zero()   
    elif e == 1 and k > 2 :
       return Hecke_Eigenvalues_Yoshida(k,j)       
    elif k == 2 and j > 33 :
       return Make_Unknown()
    elif k == 2 and  j < 33 :
       return Hecke_Eigenvalues_Yoshida(k,j)

def Hecke_Traces_Eigenvalues_Yoshida(k,j,e,prime_bound=200):
    """
    Compute                     
    """
    ws = [j+2, j+2*k-2]
    signs = ['+','-']
    fields = ['traces']
    if ws[1] != ws[0]:
        result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}, fields) for w in ws] for sign in signs] 
        pairs = [[[f,g] for f in result[1-sign_idx][1] for g in result[sign_idx][0]] for sign_idx in [0,1]]
        pairs = reduce(lambda x,y : x+y, pairs)
    else:
        result = [db.mf_newforms.search({'level' : '2', 'weight': str(ws[0]), 'atkin_lehner_string': sign}, fields) for sign in signs]
        pairs = [[[f,g] for f in result[1-sign_idx] for g in result[sign_idx]] for sign_idx in [0,1]]
        pairs = reduce(lambda x,y : x+y, pairs)
                
    hecke_types = ['lambda_p' + suffix for suffix in [''] +
                   ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    Y_func = {ht : eval('Y_' + ht) for ht in hecke_types}
    exp = {ht : 1 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}
    ranges = {ht : prime_range(bound[ht]) for ht in hecke_types}  
    L = { ht : {p : 0 for p in ranges[ht]}  for ht in hecke_types }
    if (pairs and (e > 0) and is_even(j) and (k >= 2)):
        for ht in hecke_types:
            for p in ranges[ht]:
                Tr_a = 0
                Tr_b = 0
                for pair in pairs:
                    Tr_a += pair[0]['traces'][p^exp[ht]-1]
                    Tr_b += pair[1]['traces'][p^exp[ht]-1]
                L[ht][p] = Y_func[ht](k,j,p,Tr_a,Tr_b) 
    return L   
   
def Hecke_Eigenvalues_Yoshida_all_forms(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = []
    if (e == 0) or is_odd(j) or (k < 2):
        return forms
    ws = [j+2, j+2*k-2]
    signs = ['+','-']
    fields = ['label', 'dim', 'nf_label', 'hecke_ring_index', 'field_poly_is_cyclotomic', 'field_poly_root_of_unity', 'field_poly_is_real_cyclotomic', 'hecke_ring_index_proved', 'hecke_ring_generator_nbound', 'field_disc', 'field_disc_factorization', 'field_poly', 'hecke_ring_index_factorization', 'relative_dim', 'traces', 'is_polredabs', 'atkin_lehner_string']
    if ws[1] != ws[0]:
        result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}, fields) for w in ws] for sign in signs] 
        pairs = [[[f,g] for f in result[1-sign_idx][1] for g in result[sign_idx][0]] for sign_idx in [0,1]]
        pairs = reduce(lambda x,y : x+y, pairs)
    else:
        result = [db.mf_newforms.search({'level' : '2', 'weight': str(ws[0]), 'atkin_lehner_string': sign}, fields) for sign in signs]
        pairs = [[[f,g] for f in result[1-sign_idx] for g in result[sign_idx]] for sign_idx in [0,1]]
        pairs = reduce(lambda x,y : x+y, pairs)
    
    hecke_types = ['lambda_p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    Y_func = {ht : eval('Y_' + ht) for ht in hecke_types}
    exp = {ht : 1 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}

    for pair in pairs:
        f,g = pair
        orbit = f
        orbit['related_objects'] = [ 'ModularForm/GL2/Q/holomorphic/' + '/'.join(orb['label'].split('.')) for orb in [f,g]]
        orbit['is_cuspidal'] = True
        orbit['aut_rep_type'] = 'Y'
        orbit['dim'] = f['dim']*g['dim']
        orbit['atkin_lehner_string'] = '-'
        for ht in hecke_types:
            orbit['trace_' + ht] = [Y_func[ht](k,j,p,f['traces'][p-1],g['traces'][p-1]) for p in prime_range(bound[ht])]
        # We don't want to accidentally use the same label or traces for the Yoshida form
        for field_name in ['label', 'traces']:
            dummy = orbit.pop(field_name)
        forms.append(orbit)
    return forms

# For now we restrict to 100 since column 'an' of mf_hecke_nf only stores up to 100
# !! TODO - the ap 'column' stores further (up to 997). W can use it to get to the a_{p^2}
# However, for that we will need to perform actual field arithmetic so we save it for later

def Hecke_Eigenvalues_Yoshida_all_evs(k,j,e,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    evs = []
    if (e == 0) or is_odd(j) or (k < 2):
        return evs
    ws = [j+2, j+2*k-2]
    signs = ['+','-']
    fields = ['hecke_orbit_code']
    if ws[1] != ws[0]:
        result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}, fields) for w in ws] for sign in signs] 
        pairs = [[[f,g] for f in result[1-sign_idx][1] for g in result[sign_idx][0]] for sign_idx in [0,1]]
        pairs = reduce(lambda x,y : x+y, pairs)
    else:
        result = [db.mf_newforms.search({'level' : '2', 'weight': str(ws[0]), 'atkin_lehner_string': sign}, fields) for sign in signs]
        pairs = [[[f,g] for f in result[1-sign_idx] for g in result[sign_idx]] for sign_idx in [0,1]]
        pairs = reduce(lambda x,y : x+y, pairs)

    hecke_types = ['lambda_p' + suffix for suffix in [''] + ['_square' + sfx for sfx in [''] + ['_' + str(i) for i in range(3)]]]
    Y_func = {ht : eval('Y_' + ht) for ht in hecke_types}
    exp = {ht : 1 for ht in hecke_types}
    bound = {ht : previous_prime(floor(prime_bound^(1/exp[ht])))+1 for ht in hecke_types}
    for orbit_pair in pairs:
        ev_pair = [db.mf_hecke_nf.lucky({'hecke_orbit_code' : orbit['hecke_orbit_code']}) for orbit in orbit_pair]
        if (ev_pair[0] and ev_pair[1]):
            basis0, inv_basis0 = get_nf_basis(ev_pair[0])
            basis1, inv_basis1 = get_nf_basis(ev_pair[1])
            F0 = basis0[-1].parent()
            F1 = basis1[-1].parent()
            F = F1.composite_fields(F0)[0]
            nu = F.gen(0)
            basis = F.integral_basis()
            mat = Matrix([list(b) for b in basis])
            ev = ev_pair[0]
            ev['hecke_ring_denominators'] = [row.denominator() for row in mat]
            ev['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
            ev['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat^(-1)]
            ev['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat^(-1)]
            inv_coeff_data = zip(ev['hecke_ring_inverse_numerators'], ev['hecke_ring_inverse_denominators'])
            inv_basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]         
            ev['field_poly'] = list(F.defining_polynomial())
            ev['hecke_ring_power_basis'] = False
            ev['hecke_ring_cyclotomic_generator'] = 0
            ev['hecke_ring_rank'] = F.degree()
            ev['maxp'] = bound['lambda_p'] -1 
            ev['maxp_square'] = bound['lambda_p_square'] - 1
            for ht in hecke_types:
                aps_lists = [ev_pair[0]['an'][p^exp[ht]-1] for p in prime_range(bound[ht])]
                bps_lists = [ev_pair[1]['an'][p^exp[ht]-1] for p in prime_range(bound[ht])]
                aps = nf_lists_to_elements(aps_lists, basis0)
                bps = nf_lists_to_elements(bps_lists, basis1)
                aps = [F(ap) for ap in aps]
                bps = [F(bp) for bp in bps]
                ps = prime_range(bound[ht])
                ev[ht] = nf_elts_to_lists([F(Y_func[ht](k,j,ps[i],aps[i],bps[i])) for i in range(len(ps))], inv_basis)
            # We don't want to accidentally use the same fields for the Yoshida lift
            for field_name in ['label', 'hecke_orbit_code', 'an', 'ap', 'weight']:
                dummy = ev.pop(field_name)
            evs.append((ev))
    return evs
