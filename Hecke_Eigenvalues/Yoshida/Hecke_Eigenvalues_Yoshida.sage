import os
os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db

def Y_lambda_p(k,j,p,a_p,b_p):
    return a_p+p^(k-2)*b_p

def Y_lambda_p_square(k,j,p,a_p,b_p):
    return a_p^2+p^(k-2)*a_p*b_p+p^(2*k-4)*b_p^2-p^(2*k+j-4)*(1+2*p)

def Y_lambda_p_square_0(k,j,p,a_p,b_p):
    return a_p^2+p^(2*k-4)*b_p^2+a_p*b_p*p^(k-3)*(p-1)-2*p^(2*k+j-4)*(1+p)

def Y_lambda_p_square_1(k,j,p,a_p,b_p):
    return a_p*b_p*p^(k-3)+(p^2-1)*p^(2*k+j-6)

def Y_lambda_p_square_2(k,j,p,a_p,b_p):
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
    exp = {ht : 1 if ht == 'lambda_p' else 2 for ht in hecke_types}
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
