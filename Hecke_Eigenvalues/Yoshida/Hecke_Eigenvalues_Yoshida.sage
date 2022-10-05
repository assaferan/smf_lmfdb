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


