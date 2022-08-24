import os
os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db


def SK_Hecke_p(k):
    w = 2*k-2 
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(200): 
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(200):
            L[p] = p^(k-2)+ p^(k-1)+Tr[p-1]
    return L 

def SK_Hecke_p_square(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = Tr[p-1]^2+(p+1)*p^(k-2)*Tr[p-1]+p^(2*k-2)
    return L
    
def SK_Hecke_p_square_0(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = Tr[p-1]^2+p^(k-3)*(p^2-1)*(Tr[p-1]+p^(k-1))
    return L

def SK_Hecke_p_square_1(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = p^(k-3)*(p+1)*Tr[p-1]+p^(2*k-6)*(p^2-1)
    return L

def SK_Hecke_p_square_2(k):
    w = 2*k-2
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = p^(2*k-6)
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
   

