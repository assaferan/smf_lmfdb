import os
os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db


def klingen_eis_Hecke_p(k,j):
    w = j+k 
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(200): 
        L[p]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(200):
            L[p] = (1+p^(k-2))*Tr[p-1]
    return L 

def klingen_eis_Hecke_p_square(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = (1+p^(k-2)+p^(2*k-4))*Tr[p^2-1]+p^(2*k+j-4)*(p-1)
    return L
    
def klingen_eis_Hecke_p_square_0(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = (1+(p-1)*p^(k-3)+p^(2*k-4))*Tr[p^2-1]+p^(k+j-2)*(p^(k-2)*(p-1)-p^(2*k-4)-1)
    return L

def klingen_eis_Hecke_p_square_1(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = p^(k-3)*Tr[p^2-1]+p^(k+j-2)*(1-p^(k-4)+p^(2*k-4))
    return L

def klingen_eis_Hecke_p_square_2(k,j):
    w = j+k
    results = db.mf_newforms.search({'level' : '1', 'weight': str(w)})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for res in results:
        Tr = res['traces']
        for p in prime_range(sqrt(200)):
            L[p^2] = p^(2*k+j-6)
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
        L[p^2]=0
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

def Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac(k,j,e):
    if e == 1 : 
       return Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_charac()    
    elif e == 0 : 
       return Hecke_Eigenvalues_Klingen_Eisenstein_Series(k,j)
