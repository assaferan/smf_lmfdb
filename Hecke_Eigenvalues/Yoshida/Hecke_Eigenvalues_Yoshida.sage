import os
os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db

def Yoshida_Hecke_p(k,j):
    w1 = j+2
    w2 = j+2*k-2
    resultsplus1 = db.mf_newforms.search({'level' : '2', 'weight': str(w1), 'atkin_lehner_string': '+'})
    resultsminus1 = db.mf_newforms.search({'level' : '2', 'weight': str(w1), 'atkin_lehner_string': '-'})
    resultsplus2 = db.mf_newforms.search({'level' : '2', 'weight': str(w2), 'atkin_lehner_string': '+'})
    resultsminus2 = db.mf_newforms.search({'level' : '2', 'weight': str(w2), 'atkin_lehner_string': '-'})
    L1={}
    for p in prime_range(200):
        L1[p]=0
    for resplus1 in resultsplus1:
        Trplus1 = resplus1['traces']
        for resminus1 in resultsminus1:
            Trminus1 = resminus['traces']
            for p in prime_range(200):
                L[p] = Trplus[p-1]+ Trminus[p-1]
    return L

def Yoshida_Hecke_p_square(k,j):
    w = j+k
    resultsplus = db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': '+'})
    resultsminus = db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': '-'})
    L={}
    for p in prime_range(sqrt(200)):
        L[p^2]=0
    for resplus in resultsplus:
        Trplus = resplus['traces']
        for resminus in resultsminus:
            Trminus = resminus['traces']
            for p in prime_range(sqrt(200)):
                L[p^2] = Trplus[p-1]^2+Trminus[p-1]^2+Trplus[p-1]*Trminus[p-1]-2*(p+1)*p^j
    return L
