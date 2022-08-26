import os
os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db

def Yoshida_Hecke_p(k,j):
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


def Yoshida_Hecke_p_square(k,j):
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
