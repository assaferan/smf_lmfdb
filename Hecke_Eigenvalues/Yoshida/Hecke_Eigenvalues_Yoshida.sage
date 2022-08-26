import os
os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db

def Yoshida_Hecke_p(k,j):
    signs = ['+','-']
    ws = [j+2, j+2*k-2]
    result = [[db.mf_newforms.search({'level' : '2', 'weight': str(w), 'atkin_lehner_string': sign}) for w in ws] for sign in signs]
    L = { p : 0 for p in prime_range(200)}
    for sign_idx in [0,1]:
        for g in result[signs[sign_idx]][0]:
              tr_g = g['traces']
              for f in result[signs[1-sign_idx]][1]:
                  tr_f = f['traces']
                  for p in prime_range(200):
                       L[p] += tr_f+p^(k-2)*tr_g
    return L