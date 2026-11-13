import os
# os.chdir("/scratch/home/fclery/lmfdb")
from lmfdb import db
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import (
    SK_query, SK_space_traces, SK_all_forms, SK_all_evs, HECKE_TYPES)

# The lifts of level 1 and of level 2 (family P) are computed as the paramodular ones
# are, in Hecke_Eigenvalues/paramodular/Hecke_Eigenvalues_Saito_Kurokawa.py, which
# documents the formulas and conventions.

def ms_label(k,N):
    return '.'.join([str(N), str(k), 'a'])

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

def Make_Zero():
    """
    Compute
    """
    L={}
    L['lambda_p'] = Make_Zero_p()
    L['lambda_p_square'] = Make_Zero_p_square()
    L['lambda_p_square_0'] = Make_Zero_p_square()
    L['lambda_p_square_1'] = Make_Zero_p_square()
    L['lambda_p_square_2'] = Make_Zero_p_square()
    return L


def SK_level_query(k, e):
    """The classical forms lifted at level 1 (e = 0) and level 2 (e = 1)."""
    if e == 0:
        return SK_query(k, 1)
    return {'level' : 2, 'weight': 2*k-2, 'atkin_lehner_string': '+'}

def Hecke_Traces_Eigenvalues_Saito_Kurokawa(k,j,e,prime_bound=200):
    """
    Traces of the Hecke eigenvalues on the space of SK lifts, { ht : { p : trace } }
    """
    if j > 0:
        return { ht : {p : 0 for p in prime_range(prime_bound)} for ht in HECKE_TYPES }
    return SK_space_traces(k, SK_level_query(k, e), prime_bound)

def Hecke_Eigenvalues_Saito_Kurokawa_num_forms(k,j,e,prime_bound=200):
    if (j > 0):
        return 0
    return db.mf_newforms.count(SK_level_query(k, e))

def Hecke_Eigenvalues_Saito_Kurokawa_all_forms(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    if (j > 0):
        return []
    return SK_all_forms(k, SK_level_query(k, e), True, prime_bound)

def Hecke_Eigenvalues_Saito_Kurokawa_all_evs(k,j,e,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    if (j > 0):
        return []
    return SK_all_evs(k, SK_level_query(k, e), prime_bound)
