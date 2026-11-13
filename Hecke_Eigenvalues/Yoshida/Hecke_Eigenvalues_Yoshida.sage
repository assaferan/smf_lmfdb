from smf_lmfdb.Hecke_Eigenvalues.Klingen_Yoshida import (
    Y_space_traces, Y_num_forms, Y_all_forms, Y_all_evs)

# The Yoshida lifts of level 2 (e = 1) are computed in
# Hecke_Eigenvalues/Klingen_Yoshida.py, which documents the formulas and conventions.

def Hecke_Traces_Eigenvalues_Yoshida(k,j,e,prime_bound=200):
    """
    Traces of the Hecke eigenvalues on the Yoshida lifts of weight (k,j), { ht : { p : trace } }
    """
    return Y_space_traces(k, j, e, prime_bound)

def Hecke_Eigenvalues_Yoshida_num_forms(k,j,e):
    return Y_num_forms(k, j, e)

def Hecke_Eigenvalues_Yoshida_all_forms(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    return Y_all_forms(k, j, e, prime_bound)

def Hecke_Eigenvalues_Yoshida_all_evs(k,j,e,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    return Y_all_evs(k, j, e, prime_bound)
