from smf_lmfdb.Hecke_Eigenvalues.Klingen_Yoshida import (
    KE_space_traces, KE_num_forms, KE_all_forms, KE_all_evs)

# The Klingen-Eisenstein series of level 1 are computed in
# Hecke_Eigenvalues/Klingen_Yoshida.py, which documents the formulas and conventions.
# Only level 1 (e = 0) has them here.

def Hecke_Traces_Eigenvalues_Klingen_Eisenstein_Series_Fast(k,j, prime_bound=200):
    """
    Traces of the Hecke eigenvalues on the Klingen-Eisenstein series of weight (k,j)
    and level 1, { ht : { p : trace } }
    """
    return KE_space_traces(k, j, prime_bound)

def Hecke_Eigenvalues_Klingen_Eisenstein_num_forms(k,j,e,prime_bound=200):
    if (e == 1):
        return 0
    return KE_num_forms(k, j)

def Hecke_Eigenvalues_Klingen_Eisenstein_all_forms(k,j,e,prime_bound=200):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    if (e == 1):
        return []
    return KE_all_forms(k, j, prime_bound)

def Hecke_Eigenvalues_Klingen_Eisenstein_all_evs(k,j,e,prime_bound=100):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    if (e == 1):
        return []
    return KE_all_evs(k, j, prime_bound)
