from sage.all import load

load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')

def Hecke_Eigenvalue_Traces_Klingen_Eisenstein(k,j,e,prime_bound=200):
    return Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac(k,j,e,prime_bound=prime_bound)
