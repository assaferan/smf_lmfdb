from sage.all import load

load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')

def Hecke_Eigenvalues_Siegel_Eisenstein(k,j,e,prime_bound=200):
    return Hecke_Eigenvalues_Siegel_Eisenstein_Series_All(k,j,e,prime_bound=prime_bound)

