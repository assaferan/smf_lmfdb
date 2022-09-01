from sage.all import *
from smf_lmfdb.db_tables.common_populate import MAX_P
import os
cwd = os.getcwd()
os.chdir('smf_lmfdb/Dimension_formulas')
load('dimformula_smf_degree2_level_1.sage')
os.chdir("../Hecke_Eigenvalues/Siegel_Eisenstein_series")
load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')
os.chdir("../Klingen_Eisenstein_series")
load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')
os.chdir(cwd)

def smf_dims_degree_2_level_1(j,k,e):
    return dim_splitting_smf_degree_2_level_1(j,k,e)

def Hecke_Eigenvalue_Traces_Klingen_Eisenstein(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac(k,j,e,prime_bound=prime_bound)

def Hecke_Eigenvalues_Traces_Siegel_Eisenstein(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Eigenvalues_Siegel_Eisenstein_Series_All(k,j,e,prime_bound=prime_bound)

def Hecke_Eigenforms_Siegel_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Siegel_Eisenstein_all_forms(k,j,e)

def Hecke_Eigenforms_Klingen_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Klingen_Eisenstein_all_forms(k,j,e)

def Hecke_Eigenvalues_Siegel_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Siegel_Eisenstein_all_evs(k,j,e)

def Hecke_Eigenvalues_Klingen_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Klingen_Eisenstein_all_evs(k,j,e)


# TODO : this is highly inefficient, can sieve it through, and can compue the exact recurrence relation for powers of p
def Get_All_Hecke_Eigenvalues_Up_To(prec, ap, ap2, wt):
    k,j = wt
    if (k == 0) and (j == 0):
        return [0 for n in range(prec)]
    a = []
    primes = prime_range(prec+1)
    prime_idx = {primes[i] : i for i in range(len(primes))}
    for n in range(1,prec+1):
        if n == 1:
            an = 1
        elif is_prime(n):
            an = ap[prime_idx[n]]
        else:
            fac = factor(n) 
            if len(fac) == 1:
                # computing a_{p^r} 
                p  = fac[0][0]
                r  = fac[0][1]
                i = prime_idx[p]
                if r == 2:
                    an = ap2[i]
                else:
                    mu = j + 2*k - 3
                    ZZ = ap[i].parent()
                    R = PowerSeriesRing(ZZ, 't')
                    t = R.gen()
                    P = 1 - p**(mu - 1) * t**2
                    Q = 1 - ap[i] * t + (ap[i]**2 - ap2[i] - p**(mu-1))*t**2 - p**mu * ap[i]*t**3 + p**(2*mu)*t**4
                    an = (P/Q).coefficients()[r]
            else:  # a_m*a_r := a_{mr} and we know all a_i for i<n.
                m  = fac[0][0]**fac[0][1]
                an = a[m]*a[n//m]
        a += [an]
    return a
