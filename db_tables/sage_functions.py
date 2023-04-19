#from sage.all import (factor, load, Integer, is_prime, PowerSeriesRing, previous_prime, prime_range, ZZ)
from sage.all import *
from smf_lmfdb.db_tables.common_populate import MAX_P
import os
cwd = os.getcwd()
os.chdir('smf_lmfdb/Dimension_formulas')
load('dimformula_smf_degree2_level_1.sage')
os.chdir("Principal_Congruence_Subgroup")
load('DimFormulaSMFPrincipalCogruenceSubgroup.sage')
os.chdir("../../Hecke_Eigenvalues/Siegel_Eisenstein_series")
load('Hecke_Eigenvalues_Siegel_Eisenstein.sage')
os.chdir("../Klingen_Eisenstein_series")
load('Hecke_Eigenvalues_Klingen_Eisenstein.sage')
os.chdir("../Saito_Kurokawa")
load('Hecke_Eigenvalues_Saito_Kurokawa.sage')
os.chdir("../Yoshida")
load('Hecke_Eigenvalues_Yoshida.sage')
os.chdir(cwd)

def smf_dims_degree_2_level_1(j,k,e):
    return dim_splitting_smf_degree_2_level_1(j,k,e)

def smf_dims_degree_2_level_2(k,j):
    entry = {}
    all_dims = All_List_Mult_Irrep_VV(k,j)
    # entry['degree'] = 2
    # entry['family'] = 'P'
    # entry['level'] = 2
    # entry['weight'] = [k,j]
    # entry['char_orbit'] = 0
    isotypes = ['total_dim', 'cusp_dim', 'eis_dim', 'eis_F_dim', 'eis_Q_dim',
                'cusp_P_dim', 'cusp_Y_dim', 'cusp_G_dim']
    for i in range(8):
        entry[isotypes[i]] = all_dims[i][1]
        entry['old_' + isotypes[i]] = all_dims[i][0][0]
        # for now, new is simply the ones that are not in level 1
        entry['new_' + isotypes[i]] = all_dims[i][1] - all_dims[i][0][0]
    # dimensions of the irreps according to this ordering, for later
    dims = [1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1]
    return entry

def Hecke_Eigenvalue_Traces_Klingen_Eisenstein(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Traces_Eigenvalues_Klingen_Eisenstein_Series_Fast(k,j,prime_bound=prime_bound)
    # return Hecke_Eigenvalues_Klingen_Eisenstein_Series_with_or_without_charac(k,j,e,prime_bound=prime_bound)

def Hecke_Eigenvalues_Traces_Siegel_Eisenstein(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Eigenvalues_Siegel_Eisenstein_Series_All(k,j,e,prime_bound=prime_bound)

def Hecke_Eigenvalue_Traces_Saito_Kurokawa(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Traces_Eigenvalues_Saito_Kurokawa(k,j,e,prime_bound=prime_bound)

def Hecke_Eigenvalue_Traces_Yoshida(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Traces_Eigenvalues_Yoshida(k,j,e,prime_bound=prime_bound)

def num_forms_Siegel_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Siegel_Eisenstein_num_forms(k,j,e)

def num_forms_Klingen_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Klingen_Eisenstein_num_forms(k,j,e)

def num_forms_Saito_Kurokawa(k,j,e):
    return Hecke_Eigenvalues_Saito_Kurokawa_num_forms(k,j,e)

def num_forms_Yoshida(k,j,e):
    return Hecke_Eigenvalues_Yoshida_num_forms(k,j,e)

def Hecke_Eigenforms_Siegel_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Siegel_Eisenstein_all_forms(k,j,e)

def Hecke_Eigenforms_Klingen_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Klingen_Eisenstein_all_forms(k,j,e)

def Hecke_Eigenforms_Saito_Kurokawa(k,j,e):
    return Hecke_Eigenvalues_Saito_Kurokawa_all_forms(k,j,e)

def Hecke_Eigenforms_Yoshida(k,j,e):
    return Hecke_Eigenvalues_Yoshida_all_forms(k,j,e)

def Hecke_Eigenvalues_Siegel_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Siegel_Eisenstein_all_evs(k,j,e)

def Hecke_Eigenvalues_Klingen_Eisenstein(k,j,e):
    return Hecke_Eigenvalues_Klingen_Eisenstein_all_evs(k,j,e)

def Hecke_Eigenvalues_Saito_Kurokawa(k,j,e):
    return Hecke_Eigenvalues_Saito_Kurokawa_all_evs(k,j,e)

def Hecke_Eigenvalues_Yoshida(k,j,e):
    return Hecke_Eigenvalues_Yoshida_all_evs(k,j,e)

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
            idx = 0
            while (ap[idx] == 'NULL'):
                idx += 1
            if (type(ap[idx]) == int):
                F = ZZ
            else:
                F = ap[idx].parent()
            an = F(1)
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
                    if (ap[i] == 'NULL'):
                        an = 'NULL'
                    else:
                        if (type(ap[i]) == int):
                            F = ZZ
                        else:
                            F = ap[i].parent()
                        R = PowerSeriesRing(F, 't')
                        t = R.gen()
                        P = 1 - p**(mu - 1) * t**2
                        Q = 1 - ap[i] * t + (ap[i]**2 - ap2[i] - p**(mu-1))*t**2 - p**mu * ap[i]*t**3 + p**(2*mu)*t**4
                        an = (P/Q).coefficients()[r]
            else:  # a_m*a_r := a_{mr} and we know all a_i for i<n.
                m  = fac[0][0]**fac[0][1]
                if (a[m] == 'NULL') or (a[n//m] == 'NULL'):
                    an = 'NULL'
                else:
                    an = a[m]*a[n//m]
        a += [an]
    return a
