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

def Hecke_Eigenvalues_Siegel_Eisenstein(k,j,e,prime_bound=MAX_P+1):
    return Hecke_Eigenvalues_Siegel_Eisenstein_Series_All(k,j,e,prime_bound=prime_bound)
