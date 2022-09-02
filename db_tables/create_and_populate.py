from sage.all import (is_squarefree)
from smf_lmfdb.db_tables.smf_newspaces_create_table import create_table_smf_newspaces
from smf_lmfdb.db_tables.smf_newspaces_populate import populate_smf_newspaces

def N_bound(k):
    '''
    The maximal N for which the lmfdb stores data for S_k(N) (classical modular forms)
    '''
    triv_bound = 40000 // k**2 
    small_bound = min(100000 // k**2, 10) 
    return max(triv_bound, small_bound) 

def get_box_triples(max_k = 20, max_j = 20):
    '''
    Gets all the triples (k,j,N) with k <= max_k, j <= max_j.
    Note that for N > 1 the paramodular formula only holds for k >= 3
    '''
    triple_list = [(k,j,1) for k in range(max_k+1) for j in range(max_j+1)]
    triple_list += [(k,j,N) for k in range(3,max_k+1) for j in range(max_j+1) for N
                    in range(2,N_bound(2*k+j-2)+1) if is_squarefree(N)]
    return triple_list

def create_and_populate_smf_newspaces():
    create_table_smf_newspaces()
    populate_smf_newspaces(get_box_triples())
    return
