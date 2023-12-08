# import pickle
from sage.all import (Matrix, NumberField, nth_prime, pari, PolynomialRing, prime_divisors, prime_range, QQ, primes_first_n, sqrt, factor, floor, divisors, is_squarefree, prod)
from smf_lmfdb.db_tables.common_create_table import SUBSPACE_TYPES, HECKE_TYPES
from smf_lmfdb.db_tables.nf_elt import nf_elts_to_lists
from smf_lmfdb.Dimension_formulas.paramodular.DimensionFormulas import Saito_Kurokawa_lift_dim
from smf_lmfdb.Hecke_Eigenvalues.paramodular.Hecke_Eigenvalues_Saito_Kurokawa import Hecke_Eigenvalues_PM_SK_all_evs, Hecke_Eigenvalues_PM_SK_all_forms, Hecke_Eigenvalues_PM_SK_new_num_forms

from lmfdb import db

def parse_omf5(k,j,N):
    if (k == 3) and (j == 0):
        folder = "smf_lmfdb/Hecke_Eigenvalues/paramodular/omf5_data/hecke_evs_3_0/data/"
        fname = folder + "hecke_ev_%d_%d_%d.dat" %(k,j,N)
    else:
        weight_str = str(k) + "_" + str(j)
        folder = "smf_lmfdb/Eigenforms_Weight" + weight_str + "_omf/"
        fname = folder + "hecke_ev_" + weight_str + "_" + str(N) + ".dat"
    fl = open(fname)
    r = fl.read()
    fl.close()
    # handling long numbers
    for num in range(10):
      r = r.replace(str(num) + "L", str(num))
    forms = eval(r)
    for f in forms:
        f['dim'] = len(f['field_poly'])-1
    forms = sorted(forms, key=lambda f : [f['dim']] + f['trace_lambda_p'])
    return forms

def al_str_to_num(al_str, N):
    ps = prime_divisors(N)
    return prod([ps[i] for i in range(len(al_str)) if al_str[i] == '-'])   

def Hecke_Eigenvalues_Traces_paramodular(k,j,N):
    """
    Return traces of the Hecke eigenvalues on each of the spaces of paramodular forms              
    """
    forms = parse_omf5(k,j,N)
    hecke_types = ['lambda_' + suff for suff in ['p', 'p_square']]
    num_ps = { ht : max([0] + [len(f['trace_'+ht]) for f in forms]) for ht in hecke_types}
    aut_types = {'F' : 'eis_F', 'Q' : 'eis_Q', 'P' : 'cusp_P', 'Y' : 'cusp_Y', 'G' : 'cusp_G'}
    traces = { aut_types[aut] + '_' + ht : [0 for t in range(num_ps[ht])]
               for aut in aut_types for ht in hecke_types}
    divs = [d for d in divisors(N) if is_squarefree(d)]
    new_cusp_G_dim = 0
    al_dims_G = [0 for d in divs]
    for f in forms:
        if f['aut_rep_type'] != 'G':
            continue
        f_dim = len(f['field_poly'])-1
        new_cusp_G_dim += f_dim
        div_idx = divs.index(al_str_to_num(f['atkin_lehner_string'], N))
        al_dims_G[div_idx] += f_dim
        for ht in hecke_types:
            for i in range(len(f['trace_' + ht])):
                if type(f['trace_' + ht][i]) == str:
                    traces[aut_types[f['aut_rep_type']] + '_' + ht][i] = 'NULL'
                else:
                    traces[aut_types[f['aut_rep_type']] + '_' + ht][i] += f['trace_' + ht][i]

    traces_P = Hecke_Traces_Eigenvalues_Saito_Kurokawa(k,j,N)

    for k in traces.keys():
        for i in range(len(traces[k])):
            traces[k][i] += traces_P[k][i]
    
    return traces, new_cusp_G_dim, al_dims_G

def num_forms_paramodular(k,j,N):
    forms = parse_omf5(k,j,N)
    num_G_forms = len([f for f in forms if f['aut_rep_type'] == 'G'])
    dim_G_forms = sum([len(f['field_poly'])-1 for f in forms
                       if f['aut_rep_type'] == 'G'])
    num_P_forms = Hecke_Eigenvalues_PM_SK_new_num_forms(k,j,N)
    return num_P_forms + num_G_forms, dim_G_forms

def Hecke_Eigenforms_paramodular(k,j,N):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = parse_omf5(k,j,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]

    forms = [f for f in forms if f['aut_rep_type'] == 'G']
    for orbit in forms:
        # if we have not saved the eigenvalues
        orbit['is_cuspidal'] = True
        orbit['dim'] = len(orbit['field_poly']) - 1
        orbit['is_polredabs'] = False
        # For now, all our fields are absolute. Change that in the future
        orbit['relative_dim'] = orbit['dim']
        pol = Qx([int(c) for c in orbit['field_poly']])
        orbit['field_poly_is_cyclotomic'] = pol.is_cyclotomic()
        # !! TODO - check if that actually happens to be true
        orbit['field_poly_is_real_cyclotomic'] = False
        orbit['field_poly_root_of_unity'] = pol.is_cyclotomic(certificate=True)
        if ('hecke_ring_index' in orbit) and (orbit['dim'] <= 20):
            pol = Qx(pari(pol).polredbest().polredabs())
            F = NumberField(pol, name = "a")
            coeffs = [int(c) for c in pol.coefficients(sparse=False)]
            orbit['nf_label'] = db.nf_fields.lucky({'coeffs' : coeffs}, 'label')
            orbit['field_disc'] = F.disc()
            orbit['field_disc_factorization'] = [list(fac) for fac in F.disc().factor()]
            if F.disc() < 0:
                orbit['field_disc_factorization'] = [[-1,1]] + orbit['field_disc_factorization']
            orbit['hecke_ring_index_factorization'] = [list(fac) for fac in
                                                       factor(orbit['hecke_ring_index'])]
        orbit['hecke_ring_index_proved'] = False
        if 'related_objects' not in orbit:
            orbit['related_objects'] = []
            
        # for field_name in ['hecke_ring', 'lambda_p', 'lambda_p_square']:
        for field_name in ['lambda_p', 'lambda_p_square']:
            if field_name in orbit:
                dummy = orbit.pop(field_name)

    P_forms = Hecke_Eigenvalues_PM_SK_all_forms(k,j,N)
    forms += P_forms
                
    return forms


def Hecke_Eigenvalues_paramodular(k,j,N):
    '''
    Returns a list of dictionaries.
    Each dictionary is an entry for the Hecke eigenvalues over a number field db,
    meant to be uploaded to smf_hecke_nf
    '''
    evs = parse_omf5(k,j,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]

    # we will only use those for which we have a representation of the hecke ring
    evs = [ev for ev in evs if ev['aut_rep_type'] == 'G']
    
    for ev in evs:
        F = NumberField(Qx(ev['field_poly']), name = "nu")
        nu = F.gen(0)
        
        if 'hecke_ring_numerators' in ev:
            ev['hecke_ring_cyclotomic_generator'] = 0
            ev['hecke_ring_rank'] = F.degree()
            ev['maxp'] = nth_prime(len(ev['trace_lambda_p']))
            ev['maxp_square'] = nth_prime(len(ev['trace_lambda_p_square']))
            # ev['lambda_p'] = nf_elts_to_lists(ev['lambda_p'], inv_basis)
            # ev['lambda_p_square'] = nf_elts_to_lists(ev['lambda_p_square'], inv_basis)

        for field_name in ['atkin_lehner_eigenvals',
                           'atkin_lehner_string',
                           'aut_rep_type',
                           # 'hecke_ring',
                           'hecke_ring_index',
                           'hecke_ring_generator_nbound']:
            if field_name in ev:
                dummy = ev.pop(field_name)

    evs_P = Hecke_Eigenvalues_PM_SK_all_evs(k,j,N)
    evs += evs_P
    
    return evs
