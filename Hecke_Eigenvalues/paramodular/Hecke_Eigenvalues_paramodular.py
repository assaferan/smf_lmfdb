# import pickle
from sage.all import (Matrix, NumberField, nth_prime, pari, PolynomialRing, prime_divisors, prime_range, QQ, primes_first_n, sqrt, factor, floor, divisors, is_squarefree, prod)
from smf_lmfdb.db_tables.common_create_table import SUBSPACE_TYPES, HECKE_TYPES
from smf_lmfdb.db_tables.nf_elt import nf_elts_to_lists

from lmfdb import db

def parse_omf5(k,j,N,hecke_ring=True):
    folder = "smf_lmfdb/Hecke_Eigenvalues/paramodular/omf5_data/hecke_evs_3_0/data/"
    fname = folder + "hecke_ev_%d_%d_%d.dat" %(k,j,N)
    # Qx = PolynomialRing(QQ, name="x")
    # x = Qx.gens()[0]
    # pickled = open(fname, "rb").read()
    # forms = pickle.loads(pickled)
    fl = open(fname)
    r = fl.read()
    fl.close()
    # handling long numbers
    for num in range(10):
      r = r.replace(str(num) + "L", str(num))
    forms = eval(r)
    return forms
    # ret = []
    # if (len(forms) > 0):
    #   f = forms[0]
    #    bad_ps = [p for p in primes_first_n(len(f['lambda_p'])) if N % p == 0]
    #    ps = primes_first_n(len(f['lambda_p']) + len(bad_ps))
    #    good_ps = [p for p in ps if p not in bad_ps]
    #    assert len(good_ps) == len(f['lambda_p'])
        
    # for f in forms:
    #    pol = Qx(f['field_poly'])
        # F = NumberField(pol, name = "a")
    #    if len(f['lambda_p']) > 0:
            # !! TODO : represent the eigenvalues in the polredabs field, currently some things break
    #        f['lambda_p'] = ['NULL' if ps[i] in bad_ps else f['lambda_p'][good_ps.index(ps[i])] for i in range(len(ps))]
    #        f['lambda_p_square'] = ['NULL' if ps[i] in bad_ps else f['lambda_p_square'][good_ps.index(ps[i])]
    #                                for i in range(len(f['lambda_p_square']))]
    #        if (hecke_ring):
    #            print("Computing hecke ring form form no.", forms.index(f), ", N = ", N)
                # !!! This can be very slow
                # hecke_ring = F.order(orbit['lambda_p'])
    #            F = f['lambda_p'][0].parent()
    #            a = F.gens()[0]
    #            index = 0
    #            p_idx = 0
    #            nbound = 0
    #            while (type(f['lambda_p'][p_idx]) == str):
    #                p_idx += 1
    #            while ((index != 1) and (p_idx < len(f['lambda_p']))):
    #                print("p_idx = ", p_idx)
    #                H = F.order([x for x in f['lambda_p'][:p_idx+1] if type(x) != str])
    #                new_index = H.index_in(F.ring_of_integers())
    #                if (new_index != index):
    #                    index = new_index
    #                    nbound = p_idx
    #                p_idx += 1
    #            f['hecke_ring'] = H
    #            f['hecke_ring_index'] = index
    #            f['hecke_ring_generator_nbound'] = nbound            
    #        ret.append(f)
    #    else:
    #        f['trace_lambda_p'] = ['NULL' if ps[i] in bad_ps else f['trace_lambda_p'][good_ps.index(ps[i])] for i in range(len(ps))]
    #        f['trace_lambda_p_square'] = ['NULL' if ps[i] in bad_ps else f['trace_lambda_p_square'][good_ps.index(ps[i])]
    #                                for i in range(len(f['trace_lambda_p_square']))]
    # return ret

def al_str_to_num(al_str, N):
    ps = prime_divisors(N)
    return prod([ps[i] for i in range(len(al_str)) if al_str[i] == '-'])   
    
def Hecke_Eigenvalues_Traces_paramodular(k,j,N, B = 100):
    """
    Return traces of the Hecke eigenvalues on each of the spaces of paramodular forms              
    """
    forms = parse_omf5(k,j,N,False)
    hecke_types = ['lambda_' + suff for suff in ['p', 'p_square']]
    num_ps = { 'lambda_p' : len(prime_range(B)), 'lambda_p_square' : len(prime_range(floor(sqrt(B))))}
    aut_types = {'F' : 'eis_F', 'Q' : 'eis_Q', 'P' : 'cusp_P', 'Y' : 'cusp_Y', 'G' : 'cusp_G'}
    traces = { aut_types[aut] + '_' + ht : [0 for t in range(num_ps[ht])]
               for aut in aut_types for ht in hecke_types}
    divs = [d for d in divisors(N) if is_squarefree(d)]
    al_dims = {'ALdims' : [0 for d in divs], 'ALdims_G' : [0 for d in divs], 'ALdims_P' : [0 for d in divs]}
    for f in forms:
        # !! TODO - handle the old forms and classify them as well
        if f['aut_rep_type'] in ['O','F','Y']:
            continue
        div_idx = al_str_to_num(f['atkin_lehner_string'], N)
        al_dims['ALdims'][div_idx] += f['dim']
        al_dims['ALdims_' + f['aut_rep_type']] += f['dim']
        for ht in hecke_types:
            for i in range(len(f['trace_' + ht])):
                if type(f['trace_' + ht][i]) == str:
                    traces[aut_types[f['aut_rep_type']] + '_' + ht][i] = 'NULL'
                else:
                    traces[aut_types[f['aut_rep_type']] + '_' + ht][i] += f['trace_' + ht][i]
    return traces.update(al_dims)   

def num_forms_paramodular(k,j,N):
    folder = "smf_lmfdb/Hecke_Eigenvalues/paramodular/omf5_data/hecke_evs_3_0/data/"
    fname = folder + "hecke_ev_%d_%d_%d.dat" %(k,j,N)
    #pickled = open(fname, "rb").read()
    #forms = pickle.loads(pickled)
    #forms = eval(open(fname).read())
    forms = parse_omf5(k,j,N)
    # return sum([len(forms[al_sign]) for al_sign in forms])
    return len(forms), sum([len(f['field_poly'])-1 for f in forms if f['aut_rep_type'] == 'G'])

def Hecke_Eigenforms_paramodular(k,j,N):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = parse_omf5(k,j,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]

    forms = [f for f in forms if f['aut_rep_type'] not in ['O','F','Y']]
    for orbit in forms:
        # if we have not saved the eigenvalues
        orbit['is_cuspidal'] = True
        orbit['dim'] = len(orbit['field_poly']) - 1
        # orbit['trace_lambda_p'] = [x.trace() if type(x) != str else 'NULL' for x in orbit['lambda_p']]
        # orbit['trace_lambda_p_square'] = [x.trace() if type(x) != str else 'NULL' for x in orbit['lambda_p_square']]
        orbit['is_polredabs'] = False
        # For now, all our fields are absolute. Change that in the future
        orbit['relative_dim'] = orbit['dim']
        pol = Qx([int(c) for c in orbit['field_poly']])
        orbit['field_poly_is_cyclotomic'] = pol.is_cyclotomic()
        # !! TODO - check if that actually happens to be true
        orbit['field_poly_is_real_cyclotomic'] = False
        orbit['field_poly_root_of_unity'] = pol.is_cyclotomic(certificate=True)
        if 'hecke_ring_index' in orbit:
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
    evs = [ev for ev in evs if ev['aut_rep_type'] not in ['O','F','Y']]
    
    for ev in evs:
        F = NumberField(Qx(ev['field_poly']), name = "nu")
        nu = F.gen(0)
        # if 'hecke_ring' in ev:
        #    basis = ev['hecke_ring'].basis()
        #    mat = Matrix([list(b) for b in basis])
        #    ev['hecke_ring_denominators'] = [row.denominator() for row in mat]
        #    ev['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
        #    ev['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
        #    ev['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat**(-1)]
            
        # inv_coeff_data = zip(ev['hecke_ring_inverse_numerators'], ev['hecke_ring_inverse_denominators'])
        # inv_basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]         
        # ev['hecke_ring_power_basis'] = False
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
    return evs
