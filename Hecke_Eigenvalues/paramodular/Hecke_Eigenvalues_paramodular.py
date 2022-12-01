from sage.all import (Matrix, NumberField, nth_prime, pari, PolynomialRing, prime_divisors, QQ)
from smf_lmfdb.db_tables.common_create_table import SUBSPACE_TYPES, HECKE_TYPES
from smf_lmfdb.db_tables.nf_elt import nf_elts_to_lists

from lmfdb import db

def parse_omf5(k,j,N,hecke_ring=True):
    folder = "smf_lmfdb/Hecke_Eigenvalues/paramodular/"
    fname = folder + "hecke_ev_%d_%d_%d.dat" %(k,j,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]
    al_signs = eval(open(fname).read())
    ret = []
    for al_sign in al_signs:
        forms = al_signs[al_sign]
        for f in forms:
            pol = Qx([int(c) for c in f['field_poly'].split()[1:]])
            F = NumberField(pol, name = "a")
            a = F.gens()[0]
            f['lambda_p'] = [F(lamda) for lamda in f['lambda_p']]
            f['lambda_p_square'] = [F(lamda) for lamda in f['lambda_p_square']]
            # !! TODO : represent the eigenvalues in the polredabs field, currently some things break
            coeffs = [int(c) for c in pol.coefficients(sparse=False)]
            f['field_poly'] = coeffs
            f['atkin_lehner_eigenvals'] = [[p, -1 if al_sign % p == 0 else 1] for p in prime_divisors(N)]
            f['atkin_lehner_string'] = ''.join(['-' if al_sign % p == 0 else '+' for p in prime_divisors(N)])
            if (hecke_ring):
                # !!! This can be very slow
                # hecke_ring = F.order(orbit['lambda_p'])
                index = 0
                p_idx = 0
                nbound = 0
                while ((index != 1) and (p_idx < len(f['lambda_p']))):
                    H = F.order(f['lambda_p'][:p_idx+1])
                    new_index = H.index_in(F.ring_of_integers())
                    if (new_index != index):
                        index = new_index
                        nbound = p_idx
                        p_idx += 1
                f['hecke_ring'] = H
                f['hecke_ring_index'] = index
                f['hecke_ring_generator_nbound'] = nbound
                        
            ret.append(f)
    return ret

def Hecke_Eigenvalues_Traces_paramodular(k,j,N):
    """
    Compute                     
    """
    forms = parse_omf5(k,j,N,False)
    hecke_types = ['lambda_' + suff for suff in ['p', 'p_square']]
    aut_types = {'F' : 'eis_F', 'Q' : 'eis_Q', 'P' : 'cusp_P', 'Y' : 'cusp_Y', 'G' : 'cusp_G'}
    traces = { aut_types[aut] + '_' + ht : [0 for t in forms[0][ht]]
               for aut in aut_types for ht in hecke_types}
    for f in forms:
        for ht in hecke_types:
            for i in range(len(f[ht])):
                traces[aut_types[f['aut_rep_type']] + '_' + ht][i] += f[ht][i].trace()
    
    return traces   

def num_forms_paramodular(k,j,N):
    folder = "smf_lmfdb/Hecke_Eigenvalues/paramodular/"
    fname = folder + "hecke_ev_%d_%d_%d.dat" %(k,j,N)
    forms = eval(open(fname).read())
    return sum([len(forms[al_sign]) for al_sign in forms])

def Hecke_Eigenforms_paramodular(k,j,N):
    '''
    Returns a list of dictionaries.
    Each dictionary is a newform orbit, meant to be uploaded to smf_newforms
    '''
    forms = parse_omf5(k,j,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]
    
    for orbit in forms:
        orbit['is_cuspidal'] = True
        orbit['dim'] = len(orbit['field_poly']) - 1
        orbit['trace_lambda_p'] = [x.trace() for x in orbit['lambda_p']]
        orbit['trace_lambda_p_square'] = [x.trace() for x in orbit['lambda_p_square']]
        orbit['is_polredabs'] = False
        # For now, all our fields are absolute. Change that in the future
        orbit['relative_dim'] = orbit['dim']
        pol = Qx([int(c) for c in orbit['field_poly']])
        orbit['field_poly_is_cyclotomic'] = pol.is_cyclotomic()
        # !! TODO - check if that actually happens to be true
        orbit['field_poly_is_real_cyclotomic'] = False
        orbit['field_poly_root_of_unity'] = pol.is_cyclotomic(certificate=True)
        pol = Qx(pari(pol).polredbest().polredabs())
        coeffs = [int(c) for c in pol.coefficients(sparse=False)]
        orbit['nf_label'] = db.nf_fields.lucky({'coeffs' : coeffs}, 'label')
        F = NumberField(pol, name = "a")
        orbit['field_disc'] = F.disc()
        orbit['field_disc_factorization'] = [list(fac) for fac in F.disc().factor()]
        if F.disc() < 0:
            orbit['field_disc_factorization'] = [[-1,1]] + orbit['field_disc_factorization']
        orbit['hecke_ring_index_factorization'] = [list(fac) for fac in
                                                   orbit['hecke_ring_index'].factor()]
        orbit['hecke_ring_index_proved'] = False
            
        for field_name in ['hecke_ring', 'lambda_p', 'lambda_p_square']:
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
    
    for ev in evs:
        F = NumberField(Qx(ev['field_poly']), name = "nu")
        basis = ev['hecke_ring'].basis()
        # basis = F.integral_basis()
        mat = Matrix([list(b) for b in basis])
        ev['hecke_ring_denominators'] = [row.denominator() for row in mat]
        ev['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
        ev['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
        ev['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat**(-1)]
        inv_coeff_data = zip(ev['hecke_ring_inverse_numerators'], ev['hecke_ring_inverse_denominators'])
        inv_basis = [sum([nums[i] * basis[i] for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]         
        ev['hecke_ring_power_basis'] = False
        ev['hecke_ring_cyclotomic_generator'] = 0
        ev['hecke_ring_rank'] = F.degree()
        ev['maxp'] = nth_prime(len(ev['lambda_p']))
        ev['maxp_square'] = nth_prime(len(ev['lambda_p_square']))
        ev['lambda_p'] = nf_elts_to_lists(ev['lambda_p'], inv_basis)
        ev['lambda_p_square'] = nf_elts_to_lists(ev['lambda_p_square'], inv_basis)

        for field_name in ['atkin_lehner_eigenvals',
                           'atkin_lehner_string',
                           'aut_rep_type',
                           'hecke_ring',
                           'hecke_ring_index',
                           'hecke_ring_generator_nbound']:
            dummy = ev.pop(field_name)
    return evs
