from sage.all import (PolynomialRing, NumberField, ZZ, QQ, gcd, lcm, vector)
from lmfdb.utils import integer_squarefree_part
from smf_lmfdb.db_tables.common_populate import make_space_label, table_reload, get_hecke, common_entry_values, base_26, MAX_P, MAX_P_SQUARE
from smf_lmfdb.db_tables.smf_newforms_populate import make_orbit_code
from smf_lmfdb.db_tables.sage_functions import Hecke_Eigenvalues_Siegel_Eisenstein, Hecke_Eigenvalues_Klingen_Eisenstein, Get_All_Hecke_Eigenvalues_Up_To

from lmfdb import db

def get_nf_basis(ev):
    R = PolynomialRing(QQ, 'x')
    K = NumberField(R(ev['field_poly']), 'nu')
    nu = K.gen(0)
    dim = K.degree()

    if ev['field_poly'] == [1,0,1]:
        basis = inv_basis = [1, nu]
        return basis, inv_basis

    # Checking for the two other cyclotomic polynomials of degree 2
    if ev['hecke_ring_power_basis'] and ev['field_poly'] in [[1,-1,1],[1,1,1]]:
        basis = inv_basis = [1, nu]
        return basis, inv_basis
    
    if dim == 2:
        c, b, a = map(ZZ, ev['field_poly'])
        D = b**2 - 4*a*c
        d = integer_squarefree_part(D)
        s = (D//d).isqrt()
        if ev['hecke_ring_power_basis']:
            k, l = ZZ(0), ZZ(1)
        else:
            k, l = map(ZZ, ev['hecke_ring_numerators'][1])
            k = k / ev['hecke_ring_denominators'][1]
            l = l / ev['hecke_ring_denominators'][1]
        beta = vector((k - (b*l)/(2*a), ((s*l)/(2*a)).abs()))
        den = lcm(beta[0].denom(), beta[1].denom())
        beta *= den
        sqrt = K(d).sqrt()
        num = beta[0] + beta[1]*sqrt
        frac = num / den
        basis = inv_basis = [1, frac]
        return basis, inv_basis
    
    if ev['hecke_ring_power_basis']:
        basis = inv_basis = [nu**i for i in range(dim)]
        return basis, inv_basis
    
    coeff_data = zip(ev['hecke_ring_numerators'], ev['hecke_ring_denominators'])
    basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in coeff_data]

    inv_coeff_data = zip(ev['hecke_ring_inverse_numerators'], ev['hecke_ring_inverse_denominators'])
    inv_basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]
    
    return basis, inv_basis

def nf_lists_to_elements(coeffs, basis):
    d = len(basis)
    return [sum([coeff[i]*basis[i] for i in range(d)]) for coeff in coeffs]

def nf_elts_to_lists(elts, inv_basis):
    d = len(inv_basis)
    def to_list(elt):
        if type(elt) == int:
            return [elt] + [0 for i in range(d-1)]
        return list(elt)
    return [ list(sum([to_list(elt)[i]*inv_basis[i] for i in range(d)])) for elt in elts]

def entry_add_columns(e, ext_data):
    e['id'] = ext_data['id']
    e['char_orbit_label'] = base_26(e['char_orbit_index'])
    space_label = make_space_label(e)
    dummy = e.pop('char_orbit_label')
    hecke_orbit = ext_data['num_forms']
    e['label'] = space_label + '.' + base_26(hecke_orbit)
    e['hecke_orbit_code'] = make_orbit_code(e['degree'], e['family'], e['level'], e['weight'][0], e['weight'][1], e['char_orbit_index'], hecke_orbit)
    basis, inv_basis = get_nf_basis(e)
    e['an'] = nf_elts_to_lists(Get_All_Hecke_Eigenvalues_Up_To(e['maxp']+1, nf_lists_to_elements(e['lambda_p'], basis),
                                                               nf_lists_to_elements(e['lambda_p_square'], basis), e['weight']),
                               inv_basis)
    return e

def populate_smf_hecke_nf(triple_list):
    table = db.smf_hecke_nf
    aux_fname = "smf_hecke_nf_table.dat"
    entries = []
    for triple in triple_list:
       k,j,e = triple
       if (j % 2 == 1) or (k == 1):
           continue
       entry = common_entry_values(k,j,e)
       sub_funcs = {'eis_F' : Hecke_Eigenvalues_Siegel_Eisenstein,
                    'eis_Q' : Hecke_Eigenvalues_Klingen_Eisenstein}
       for sub in sub_funcs.keys():
           evs = sub_funcs[sub](k,j,e)
           for ev in evs:
               entry_sub = entry.copy()
               entry_sub.update(ev)
               entries.append(entry_sub)
    table_reload(table, entries, entry_add_columns, aux_fname)
    return
