from sage.all import (PolynomialRing, NumberField, ZZ, QQ, gcd, lcm, vector)
from lmfdb.utils import integer_squarefree_part

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
        elif type(elt) == str:
            return elt
        return list(elt)
    return [ list(sum([to_list(elt)[i]*inv_basis[i] for i in range(d)])) if type(elt) != str else elt for elt in elts]
