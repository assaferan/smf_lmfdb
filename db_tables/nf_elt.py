from sage.all import (Integer, Matrix, PolynomialRing, NumberField, ZZ, QQ, gcd, lcm, vector)
from lmfdb.utils import integer_squarefree_part

def get_nf_basis(ev, two=False):
    """
    Return (basis, inv_basis) for the Hecke ring of ev: basis[i] is the i-th
    generator (hecke_ring_numerators[i] / hecke_ring_denominators[i], or the
    power basis), and inv_basis[i] is the i-th row of the inverse change of
    basis, as nf_elts_to_lists expects.

    Coordinates stored in the database are with respect to this stored basis,
    so that is the default.  With two=True a quadratic field instead gets the
    basis 1, beta with beta = (a + b sqrt(d))/den, b > 0; this is only allowed
    when it is the stored basis (B_0 = 1 and beta = B_1), since otherwise stored
    coordinates would be misread (the cause of the wrong degree-2 traces in
    omf5_data, see assaferan/omf5_data#2).
    """
    R = PolynomialRing(QQ, 'x')
    K = NumberField(R(ev['field_poly']), 'nu')
    nu = K.gen(0)
    dim = K.degree()

    # The cyclotomic polynomials of degree 2, in the power basis
    if ev['hecke_ring_power_basis'] and ev['field_poly'] in [[1,0,1],[1,-1,1],[1,1,1]]:
        basis = inv_basis = [1, nu]
        return basis, inv_basis

    if (dim == 2) and (two):
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
        basis = [K(1), frac]
        if not ev['hecke_ring_power_basis']:
            stored, _ = get_nf_basis(ev, False)
            if basis != stored:
                raise ValueError("two=True: the basis 1, %s is not the stored Hecke ring basis %s" % (frac, stored))
        elif basis != [K(1), nu]:
            raise ValueError("two=True: the basis 1, %s is not the stored power basis" % frac)
        # rows of the inverse of the matrix whose rows are the basis in the power basis
        Minv = Matrix(QQ, [list(b) for b in basis]).inverse()
        inv_basis = [K(list(row)) for row in Minv]
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
    return [sum([coeff[i]*basis[i] for i in range(d)]) if type(coeff) != str else coeff for coeff in coeffs]

def nf_elts_to_lists(elts, inv_basis):
    d = len(inv_basis)
    def to_list(elt):
        if (type(elt) == int) or (type(elt) == Integer):
            return [elt] + [0 for i in range(d-1)]
        elif type(elt) == str:
            return elt
        return list(elt)
    return [ list(sum([to_list(elt)[i]*inv_basis[i] for i in range(d)])) if type(elt) != str else elt for elt in elts]

def apply_to_nf_elt(func):
    def wrapper(*args):
        ev = args[-1]
        elt = args[-2]
        basis, inv_basis = get_nf_basis(ev)
        nf_elt = nf_lists_to_elements([elt], basis)[0]
        ret = func(*(list(args[:-2]) + [nf_elt]))
        return nf_elts_to_lists([ret], inv_basis)[0]
    return wrapper
