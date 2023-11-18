from sage.all import (Matrix, PolynomialRing, QQ, NumberField, prime_range, prime_divisors, divisors, is_square, ZZ, sqrt, sign)

def integer_squarefree_part(n):
    """ returns the squarefree part of the integer n (uses factor rather than calling pari like sage 9.3+ does) """
    return sign(n)*prod([p**(e%2) for p, e in ZZ(n).factor()])

def get_nf_basis(ev, quad_special=False):
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
    
    if (dim == 2) and quad_special:
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
    return [sum([coeff[i]*basis[i] for i in range(d)]) if type(coeff) != str else coeff for coeff in coeffs]

def nf_elts_to_lists(elts, inv_basis):
    d = len(inv_basis)
    def to_list(elt):
        if type(elt) == int:
            return [elt] + [0 for i in range(d-1)]
        elif type(elt) == str:
            return elt
        return list(elt)
    return [ list(sum([to_list(elt)[i]*inv_basis[i] for i in range(d)])) if type(elt) != str else elt for elt in elts]

def add_nf_data(f):
    pol = Qx([int(c) for c in f['field_poly'].split()[1:]])
    F = NumberField(pol, name = "a")
    V, from_V, to_V = F.vector_space()
    a = F.gens()[0]
    coeffs = [int(c) for c in pol.coefficients(sparse=False)]
    f['field_poly'] = coeffs
    max_p_idx = len(f['lambda_p'])
    if (F.degree() <= max_deg):
        # Tnis take too long
        # hecke_ring = F.order(f['lambda_p'])
        index = 0
        p_idx = 0
        nbound = 0
        num_gens = 0
        is_init_H = False
        while (type(f['lambda_p'][p_idx]) == str):
            p_idx += 1                    
            while ((index != 1) and (p_idx < max_p_idx)):
                print("p_idx = ", p_idx)
                gens = [x for x in f['lambda_p'][:p_idx+1] if type(x) != str]              
                mod_gens = [to_V(x) for x in gens]
                ambient = ZZ**V.dimension()
                W = ambient.span(mod_gens)
                if (W.rank() == F.degree()):
                    is_init_H = True
                    H = F.order(gens)
                    new_index = H.index_in(F.ring_of_integers())
                    if (new_index != index):
                        index = new_index
                        nbound = p_idx
                p_idx += 1
            if (is_init_H):
                f['hecke_ring'] = H
                f['hecke_ring_index'] = index
                f['hecke_ring_generator_nbound'] = nbound
            else:
                f['hecke_ring'] = F.ring_of_integers()
                f['hecke_ring_index'] = 1
                f['hecke_ring_generator_nbound'] = 0

            basis = f['hecke_ring'].basis()
            _ = f.pop('hecke_ring')
            f['hecke_ring_power_basis'] = False
            mat = Matrix([list(b) for b in basis])
            f['hecke_ring_denominators'] = [row.denominator() for row in mat]
            f['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
            f['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
            f['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat**(-1)]   
            inv_coeff_data = zip(f['hecke_ring_inverse_numerators'], f['hecke_ring_inverse_denominators'])
            inv_basis = [sum([nums[i] * a**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]
            f['lambda_p'] = nf_elts_to_lists(f['lambda_p'], inv_basis)
            f['lambda_p_square'] = nf_elts_to_lists(f['lambda_p_square'], inv_basis)
        else:
            f['trace_lambda_p'] = ['NULL' if f['lambda_p'][i] == 'NULL' else f['lambda_p'][i].trace() for i in range(len(f['lambda_p']))]
            f['trace_lambda_p_square'] = ['NULL' if f['lambda_p_square'][i] == 'NULL' else f['lambda_p_square'][i].trace() for i in range(len(f['lambda_p_square']))]
            f['lambda_p'] = []
            f['lambda_p_square'] = []
