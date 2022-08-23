def dim_VV_sp4Z_j_3_without_charac(j):
    """
    Compute the dimension of S_{3,j}(Sp(4,Z))
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=j+1)
    num = t^36
    denom = (1-t^6) * (1-t^8) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(j)
    d = f.list()[j]
    return d

def dimension_cusp_forms_SP4Z(J, k):
    """
    Uses Tsushima's formula from Theorem 4 of 'An explicit dimension formula
    for the spaces of generalized automorphic forms with respect to Sp(2,Z)'
    Proc. Japan Acad. 59 Ser A (1983).

    Tsushima proves the correctness of his formula for (j = 0 and  k >= 4)
    or (j >= 1 and k >= 5), but Bergstroem-Faber-van der Geer prove that it
    holds for (j >= 0 and k >= 4), see page 97 of 'Siegel modular forms of
    degree three and the cohomology of local systems' Sel. Math. New Ser.
    (2014) 20:83-124.
    """
    #if (J % 2) == 1:
    #    return ZZ(0)
    j = ZZ(J/2)
    #if j < 0:
    #    raise ValueError("j cannot be negative")

    k = ZZ(k)
    #if k < 4:
    #    raise ValueError("not implemented for k < 4")

    res = 2^(-7) * 3^(-3) * 5^(-1) * (2*j+1) * (k-2) * (2*j+k-1) * (2*j+2*k-3)
    res += - 2^(-5) * 3^(-2) * (2*j+1) * (2*j+2*k-3)
    res += 2^(-4) * 3^(-1) * (2*j+1)
    res += (-1)^k * (2^(-7) * 3^(-2) * 7 * (k-2) * (2*j+k-1) - 2^(-4) * 3^(-1) * (2*j+2*k-3) + 2^(-5) * 3)
    res += (-1)^j * (2^(-7) * 3^(-1) * 5 * (2*j+2*k-3) - 2^(-3))
    res += (-1)^k * (-1)^j * 2^(-7) * (2*j+1)

    i = CyclotomicField(4).gen()
    rho = CyclotomicField(3).gen()
    omega = CyclotomicField(5).gen()
    sigma = CyclotomicField(12).gen()

    res += (i^k * (2^(-6) * 3^(-1) * i * (2*j+k-1) - 2^(-4) * i)).trace()
    res += ((-1)^k * i^j * 2^(-5) * (i+1)).trace()
    res += (i^k * (-1)^j * (2^(-6) * 3^(-1) * (k-2) - 2^(-4))).trace()
    res += ((-i)^k * i^j * 2^(-5) * (i+1)).trace()
    res += ((-1)^k * rho^j * 3^(-3) * (rho+1)).trace()
    res += (rho^k * rho^j * 2^(-2) * 3^(-4) * (2*rho+1) * (2*j+1)).trace()
    res += - (rho^k * (-rho)^j * 2^(-2) * 3^(-2) * (2*rho+1)).trace()
    res += ((-rho)^k * rho^j * 3^(-3)).trace()
    res += (rho^j * (2^(-1) * 3^(-4) * (1-rho) * (2*j+2*k-3) - 2^(-1) * 3^(-2) * (1-rho))).trace()
    res += (rho^k * (2^(-3) * 3^(-4) * (rho+2) * (2*j+k-1) - 2^(-2) * 3^(-3) * (5*rho+6))).trace()
    res += - ((-rho)^k * (2^(-3) * 3^(-3) * (rho+2) * (2*j+k-1) - 2^(-2) * 3^(-2) * (rho+2))).trace()
    res += (rho^k * (rho^2)^j * (2^(-3) * 3^(-4) * (1-rho) * (k-2) + 2^(-2) * 3^(-3) * (rho-5))).trace()
    res += ((-rho)^k * (rho^2)^j * (2^(-3) * 3^(-3) * (1-rho) * (k-2) - 2^(-2) * 3^(-2) * (1-rho))).trace()
    res += (omega^k * (omega^4)^j * 5^(-2)).trace()
    res += - (omega^k * (omega^3)^j * 5^(-2) * omega^2).trace()
    res += ((sigma^7)^k * (-1)^j * 2^(-3) * 3^(-2) * (sigma^2+1)).trace()
    res += - ((sigma^7)^k * (sigma^8)^j * 2^(-3) * 3^(-2) * (sigma+sigma^3)).trace()
    return ZZ(res)

    




def dim_splitting_VV_All_weight(k,j):
    """
    Put everything together
    """
    k = ZZ(k)
    j = ZZ(j)
    
    if (j % 2) == 1 or k == 0 or k == 1:
        return {
                      'degree': 2, 
                      'type': 'S', 
                      'level': 1, 
                      'weight': [k, j], 
                      'char_orbit' : 0, 
                      'total_dim': 0, 
                      'cusp_dim': 0,
                      'eis_dim': 0,
                      'eis_P_dim': 0,
                      'eis_Q_dim': 0,
                      'cusp_P_dim': 0, 
                      'cusp_Y_dim': 0, 
                      'cusp_G_dim': 0}




    elif k == 2 :
        if j > 53:
                  raise ValueError("we do not know")
        return {
                'degree': 2, 
                'type': 'S', 
                'level': 1, 
                'weight': [2, j], 
                'char_orbit' : 0, 
                'total_dim': 0, 
                'cusp_dim': 0, 
                'eis_dim': 0, 
                'eis_P_dim': 0,
                'eis_Q_dim': 0,
                'cusp_P_dim': 0, 
                'cusp_Y_dim': 0, 
                'cusp_G_dim': 0}

    elif k == 3 and (j % 2) == 0:  
        return {
                'degree': 2, 
                'type': 'S', 
                'level': 1, 
                'weight': [3, j], 
                'char_orbit' : 0, 
                'total_dim': dim_VV_sp4Z_j_3_without_charac(j), 
                'cusp_dim': dim_VV_sp4Z_j_3_without_charac(j), 
                'eis_dim': 0, 
                'eis_P_dim': 0,
                'eis_Q_dim': 0,
                'cusp_P_dim': 0,
                'cusp_Y_dim': 0, 
                'cusp_G_dim': dim_VV_sp4Z_j_3_without_charac(j)}

    elif k > 3 :
        return {
'degree': 2, 
'type': 'S', 
'level': 1, 
'weight': [k, j], 
'char_orbit' : 0, 
'total_dim': dimension_cusp_forms_SP4Z(j,k)+ dimension(CuspForms(SL2Z,k+j)), 
'cusp_dim': dimension_cusp_forms_SP4Z(j,k), 
'eis_dim': dimension(CuspForms(SL2Z,k+j)), 
'eis_P_dim': 0,
'eis_Q_dim': dimension(CuspForms(SL2Z,k+j)),
'cusp_P_dim': 0,
'cusp_Y_dim': 0, 
'cusp_G_dim': dimension_cusp_forms_SP4Z(j,k)
}




        
  