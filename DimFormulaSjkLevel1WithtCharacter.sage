load('DimFormulaSjkLevel1WithoutCharacter.sage')

def dimension_V(j, k):
    if (k % 2) == 0:
        d1 = 2^(-1) * dimension_cusp_forms(1, k+j/2)
        d2 = 2^(-1) * sum([dimension_cusp_forms(1, k+j-2*a) * dimension_cusp_forms(1, k+2*a) for a in range(j/2+1)])

    if (k % 2) == 1:
        d1 = - 2^(-1) * dimension_cusp_forms(1, k+j/2)
        d2 = 2^(-1) * sum([dimension_cusp_forms(1, k-1+j-2*a) * dimension_cusp_forms(1, k+1+2*a) for a in range(j/2)])

    return d1 + d2

def dimension_cusp_forms_sp4z_with_charac(j, k):
    d1 = dimension_cusp_forms_sp4z(j, k+5)
    d2 = dimension_V(j, k+5)

    return d1-d2

    

    

    

