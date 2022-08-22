def dimension_new_cusp_forms_plus(N,k):
    """
    Compute the dimension of the plus space of 
    new cusp forms on Gamma_0(N)
    """
    k = ZZ(k)

    if k==2:
        raise ValueError("we don't know")

    if (k % 8) == 0 or (k % 8)== 2:
        return (1/2)*dimension_new_cusp_forms(Gamma0(N),k)+1/2 

    elif (k % 8) == 4 or (k % 8)== 6:
        return (1/2)*dimension_new_cusp_forms(Gamma0(N),k)


def dimension_new_cusp_forms_minus(N,k):
    """
    Compute the dimension of the minus space of 
    new cusp forms on Gamma_0(N)
    """
    k = ZZ(k)

    if k==2:
        raise ValueError("we don't know")

    if (k % 8) == 0 or (k % 8)== 2:
        return (1/2)*dimension_new_cusp_forms(Gamma0(N),k)-1/2 

    elif (k % 8) == 4 or (k % 8)== 6:
        return (1/2)*dimension_new_cusp_forms(Gamma0(N),k)