def siegel_eis_Hecke_p(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series at p
    """
    k = ZZ(k)
    p = ZZ(p)
    l=(1+p^(k-1))*(1+p^(k-2))
    return l

def siegel_eis_Hecke_p_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(200):
        evs[p] = siegel_eis_Hecke_p(p,k)
    return evs

def siegel_eis_Hecke_p_square(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series at p^2
    """
    k = ZZ(k)
    p = ZZ(p)
    l=1+p^(4*k-6)+p^(k-1)+p^(k-2)+2*p^(2*k-3)+p^(3*k -5)+p^(3*k-4)+p^(2*k-2)
    return l

def siegel_eis_Hecke_p_square_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square(p,k)
    return evs

def siegel_eis_Hecke_p_square_0(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series for the Hecke operator 
    T_{2,0}(p^2)
    """
    k = ZZ(k)
    p = ZZ(p)
    l=1+p^(4*k-6)+2*p^(2*k-3)+p^(3*k-4)+p^(k-1)+p^(2*k-2)-p^(2*k-4)-p^(3*k-6)-p^(k-3)
    return l

def siegel_eis_Hecke_p_square_0_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square_0(p,k)
    return evs

def siegel_eis_Hecke_p_square_1(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series for the Hecke operator 
    T_{2,1}(p^2)
    """
    k = ZZ(k)
    p = ZZ(p)
    l=p^(2*k-4)-p^(2*k-6)+p^(3*k-6)+p^(3*k-5)+p^(k-3)+p^(k-2)
    return l

def siegel_eis_Hecke_p_square_1_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square_1(p,k)
    return evs

def siegel_eis_Hecke_p_square_2(p,k):
    """
    Compute the Hecke eigenvalues of the Siegel Eisenstein series for the Hecke operator 
    T_{2,2}(p^2)
    """
    k = ZZ(k)
    p = ZZ(p)
    l=p^(2*k-6)
    return l

def siegel_eis_Hecke_p_square_2_list(k):
    k = ZZ(k)
    evs = {}
    for p in prime_range(sqrt(200)):
        evs[p] = siegel_eis_Hecke_p_square_2(p,k)
    return evs


def Hecke_Eigenvalues_Siegel_Eisenstein_Series(k):
    """
    Compute                     
    """
    k=ZZ(k)
    L={}
    L['lambda_p'] = siegel_eis_Hecke_p_list(k)
    L['lambda_p_square'] = siegel_eis_Hecke_p_square_list(k)
    L['lambda_p_square_0'] = siegel_eis_Hecke_p_square_0_list(k)
    L['lambda_p_square_1'] = siegel_eis_Hecke_p_square_1_list(k)
    L['lambda_p_square_2'] = siegel_eis_Hecke_p_square_2_list(k)
    return L 

def Hecke_Eigenvalues_Siegel_Eisenstein_Series_All(k,j,e,prime_bound=200):
    empty_dic = {p : 0 for p in prime_range(prime_bound)}
    hecke_types = ['p', 'p_square', 'p_square_0', 'p_square_1', 'p_square_2']
    if (j != 0) or (e != 0) or (k < 4):
       return {'lambda_' + hecke : empty_dic for hecke in hecke_types}
    return Hecke_Eigenvalues_Siegel_Eisenstein_Series(k)
	



