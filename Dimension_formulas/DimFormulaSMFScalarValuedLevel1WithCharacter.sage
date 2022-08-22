load('DimFormulaPlusMinusNewFormGamma0.sage')

def dim_SV_sp4Z_even_weight_with_charac(k):
    """
    Compute the dimension of M_{0,k}(Sp(4,Z),eps) for
    k even
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = t^30
    denom = (1-t^4) * (1-t^6) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]
    return d

def dim_splitting_SV_even_weight_charac(k):
    """
    Compute the dimension of M_{0,k}(Sp(4,Z),eps),
                             S_{0,k}(Sp(4,Z),eps),
                             E_{0,k}(Sp(4,Z),eps),
                             KE_{0,k}(Sp(4,Z),eps),
                             SK_{0,k}(Sp(4,Z),eps),
                             S^{gen}_{0,k}(Sp(4,Z),eps)
    for k even.                         
    """
    k=ZZ(k)
    dtotal = dim_SV_sp4Z_even_weight_with_charac(k)
    deis = 0
    dklingeneis = 0
    dsaitokurokawa = 0
    dgenuine = dtotal
    dcusp = dsaitokurokawa+dgenuine
    dnoncusp = deis + dklingeneis
    L=[dtotal,dnoncusp,dcusp,deis,dklingeneis,dsaitokurokawa,dgenuine] 
    return L 

def dim_SV_sp4Z_odd_weight_with_charac(k):
    """
    Compute the dimension of M_{0,k}(Sp(4,Z),eps) for
    k odd
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = t^5
    denom = (1-t^4) * (1-t^6) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]
    return d

def dim_splitting_SV_odd_weight_charac(k):
    """
    Idem dim_splitting_SV_even_weight_Charac but for k odd
    """
    k=ZZ(k)
    dtotal = dim_SV_sp4Z_odd_weight_with_charac(k)
    deis = 0
    dklingeneis = 0
    dsaitokurokawa = dimension_new_cusp_forms_plus(2,2*(k-1))
    dgenuine = dtotal-dsaitokurokawa
    dcusp = dtotal
    dnoncusp = deis + dklingeneis
    L=[dtotal,dnoncusp,dcusp,deis,dklingeneis,dsaitokurokawa,dgenuine] 
    return L 

def dim_splitting_SV_All_weight_charac(k):
    """
    Put everything together
    """
    k = ZZ(k)
    if k == 0 :
        return [0,0,0,0,0,0,0]
    
    if k == 1 :
        return [0,0,0,0,0,0,0]
    
    if k == 2 :
        return [0,0,0,0,0,0,0]

    if (k % 2) == 1:
        return dim_splitting_SV_odd_weight_charac(k)

    if (k % 2) == 0 and  k > 2:
        return dim_splitting_SV_even_weight_charac(k)




        
    

    

