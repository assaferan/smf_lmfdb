def dim_SV_sp4Z_even_weight_trivial_charac(k):
    """
    Compute the dimension of M_{0,k}(Sp(4,Z)) for
    k even
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = 1
    denom = (1-t^4) * (1-t^6) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]
    return d

def dim_SV_sp4Z_even_weight_klingen_eis(k):
    """
    Compute the dimension of KE_{0,k}(Sp(4,Z)) for
    k even where KE_{0,k}(Sp(4,Z)) denotes the
    subspace of M_{0,k}(Sp(4,Z)) generated by 
    Klingen-Eisenstein series
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = t^12
    denom = (1-t^4) * (1-t^6)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]   
    return d

def dim_SV_sp4Z_even_weight_saito_kurokawa(k):
    """
    Compute the dimension of SK_{0,k}(Sp(4,Z)) for
    k even where SK_{0,k}(Sp(4,Z)) denotes the
    subspace of M_{0,k}(Sp(4,Z)) generated by 
    Saito-Kurokawa lifts
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = t^10+t^12
    denom = (1-t^4) * (1-t^6)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]   
    return d  

def dim_SV_sp4Z_even_weight_genuine(k):
    """
    Compute the dimension of S^{gen}_{0,k}(Sp(4,Z)) for
    k even where S^{gen}_{0,k}(Sp(4,Z)) denotes the orthogonal
    complement (Peterssom) of S_{0,k}(Sp(4,Z)) 
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = t^20*(1+t^2+t^4-t^12-t^14)
    denom = (1-t^4) * (1-t^6)* (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]   
    return d

def dim_splitting_SV_even_weight(k):
    """
    Compute the dimension of M_{0,k}(Sp(4,Z)),
                             S_{0,k}(Sp(4,Z)),
                             E_{0,k}(Sp(4,Z)),
                             KE_{0,k}(Sp(4,Z)),
                             SK_{0,k}(Sp(4,Z)),
                             S^{gen}_{0,k}(Sp(4,Z))
    for k even.                         
    """
    k=ZZ(k)
    dtotal = dim_SV_sp4Z_even_weight_trivial_charac(k)
    deis = 1
    dklingeneis = dim_SV_sp4Z_even_weight_klingen_eis(k)
    dsaitokurokawa = dim_SV_sp4Z_even_weight_saito_kurokawa(k)
    dgenuine = dim_SV_sp4Z_even_weight_genuine(k)
    dcusp = dsaitokurokawa+dgenuine
    dnoncusp = deis + dklingeneis
    L={}
    L['degree'] = 2
    L['family'] = 'S'
    L['level'] = 1
    L['weight'] = [k,0]
    L['char_orbit'] = 0
    L['total_dim'] = dtotal
    L['cusp_dim'] = dcusp
    L['eis_dim'] = dnoncusp  
    L['eis_F_dim'] = deis  
    L['eis_Q_dim'] = dklingeneis 
    L['cusp_P_dim'] = dsaitokurokawa
    L['cusp_Y_dim'] = 0
    L['cusp_G_dim'] = dgenuine
    return L 

def dim_SV_sp4Z_odd_weight_trivial_charac(k):
    """
    Compute the dimension of M_{0,k}(Sp(4,Z)) for
    k odd
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=k+1)
    num = t^35
    denom = (1-t^4) * (1-t^6) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(k)
    d = f.list()[k]
    return d

def dim_splitting_SV_odd_weight(k):
    """
    idem as dim_splitting_SV_even_weight(k) but for k odd
    """
    k=ZZ(k)
    dtotal = dim_SV_sp4Z_odd_weight_trivial_charac(k)
    deis = 0
    dklingeneis = 0
    dsaitokurokawa = 0
    dgenuine = dtotal
    dcusp = dgenuine
    dnoncusp = deis + dklingeneis
    L={}
    L['degree'] = 2
    L['family'] = 'S'
    L['level'] = 1
    L['weight'] = [k,0]
    L['char_orbit'] = 0
    L['total_dim'] = dtotal
    L['cusp_dim'] = dcusp
    L['eis_dim'] = dnoncusp 
    L['eis_F_dim'] = deis  
    L['eis_Q_dim'] = dklingeneis 
    L['cusp_P_dim'] = dsaitokurokawa
    L['cusp_Y_dim'] = 0
    L['cusp_G_dim'] = dgenuine
    return L 



def dim_splitting_SV_All_weight(k):
    """
    Put everything together
    """
    k = ZZ(k)
    if k == 0 :
        return {'degree': 2, 'family': 'S', 'level': 1, 'weight': [0, 0], 'char_orbit' : 0, 'total_dim': 1, 'cusp_dim': 0, 'eis_dim': 1, 'eis_F_dim': 1,
 'eis_Q_dim': 0,
 'cusp_P_dim': 0, 'cusp_Y_dim': 0 , 'cusp_G_dim': 0}

    if k == 2 :
         return {'degree': 2, 'family': 'S', 'level': 1, 'weight': [2, 0], 'char_orbit' : 0, 'total_dim': 0, 'cusp_dim': 0, 'eis_dim': 0, 'eis_F_dim': 0,
 'eis_Q_dim': 0,
 'cusp_P_dim': 0, 'cusp_Y_dim': 0 , 'cusp_G_dim': 0}

    if (k % 2) == 1:
        return dim_splitting_SV_odd_weight(k)

    if (k % 2) == 0 and  k > 2:
        return dim_splitting_SV_even_weight(k) 


        
  