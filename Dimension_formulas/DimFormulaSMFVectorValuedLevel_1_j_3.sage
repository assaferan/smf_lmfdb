def dim_VV_sp4Z_j_3_without_charac(j):
    """
    Compute the dimension of S_{j,3}(Sp(4,Z))
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=j+1)
    num = t^36
    denom = (1-t^6) * (1-t^8) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(j)
    d = f.list()[j]
    return d

def dim_VV_sp4Z_j_3_with_charac(j):
    """
    Compute the dimension of S_{j,3}(Sp(4,Z),eps)
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=j+1)
    num = t^6
    denom = (1-t^6) * (1-t^8) * (1-t^10) * (1-t^12)
    f = num / denom
    k = ZZ(j)
    d = f.list()[j]
    return d    