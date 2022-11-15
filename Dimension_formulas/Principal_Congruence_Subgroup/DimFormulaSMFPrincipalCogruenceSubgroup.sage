"""
We implement the dimension formulas for modular forms on the 
principal congruence subgroup of level N of Sp(4,Z)
For regular wieghts i.e. j=0 and k>=4 or j>=1 and k>=5, 
we use the dimension formula due to Tsushima (see the note).
"""


def dim_cusp_form_GammaN(k,j,N):
    """
    Compute the dimension of S_{k,j}(Gamma[N]) for
     N>=3 and  
    j=0 and k>=4 or j>=1 and k>=5
    """
    L = divisors(N)
    LP = [x for x in L if x.is_prime() == true]
    V = [(1-x^(-2))*(1-x^(-4)) for x in LP]
    mu = prod(V)
    """
    Notice that mu*N^10 is the index of Gamma[N] in Sp(4,Z)
    """
    A = (j+1)*(j+k-1)*(j+2*k-3)*(k-2)*N^10/1080
    B = (j+1)*(j+2*k-3)*N^8/576
    C = (j+1)*N^7/96
    d = (A-B+C)*mu
    """
    d1 = (j+1)*(8*(j+k-1)*(j+2*k-3)*(k-2)+15*(6-(j+2*k-3)*N)/N^3)*mu*N^10/8640 
    d1 = d (see the note)
    """
    return d

"""
For N=2, we use the paper of Berström, Faber and van der Geer.
The conjectures there are now proved thanks to the work of Rösner (see his PhD)    
"""

def dim_cusp_form_Gamma2(k,j):
    """
    Compute the dimension of S_{k,j}(Gamma[2]) for
    j=0 and k>=4 or j>=1 and k>=5
    and j even otherwise this is zero since -1 belongs to Gamma[2]
    """
    A = (k-2)*j^3
    B = 3*(k^2-3*k-3)*j^2
    C = (2*k^3-6*k^2+(15*(-1)^k-27)*k+77+75*(-1)^(k+1))*j
    D = 2*k^3+(15*(-1)^k-9)*k^2+(135*(-1)^(k+1)-17)*k+84+300*(-1)^k
    d = (A+B+C+D)/24
    return d

"""
Dimension of spaces of scalar-valued modular forms on Gamma[2]
"""

def Dim_Scalar_Valued_P_2(k):
    """
    Compute the dimension of M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    k = ZZ(k)
    num = 1+t^5-t^8-t^13
    denom = (1-t^2)^5
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s6_SV(k):
    """
    Compute the multiplicity of the trivial rep. of S_6 in  M_{0,k}(Gamma[2]).
    This is the same as the dimension of M_{0,k}(Sp(4,Z)).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = 1+t^35
    denom = (1-t^4) * (1-t^6) * (1-t^10) * (1-t^12)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s51_SV(k):
    """
    Compute the multiplicity of the rep. s[5,1] of S_6 in M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^11*(1+t)
    denom = ((1-t^4) * (1-t^6))^2
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s42_SV(k):
    """
    Compute the multiplicity of the rep. s[4,2] of S_6 in M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^4*(1+t^15)
    denom = (1-t^2)*((1-t^4)^2)*(1-t^10)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s411_SV(k):
    """
    Compute the multiplicity of the rep. s[4,1,1] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^11*(1+t^4)
    denom = (1-t) * (1-t^4) * (1-t^6) * (1-t^12)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s33_SV(k):
    """
    Compute the multiplicity of the rep. s[3,3] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^7*(1+t^13)
    denom = (1-t^2)*(1-t^4)*(1-t^6)*(1-t^12)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s321_SV(k):
    """
    Compute the multiplicity of the rep. s[3,2,1] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^8*(1-t^8)
    denom = ((1-t^2)^2)*(1-t^5)*((1-t^6))^2
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s3111_SV(k):
    """
    Compute the multiplicity of the rep. s[3,1,1,1] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^6*(1+t^4+t^11+t^15)
    denom = (1-t^2)*(1-t^4)*(1-t^6)*(1-t^12)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s222_SV(k):
    """
    Compute the multiplicity of the rep. s[2,2,2] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^2*(1+t^23)
    denom = (1-t^2)*(1-t^4)*(1-t^6)*(1-t^12)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s2211_SV(k):
    """
    Compute the multiplicity of the rep. s[2,2,1,1] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^9
    denom = (1-t^2)*((1-t^4)^2)*(1-t^5)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s21111_SV(k):
    """
    Compute the multiplicity of the rep. s[2,1,1,1,1] of S_6 in  M_{0,k}(Gamma[2]).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^6*(1+t^11)
    denom = ((1-t^4)^2)*((1-t^6)^2)
    f = num / denom
    d = f.list()[k]
    return d

def Mult_s111111_SV(k):
    """
    Compute the multiplicity of the rep. s[1,1,1,1,1,1] of S_6 in  M_{0,k}(Gamma[2]).
    This is the same as the dimension of M_{0,k}(Sp(4,Z),eps).
    """
    R.<t> = PowerSeriesRing(ZZ, default_prec=2*k+1)
    num = t^5*(1+t^25)
    denom = (1-t^4)*(1-t^6)*(1-t^10)*(1-t^12)
    f = num / denom
    d = f.list()[k]
    return d

def Check_1(k):
    """ 
    This checks the fact that dim M_{0,k}(Gamma[2])
    equals the sum of the multiplicities of each irrep. of
    S_6 times their dimension:
    irrep.           dim
    s[6]             1    
    s[5,1]           5    
    s[4,2]           9
    s[4,1,1]         10    
    s[3,3]           5
    s[3,2,1]         16    
    s[3,1,1,1]       10    
    s[2,2,2]         5
    s[2,2,1,1]       9    
    s[2,1,1,1,1]     5    
    s[1,1,1,1,1,1]   1
    """
    d1 = Dim_Scalar_Valued_P_2(k)
    d2 = Mult_s6_SV(k)+5*Mult_s51_SV(k)+9*Mult_s42_SV(k)+10*Mult_s411_SV(k)+5*Mult_s33_SV(k) 
    d2 += 16*Mult_s321_SV(k)+10*Mult_s3111_SV(k)+5*Mult_s222_SV(k)+9*Mult_s2211_SV(k)+5*Mult_s21111_SV(k)   
    d2 += Mult_s111111_SV(k) 
    d = [d1-d2, [Mult_s6_SV(k),
Mult_s51_SV(k),Mult_s42_SV(k),Mult_s411_SV(k),Mult_s33_SV(k),Mult_s321_SV(k),Mult_s3111_SV(k),Mult_s222_SV(k),Mult_s2211_SV(k),Mult_s21111_SV(k),Mult_s111111_SV(k)]]
    return d

    
"""
To do for N>=3
- j=0, k=0, 1, 2, 3
- j>=1 k=0, 1, 2, 3, 4
- Dimnension of the spaces of non-cusp forms.
- Dimension of Yoshida lifts
- Dimension of Saito-Kurokawa lifts
- Other lifts? (cube transfer?)
- ...
"""    
