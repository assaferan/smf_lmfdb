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




        
  