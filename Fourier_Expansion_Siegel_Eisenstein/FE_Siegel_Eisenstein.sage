def MakeSum(f,D0,k):
    L = f.divisors()
    s = sum([moebius(d)*kronecker(D0,d)*d^(k-2)*sigma(f/d,2*k-3) for d in L])
    return s
    
def MakeL(D,k):
    D0 = fundamental_discriminant(D)
    f = D/D0
    f = ZZ(f)
    sqf = sqrt(f)
    sqf = ZZ(sqf)
    if  f == 1 :
        return quadratic_L_function__exact(2-k,D)
    else:
        return quadratic_L_function__exact(2-k,D0)*MakeSum(sqf,D0,k)

def MakeCohenNumber(D,k):
    if (D % 4) == 1 or (D % 4) == 2:
       return 0
    if D ==0:
       return zeta(3-2*k)
    else:
       return MakeL(-D,k)  


def MakeCoeffEisNonZero(k,L):
    B = 2*k*(2*k-2)/(bernoulli(k)*bernoulli(2*k-2))
    m = L[0]
    r = L[1]
    n = L[2]
    g = gcd(gcd(r,m),n)
    LL = g.divisors()
    s = sum([d^(k-1)* MakeCohenNumber((4*m*n-r^2)/d^2,k) for d in LL])
    C = B*s
    R.<X,Y,u> = LaurentPolynomialRing(QQ,3)
    Co = C*X^m*u^r*Y^n
    return Co

def MakeCoeffEis(k,L):
    if L == [0,0,0]:
       return 1
    else:
       return MakeCoeffEisNonZero(k,L)



def MakeFEEis(k,N):
    L = [[m,r,n] for m in range(N+1) for r in range(-N-1,N+1) for n in range(N+1) if m*n-r^2/4>=0]
    R.<X,Y,u> = LaurentPolynomialRing(QQ,3)
    s = sum([MakeCoeffEis(k,A) for A in L])
    return s


    

