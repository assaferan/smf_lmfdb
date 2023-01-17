from sage.all import (QQ, divisors, prime_divisors, prod, is_even, is_odd, is_squarefree)

'''
Here we implement the dimenion formulas for cuspidal spaces of paramodular forms,
following the paper of Wakatsuki:
Dimension formulas for spaces of vector-valued Siegel cusp forms of degree two
Journal of Number Theory 132 (2012) 200–253
'''

def Sum_Char_Pol(p,f,ch):
      L = L = f.roots(x, ring = Integers(p))
      l = len(L)
      S = sum([ch(ZZ(L[n][0])) for n in range(0,l)])
      return S 

def IK_legendre(a,p):
    if a == -1:
        if p == 2:
            return 0
        if (p % 4 == 1):
            return 1
        return -1

def mux(a,n):
    m = len(a)
    return a[n % m]

def C8(k,j):
    arr_table = [
        [1,0,0,-1,-1,-1,-1,0,0,1,1,1],
        [-1,1,0,1,1,0,1,-1,0,-1,-1,0],
        [1,-1,0,0,-1,1,-1,1,0,0,1,-1],
        [-1,0,0,-1,1,-1,1,0,0,1,-1,1],
        [1,1,0,1,-1,0,-1,-1,0,-1,1,0],
        [-1,-1,0,0,1,1,1,1,0,0,-1,-1]
        ]
    return mux(arr_table[(j % 12)//2], k)

def C9(k,j):
    arr_table = [
        [1,0,0,-1,0,0],
        [-1,1,0,1,-1,0],
        [0,-1,0,0,1,0]
        ]
    return mux(arr_table[(j % 6) // 2], k)

def C10(k,j):
    arr_table = [
        [1,0,0,-1,0],
        [-1,1,0,0,0],
        [0,0,0,0,0],
        [0,0,0,1,-1],
        [0,-1,0,0,1]
        ]
    return mux(arr_table[(j % 10) // 2], k)

def C11(k,j):
    arr_table = [
        [1,0,0,-1],
        [-1,1,0,0],
        [-1,0,0,1],
        [1,-1,0,0]
        ]
    return mux(arr_table[(j % 8) // 2], k)

C = { n : eval('C'+str(n)) for n in range(8,12)}

def H1(k,j,p,ch):
    C = QQ(2)**(-7) * QQ(3)**(-3) * QQ(5)**(-1)
    D = (j+1)*(k-2)*(j+k-1)*(j+2*k-3)*(p+1)*(p**2+1)
    E = QQ(2)**(-4) * QQ(3)**(-2)*(j+1)*(j+2*k-3)*(p+1)
    F = QQ(2)**(-2) * QQ(3)**(-1)*(j+1)
    return C*D-E+F
    
    
def H2(k,j,p,ch):
    C = QQ(2)**(-7) * QQ(3)**(-2)
    D = (j+k-1)*(k-2)*(-1)**k*ch(-1)
    E = 57 if p == 2 else 7*(p+1)**2 
    F = QQ(2)**(-3)*QQ(3)**(-1)*(-1)**k*(j+2*k-3)*(p+1)*ch(-1)
    G = QQ(2)**(-4)*(-1)**k*ch(-1)*(7-IK_legendre(-1,p)) 
    return C*D*E-F+G

def H3(k,j,p,ch):
    S = Sum_Char_Pol(p, x**2+1, ch)
    a = [(-1)**(j//2),-1,-(-1)**(j//2), 1]
    aa = mux(a,k)
    A = QQ(2)**(-2) * aa * S 
    b = [(k-2)*(-1)**(j//2), -(j+k-1), -(k-2)*(-1)**(j//2), j+k-1]
    bb = mux(b,k)
    B = QQ(2)**(-5) * QQ(3)**(-1)* bb * S*(p+1) 
    return B-A
    

def H4(k,j,p,ch):
    S = Sum_Char_Pol(p, x**2+x+1, ch)
    a = [1,-1,0]
    aa = mux(a,k)
    b= [1,0,-1]
    bb = mux(b,j+k)
    A = QQ(2)**(-1) * QQ(3)**(-2)* (aa+bb) * S 
    c = [1,0,1]
    cc = mux(c,k)
    d = [0,-1,-1]
    dd = mux(d,j+k)
    B = QQ(2)**(1) * QQ(3)**(-2)* (cc+dd) * S
    e = [j+k-1,-(j+k-1),0]
    ee = mux(e,k)
    f = [k-2,0,-(k-2)]
    ff = mux(f,j+k) 
    C = QQ(2)**(-3) * QQ(3)**(-3)* (ee+ff) * (p+1) * S
    return C-A-B

def H5(k,j,p,ch):
    S = Sum_Char_Pol(p, x**2-x+1, ch)
    a = [-1,-1,0,1,1,0]
    aa = mux(a,k)
    b= [1,0,-1,-1,0,1]
    bb = mux(b,j+k)
    A = QQ(2)**(-1) * QQ(3)**(-1)* (aa+bb) * S 
    c = [-(j+k-1),-(j+k-1),0,j+k-1,j+k-1,0]
    cc = mux(c,k)
    d = [k-2,0,-(k-2),-(k-2),0,k-2]
    dd = mux(d,j+k)
    B = QQ(2)**(-3) * QQ(3)**(-2)* (cc+dd) * (p+1) * S
    return B-A

def H62(k,j,p,ch):
    A = QQ(2)**(-3)*(-1)**(j//2)*(2+ch(-1)*(1+IK_legendre(-1,p)))
    B = QQ(2)**(-7)*QQ(3)**(-1)*(-1)**(j//2)*(j+2*k-3)*23
    C = QQ(2)**(-7)*(-1)**(j//2+k)*(j+1)*3
    return B+C-A
    
def H6p(k,j,p,ch):
    a = ch(-1)*(1+IK_legendre(-1,p))
    A = QQ(2)**(-3)*(-1)**(j//2)*(2+a)
    B = QQ(2)**(-7)*QQ(3)**(-1)*(-1)**(j//2)*(j+2*k-3)*5*(p+1+a)
    C = QQ(2)**(-7)*(-1)**(j//2+k)*(j+1)*(p+1+a)
    return B+C-A    
 
def H6(k,j,p,ch):
     if p == 2:
            return H62(k,j,p,ch)
     return  H6p(k,j,p,ch)

def H7(k,j,p,ch):
      S = Sum_Char_Pol(p, x**2+x+1, ch)
      a = [1,-1,0]
      aa = mux(a,j)
      b = [0,1,-1]
      bb = mux(b,j+2*k)
      if p == 3:
         return -QQ(2)**(-1)*QQ(3)**(-1)*aa*(2+S)+QQ(2)**(-1)*QQ(3)**(-3)*(j+2*k-3)*aa*7+QQ(2)**(-2)*QQ(3)**(-3)*(j+1)*bb
      if (p % 3) == 1:
         return -QQ(2)**(-1)*QQ(3)**(-1)*aa*(2+S)+QQ(2)**(-1)*QQ(3)**(-3)*(j+2*k-3)*aa*(p+1+S)+QQ(2)**(-2)*QQ(3)**(-3)*(j+1)*bb*(p-1+S**2)   
      return -QQ(2)**(-1)*QQ(3)**(-1)*aa*(2+S)+QQ(2)**(-1)*QQ(3)**(-3)*(j+2*k-3)*aa*(p+1+S)+QQ(2)**(-2)*QQ(3)**(-3)*(j+1)*bb*(p+1)

def H8(k,j,p,ch): 
      S1 = Sum_Char_Pol(p, x**2+x+1, ch)
      S2 = Sum_Char_Pol(p, x**2+1, ch)
      A = QQ(2)**(-2)*QQ(3)**(-1)*C8(k,j)*S1*S2
      return A
      

def H9(k,j,p,ch):     
     S1 = Sum_Char_Pol(p, x**2+x+1, ch)
     S2 = Sum_Char_Pol(p, x**2-x+1, ch)
     if p ==2:
        return QQ(2)**(-1)*QQ(3)**(-1)*C9(k,j)
     return QQ(3)**(-2)*C9(k,j)*S1*S2

def H10(k,j,p,ch):     
     S = Sum_Char_Pol(p, x**4+x**3+x**2+x+1, ch) 
     A = QQ(5)**(-1)*C10(k,j)*S
     return A
     
def H11(k,j,p,ch):
     S = Sum_Char_Pol(p, x**2+1, ch)
     if p == 2:
        return QQ(2)**(-3)*C11(k,j)
     if (p % 8) == 1: 
        return  QQ(2)**(-3)*C11(k,j)*(2*ch(-1)+S)
     if (p % 8) == 3: 
        return  QQ(2)**(-3)*C11(k,j)*2*ch(-1)  
     if (p % 8) == 5: 
        return  QQ(2)**(-3)*C11(k,j)*S
     if (p % 8) == 7: 
        return  0
        
def H12(k,j,p,ch):
    S = Sum_Char_Pol(p, x**2+x+1, ch)
    a = [0,-1,1]
    aa = mux(a,j+2*k)
    b = ch(-1)*(1+IK_legendre(-1,p))
    A = QQ(2)**(-2) * QQ(3)**(-1)*(-1)**(j//2)*aa*(b+S) 
    return A

H = { n : eval('H'+str(n)) for n in range(1,13)}

def cusp_form_gamma_0_p_chi_dim(k,j,p,ch):
    '''
    This implements the formula of Theorem 7.4
    Returns the dimension of S_{k,j}(Gamma_0(p),chi)
    for k>=5, j>=0, j even, p prime and chi Dirichlet character modulo p
    
    Example #1 (e = DirichletGroup(3, QQ).0 so e^2 is the trivial character) cf p. 251
    >>>[[cusp_form_gamma_0_p_chi_dim(k,j,3,e^2) for k in range(4,20)] for j in range(0,9) if is_even(j)]# doctest: +NORMALIZE_WHITESPACE 
    [[1, 0, 2, 0, 5, 0, 10, 0, 16, 0, 23, 1, 35, 3, 47, 4],
     [0, 0, 2, 0, 7, 3, 16, 6, 26, 12, 44, 24, 67, 37, 92, 54],
     [1, 0, 5, 3, 14, 10, 29, 20, 49, 36, 79, 61, 116, 90, 163, 130],
     [3, 4, 11, 11, 27, 25, 51, 46, 84, 74, 128, 116, 187, 168, 258, 232],
     [5, 7, 18, 19, 42, 43, 77, 74, 123, 118, 187, 181, 269, 256, 365, 349]]
    
    Example #2 (e = DirichletGroup(3, QQ).0) cf p. 251
    Warning with DirichletGroup(p, CC).0 got some floating point numbers... not cool)
    >>>[[cusp_form_gamma_0_p_chi_dim(k,j,3,e) for k in range(4,20)] for j in range(0,9) if is_even(j)]# doctest: +NORMALIZE_WHITESPACE 
    [[0, 1, 0, 4, 0, 7, 0, 12, 0, 20, 1, 29, 1, 39, 4, 55],
     [0, 1, 0, 5, 1, 10, 3, 21, 10, 36, 17, 53, 28, 79, 47, 112],
     [0, 2, 2, 9, 6, 20, 14, 38, 29, 63, 47, 95, 74, 139, 111, 191],
     [1, 7, 7, 19, 17, 38, 33, 66, 59, 106, 94, 156, 138, 220, 199, 301],
     [3, 10, 14, 29, 30, 56, 56, 98, 97, 154, 148, 223, 214, 314, 304, 426]]
    '''
    H_part = sum([H[i](k,j,p,ch) for i in range(1,13)])
    dim = H_part 
    return dim

