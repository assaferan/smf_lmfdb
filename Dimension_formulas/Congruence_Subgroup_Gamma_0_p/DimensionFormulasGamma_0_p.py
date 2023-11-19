from sage.all import (Integers, PolynomialRing, QQ, dimension_cusp_forms, Gamma0, divisors, prime_divisors, prod, is_even, is_odd, is_squarefree, ConstantFunction)

'''
Here we implement the dimension formulas for spaces of modular forms 
on the congruence subgroup Gamma_0[p]. For cusp forms, we follow 
the paper of Wakatsuki:
Dimension formulas for spaces of vector-valued Siegel cusp forms of degree two
Journal of Number Theory 132 (2012) 200–253
'''


def Sum_Char_Pol(p,f,ch):
      Zp = Integers(p)
      ZZ = Integers()
      Zp_x = PolynomialRing(Zp,"x")
      x = Zp_x.gen()
      L = f.roots(ring = Zp)
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
    if a == -3:
        if p == 3:
            return 0
        if (p % 3 == 1):
            return 1
        return -1
    if a == 3:
        return IK_legendre(-1,p) * IK_legendre(-3,p)    

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
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
    S = Sum_Char_Pol(p, x**2+1, ch)
    a = [(-1)**(j//2),-1,-(-1)**(j//2), 1]
    aa = mux(a,k)
    A = QQ(2)**(-2) * aa * S 
    b = [(k-2)*(-1)**(j//2), -(j+k-1), -(k-2)*(-1)**(j//2), j+k-1]
    bb = mux(b,k)
    B = QQ(2)**(-5) * QQ(3)**(-1)* bb * S*(p+1) 
    return B-A
    

def H4(k,j,p,ch):
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
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
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
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
      Zp = Integers(p)
      Zp_x = PolynomialRing(Zp,"x")
      x = Zp_x.gen()
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
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
    S1 = Sum_Char_Pol(p, x**2+x+1, ch)
    S2 = Sum_Char_Pol(p, x**2+1, ch)
    A = QQ(2)**(-2)*QQ(3)**(-1)*C8(k,j)*S1*S2
    return A
      

def H9(k,j,p,ch):
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
    S1 = Sum_Char_Pol(p, x**2+x+1, ch)
    S2 = Sum_Char_Pol(p, x**2-x+1, ch)
    if p ==2:
       return QQ(2)**(-1)*QQ(3)**(-1)*C9(k,j)
    return QQ(3)**(-2)*C9(k,j)*S1*S2

def H10(k,j,p,ch):
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
    S = Sum_Char_Pol(p, x**4+x**3+x**2+x+1, ch) 
    A = QQ(5)**(-1)*C10(k,j)*S
    return A
     
def H11(k,j,p,ch):
     Zp = Integers(p)
     Zp_x = PolynomialRing(Zp,"x")
     x = Zp_x.gen()
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
    Zp = Integers(p)
    Zp_x = PolynomialRing(Zp,"x")
    x = Zp_x.gen()
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
    
    Example #2 
    e = DirichletGroup(3, QQ).0 cf p. 251
    can do 
    e = DirichletGroup(p,CyclotomicField(p-1)).0
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

'''
For j=0, k=3 or 4 and chi trivial, we use Theorems 2.2 and 2.6
of Ibukiyama:
Dimension Formulas of Siegel Modular Forms of Weight 3 
and Supersingular Abelian Surfaces(Revised version).
We start with the case k=3 (Theorem 2.2) 
'''

def case1(p):
    if (p % 8) == 1:
       return -QQ(2)**(-1)
    if (p % 8) == 3:
       return -QQ(2)**(-2)
    if (p % 8) == 5:
       return -QQ(2)**(-2)
    if (p % 8) == 7:
       return 0
                      
def case2(p):
    if (p % 5) == 1:
       return -4*QQ(5)**(-1)
    if (p % 5) == 2:
       return 0
    if (p % 5) == 3:
       return 0
    if (p % 5) == 4:
       return 0
    if p == 5:
       return -QQ(5)**(-1)
       
       
def cusp_form_gamma_0_p_k_is_3_dim(p):
    '''
    This implements the formula of Theorem 2.2
    Returns the dimension of S_{3}(Gamma_0(p))
    p prime 
    '''
    if p == 2:
        return  0
    if p == 3:
        return  0
    if p >= 5:
        a = (p+1)*(p**2+1)*QQ(2880)**(-1)-7*((p+1)**2)*QQ(576)**(-1)+55*(p+1)*QQ(288)**(-1)
        b = 1+ IK_legendre(-1,p)  
        c = 1+ IK_legendre(-3,p)
        d = (p-23)*c*QQ(36)**(-1)
        e = (2*p-25)*b*QQ(96)**(-1)
        f = -b*c*QQ(12)**(-1)
        return a+d+e+f+case1(p)+case2(p)
 
def case3(p):
       if (p % 8) == 1:
         return QQ(2)**(-1)
       if (p % 8) == 3:
         return QQ(2)**(-2)
       if (p % 8) == 5:
         return QQ(2)**(-2)
       if (p % 8) == 7:
         return 0

def case4(p):
      if (p % 12) == 1:
       return QQ(3)**(-1)
      if (p % 12) == 3:
       return QQ(2)**(-2)*QQ(3)**(-1)
      if (p % 12) == 5:
       return QQ(2)**(-1)*QQ(3)**(-1)
      if (p % 12) == 7:
       return QQ(2)**(-1)*QQ(3)**(-1)
      if (p % 12) == 11:
       return 0
       
               
def cusp_form_gamma_0_p_k_is_4_dim(p):
    '''
    This implements the formula of Theorem 2.6
    Returns the dimension of S_{4}(Gamma_0(p))
    p prime.
    Unfortunately, the formula given there is wrong, need to look at
    Hashimoto: 
    The dimension of the spaces of cusp forms on Siegel upper half-plane of degree two. I, 1983 
    Theorem 7.1 to get the right expression, thank you Ibukiyama-san for fixing it!
    '''
    if p == 2:
        return  0
    if p == 3:
        return  1
    if p == 5:
        return  1    
    if p >= 7:
        a = (p+1)*(p**2+1)*QQ(576)**(-1)+QQ(7)*((p+1)**2)*QQ(192)**(-1)-QQ(11)*(p+12)*QQ(288)**(-1)
        b = (p-1)*IK_legendre(-3,p)*QQ(36)**(-1)   
        c = (2*p-41)*IK_legendre(-1,p)*QQ(96)**(-1)
        d = IK_legendre(3,p)*QQ(12)**(-1)
        return a+b+c-d+case3(p)+case4(p)       

'''
For j=0, k=2 and chi trivial, we use Theorems 1.3 
C. POOR and D. S. YUEN
Dimensions of Cusp Forms for Gamma_0(p) in Degree Two and Small Weights
Abh. Math. Sem. Univ. Hamburg 77 (2007), 59-80
This is fo p = 2, ...,41
'''
def cusp_form_gamma_0_p_k_is_2_dim(p):
      if p == 2 or p == 3 or p == 5 or p == 7 or p == 13:
       return 0
      if p == 11 or p == 17 or p == 19: 
       return 1
      if p == 23 or p == 29 or p == 31: 
       return 3
      if p  == 37:
       return 2
      if p  == 41:
       return 6

'''
For j=0, k=1  and chi trivial, we use 
T. IBUKIYAMA and N. SKORUPPA, 
A vanishing theorem for Siegel modular forms of weight one. 
Abh. Math. Sem. Univ. Hamburg 77 (2007), 229-235
We have 
S_{1}(Gamma_0(N))=0 for any N
'''

def cusp_form_gamma_0_N_k_is_1_dim(N):
      return 0

def scalar_valued_form_gamma_0_p_cusp_dim(p,k):
    ''' Returns the dimension of the cuspidal subspace of S_{k}(Gamma_0(p)) for a prime p>2 and k>5 

      Verify example for prime 2<p<12, 4<k<21 from Hashimoto: The dimension of the spaces of cusp forms on Siegel upper half-plane of degree two. I, 1983, Table 7-11]:
      [[scalar_valued_form_gamma_0_p_cusp_dim(p,k) for k in range(5,21)] for p in range(3,12) if is_prime(p)]
      [[0, 2, 0, 5, 0, 10, 0, 16, 0, 23, 1, 35, 3, 47, 4, 61],
       [0, 5, 0, 13, 0, 25, 3, 44, 6, 66, 16, 100, 25, 136, 45, 188],
       [0, 11, 0, 26, 5, 56, 15, 95, 28, 145, 58, 222, 97, 312, 143, 417],
       [2, 31, 9, 80, 33, 164, 80, 288, 158, 462, 278, 694, 444, 991, 666, 1365]]
    '''  
    f = ConstantFunction(1)
    return cusp_form_gamma_0_p_chi_dim(k,0,p,f)

def scalar_valued_form_gamma_0_p_mod_dim(p,k):
    ''' Returns the dimension of the modular form space of M_{k}(Gamma_0(p)) for a prime p>2 for k>5 (Using Hashimoto 1983 and Wakatsuki 2012, Böcherer-Ibukiyama, 2012)

            dim M_{k}(Gamma_0(p))= S_{k}(Gamma_0(p)) for k odd
         
            For dim M_{k,j}(Gamma_0(p)), we have the follwoing example for 2<p<12 and 4<k<21.
            [[scalar_valued_form_gamma_0_p_mod_dim(p,k) for k in range(5,21)] for p in range(3,12) if is_prime(p)]
            [[0, 7, 0, 10, 0, 17, 0, 25, 0, 32, 1, 46, 3, 60, 4, 74],
             [0, 10, 0, 22, 0, 34, 3, 57, 6, 79, 16, 117, 25, 153, 45, 209],
             [0, 20, 0, 35, 5, 69, 15, 112, 28, 162, 58, 243, 97, 337, 143, 442],
             [2, 42, 9, 95, 33, 183, 80, 311, 158, 489, 278, 725, 444, 1026, 666, 1404]]
    '''
    f = ConstantFunction(1)
    if (k % 2) == 1 and k>=5: 
        return cusp_form_gamma_0_p_chi_dim(k,0,p,f)
    if k  == 3 : 
        return cusp_form_gamma_0_p_k_is_3_dim(p)
    if k == 1:
        return 0              
    if (k % 2) == 0 and k>=6:
        return cusp_form_gamma_0_p_chi_dim(k,0,p,f)+2*dimension_cusp_forms(Gamma0(p),k)+3       
    if k == 4:
        return cusp_form_gamma_0_p_k_is_4_dim(p)+2*dimension_cusp_forms(Gamma0(p),4)+3
    if k == 2:
        return cusp_form_gamma_0_p_k_is_2_dim(p)+2*dimension_cusp_forms(Gamma0(p),2)+1
          
def vector_valued_form_gamma_0_p_cusp_dim(p,k,j):
    ''' Returns the dimension of the cuspidal subspace of S_{k,j}(Gamma_0(p)) for a prime p for k>5 and j>=0 even (Uses Wakatsuki 2012) 

    The following verify the example for p=3, 4<k<20 and 0<j<10 for j even from Wakatsuki JNT 2012 (page 251):
    [[vector_valued_form_gamma_0_p_cusp_dim(3,k,j) for k in range(5,20)] for j in range(0,9) if is_even(j)]
    [[0, 2, 0, 5, 0, 10, 0, 16, 0, 23, 1, 35, 3, 47, 4],
    [0, 2, 0, 7, 3, 16, 6, 26, 12, 44, 24, 67, 37, 92, 54],
    [0, 5, 3, 14, 10, 29, 20, 49, 36, 79, 61, 116, 90, 163, 130],
    [4, 11, 11, 27, 25, 51, 46, 84, 74, 128, 116, 187, 168, 258, 232],
    [7, 18, 19, 42, 43, 77, 74, 123, 118, 187, 181, 269, 256, 365, 349]]
    '''
    f = ConstantFunction(1)
    return cusp_form_gamma_0_p_chi_dim(k,j,p,f)

def vector_valued_form_gamma_0_p_mod_dim(p,k,j):
    ''' Returns the dimension of the modular form space of M_{k,j}(Gamma_0(p)) for a prime p for k>5 and j>=0 even (Uses Wakatsuki 2012, Böcherer-Ibukiyama, 2012)
    dim M_{k,j}(Gamma_0(p))= S_{k,j}(Gamma_0(p)) for k odd
    For dim M_{k,j}(Gamma_0(p)) when k is even, we have the follwoing example for p=3,4<k<20 and 0=<j<10.
    [[vector_valued_form_gamma_0_p_mod_dim(3,k,j) for k in range(5,20) if (k % 2) == 0] for j in range(0,9) if is_even(j)]
    [[4, 7, 14, 22, 29, 43, 57],
    [4, 11, 22, 32, 52, 77, 102],
    [9, 20, 35, 57, 89, 126, 175],
    [17, 33, 59, 94, 138, 199, 272],
    [24, 50, 87, 133, 199, 283, 379]]
    '''
    f = ConstantFunction(1)
    if (k % 2) == 1 and k>=5:
        return cusp_form_gamma_0_p_chi_dim(k,j,p,f)
    if (k % 2) == 0 and k>=5:
        return cusp_form_gamma_0_p_chi_dim(k,j,p,f)+2*dimension_cusp_forms(Gamma0(p),j+k)     

