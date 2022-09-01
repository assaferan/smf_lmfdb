from sage.all import (QQ, divisors, prime_divisors, prod, is_even, is_odd, is_squarefree)
from lmfdb import db

'''
Here we implement the dimenion formulas for cupidal spaces of paramodular forms,
followinf Ibukiyama-Kitiyama
'''

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

def IK_jacobi(a,N):
    return prod([IK_legendre(a,p) for p in prime_divisors(N)])

def omega(N):
    return len(prime_divisors(N))

def mux(a,n):
    m = len(a)
    return a[n % m]

def NN(m,n,N):
    return [p for p in prime_divisors(N) if p % n == m % n]

def C3(k,j):
    a = [(k-2)*(-1)**(j//2), -(j+k-1), -(k-2)*(-1)**(j//2), j+k-1]
    return mux(a,k)

def C4(k,j):
    return (j+k-1)*mux([1,-1,0],k) + (k-2)*mux([1,0,-1],j+k)

def C5(k,j):
    return (j+k-1)*mux([-1,-1,0,1,1,0], k) + (k-2)*mux([1,0,-1,-1,0,1],j+k)

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

C = { n : eval('C'+str(n)) for n in range(3,12) if n not in [6,7]}

def H1(k,j,N):
    C = QQ(2)**(-7) * QQ(3)**(-3) * QQ(5)**(-1)
    D = (j+1)*(k-2)*(j+k-1)*(j+2*k-3)
    return C*D*prod([p**2 + 1 for p in prime_divisors(N)])

def H2(k,j,N):
    C = QQ(2)**(omega(N)-8) * QQ(3)**(-2)
    D = (j+k-1)*(k-2)*(-1)**k
    E = 11 if is_even(N) else 14
    return C*D*E

def H3(k,j,N):
    C = QQ(2)**(omega(N)-6) * QQ(3)**(-1)
    E = 5 if is_even(N) else 2
    return C * C3(k,j) * E

def H4(k,j,N):
    C = QQ(2)**(omega(N)-3) * QQ(3)**(-3)
    E = 5 if (N % 3 == 0) else 1
    return C * C4(k,j) * E

def H5(k,j,N):
    C = QQ(2)**(omega(N)-3) * QQ(3)**(-2)
    return C * C5(k,j)

def H6_1(N):
    a = 1/(QQ(2**5 * 3))
    b = 1/(QQ(2**7 * 3))
    prod_a = prod([p + IK_legendre(-1,p) for p in prime_divisors(N)])
    prod_b = prod([1 + p*IK_legendre(-1,p) for p in prime_divisors(N)])
    return a*prod_a + b * prod_b

def H6_2(N):
    a = 1/(QQ(2**5 * 3))
    b = 1/(QQ(2**7 * 3))
    prod_a = prod([p + IK_legendre(-1,p) for p in prime_divisors(N)])
    prod_b = prod([1 + p*IK_legendre(-1,p) for p in prime_divisors(N)])
    return a*prod_a - b * prod_b

def H6(k,j,N):
    a = (-1)**(j//2) * (2*k+j-3)
    b = (-1)**(j//2 + k) * (j+1)
    return a * H6_1(N) + b * H6_2(N)

def H7_1(N):
    a = 1/(QQ(2**3 * 3**2))
    b = 1/(QQ(2**3 * 3**3))
    prod_a = prod([p + IK_legendre(-3,p) for p in prime_divisors(N)])
    prod_b = prod([1 + p*IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a + b * prod_b

def H7_2(N):
    a = 1/(QQ(2**3 * 3**2))
    b = 1/(QQ(2**3 * 3**3))
    prod_a = prod([p + IK_legendre(-3,p) for p in prime_divisors(N)])
    prod_b = prod([1 + p*IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a - b * prod_b 

def H7(k,j,N):
    return H7_1(N)*mux([1,-1,0],j)*(2*k+j-3)+H7_2(N)*mux([0,1,-1],j+2*k)*(j+1)

def H8(k,j,N):
    return QQ(2)**(omega(N)-2)*QQ(3)**(-1) * C8(k,j)

def H9(k,j,N):
    E = 1 if is_even(N) else 4
    return QQ(2)**(omega(N)-2)*QQ(3)**(-2) * C9(k,j) * E

def H10(k,j,N):
    E = 2**(len(NN(1,5,N)+NN(-1,5,N))) if len(NN(2,5,N) + NN(3,5,N)) == 0 else 0
    return QQ(5)**(-1) * C10(k,j) * E

def H11(k,j,N):
    E = 2**(len(NN(1,8,N)+NN(-1,8,N))) if len(NN(3,8,N) + NN(5,8,N)) == 0 else 0
    return QQ(2)**(-3) * C11(k,j) * E

def H12_1(N):
    a = 1/(QQ(2**3 * 3))
    b = 1/(QQ(2**3 * 3))
    prod_a = prod([1 + IK_legendre(3,p) for p in prime_divisors(N)])
    prod_b = prod([IK_legendre(-1,p) + IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a - b * prod_b

def H12_2(N):
    a = 1/(QQ(2**3 * 3))
    b = 1/(QQ(2**3 * 3))
    prod_a = prod([1 + IK_legendre(3,p) for p in prime_divisors(N)])
    prod_b = prod([IK_legendre(-1,p) + IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a + b * prod_b 

def H12(k,j,N):
    return H12_1(N)*mux([1,-1,0],j)*(-1)**(j//2+k)+H12_2(N)*mux([0,-1,1],j+2*k)*(-1)**(j//2)

H = { n : eval('H'+str(n)) for n in range(1,13)}

def I1(k,j,N):
    C = QQ(2)**(-3)*QQ(3)**(-1)
    D = j+1
    return C*D*prod([p+1 for p in prime_divisors(N)])

def I2(k,j,N):
    return -QQ(2)**(omega(N)-4)*QQ(3)**(-1)*(j+1)

def I3(k,j,N):
    return -QQ(2)**(omega(N)-5)*QQ(3)**(-2)*N*(j+1)*(2*k+j-3)

def I4(k,j,N):
    return QQ(2)**(omega(N)-5)*(-1)**k * (4 - IK_jacobi(-1,N))

def I5(k,j,N):
    return -QQ(2)**(omega(N)-4)*QQ(3)**(-1)*(-1)**k*(2*k+j-3)

def I6(k,j,N):
    return -QQ(2)**(omega(N)-3)*mux([(-1)**(j//2),-1,-(-1)**(j//2),1],k)

def I7(k,j,N):
    arr_table = [
        [[3,-3,0],[3,0,-3]],
        [[5,-1,4],[1,-4,-5]],
        [[1,-5,-4],[5,4,-1]]
        ]
    arrs = arr_table[N % 3]
    return -QQ(2)**(omega(N)-2)*QQ(3)**(-2)*(mux(arrs[0],k) + mux(arrs[1],j+k))

def I8(k,j,N):
    arrs = [[-1,-1,0,1,1,0], [1,0,-1,-1,0,1]]
    return -QQ(2)**(omega(N)-2)*QQ(3)**(-1)*(mux(arrs[0],k) + mux(arrs[1],j+k))

def I9(k,j,N):
    C = -QQ(2)**(-3)*(-1)**(j//2)
    return C*prod([1 + IK_legendre(-1,p) for p in prime_divisors(N)])

def I10(k,j,N):
    C = -QQ(2)**(-1)*QQ(3)**(-1)*mux([1,-1,0],j)
    return C*prod([1 + IK_legendre(-3,p) for p in prime_divisors(N)])

I = { n : eval('I'+str(n)) for n in range(1,11)}

def paramodular_cusp_dim(k,j,N):
    '''
    This implements the formula from Main Theorem 3.1
    Returns the dimension of the cuspidal subspace S_{k,j}(K(N))

    Example #1 (total dim. p. 618):
    >>> [[paramodular_cusp_dim(k,0,N) for k in range(3,12)]
    ... for N in [6,3,2,1]] # doctest: +NORMALIZE_WHITESPACE 
    [[0, 0, 0, 1, 1, 2, 3, 4, 5],
     [0, 0, 0, 1, 0, 1, 1, 2, 1],
     [0, 0, 0, 0, 0, 1, 0, 1, 1],
     [0, 0, 0, 0, 0, 0, 0, 1, 0]]

    
    '''
    assert is_squarefree(N)
    if is_odd(j):
        return 0
    assert ((j >= 2) and (k >= 5)) or ((j == 0) and (k >= 3))
    H_part = sum([H[i](k,j,N) for i in range(1,13)])
    I_part = sum([I[i](k,j,N) for i in range(1,11)])
    delta = (j == 0)*(k == 3)
    return H_part + I_part + delta

# Too tired to implement the trace formula here as well
# for now we cheat
def classical_minus_cusp_dim(k,N):
    '''
    Returns the dimension of S_k^{new, -}(N).
    We do it by querying the LMFDB for the total dimension and the
    dimension of the plus subspace

    Example #1 (p. 618):
    >>> [[classical_minus_cusp_dim(2*k-2,N) for k in range(3,12)]
    ... for N in [6,3,2,1]] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 1, 0, 1, 1, 1],
     [0, 0, 0, 1, 0, 1, 1, 1, 1],
     [0, 0, 0, 0, 0, 1, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 1, 0]]
    '''
    label = str(N)+'.'+str(k) + '.a'
    f = db.mf_newspaces.lookup(label, ['dim', 'plus_dim'])
    if f['dim'] == 0:
        return 0
    return f['plus_dim'] if is_odd(k//2) else f['dim']-f['plus_dim']

def paramodular_new_cusp_dim(k,j,N):
    '''
    This implements the formula in Proposition 4.4
    Returns the dimension of the cuspidal subspace S_{k,j}^{new}(K(N))

    Example #1 (new dim. p. 619):
    >>> [[paramodular_new_cusp_dim(k,0,N) for k in range(3,12)]
    ... for N in [6,3,2,1]] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 1, 0, 2, 2, 3],
     [0, 0, 0, 1, 0, 1, 1, 1, 1],
     [0, 0, 0, 0, 0, 1, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 1, 0]]
    
    '''
    pcd = paramodular_cusp_dim
    cmcd = classical_minus_cusp_dim
    main_term = sum([(-2)**omega(M)*pcd(k,j,N//M) for M in divisors(N)])
    # still have to deal with the j = 0 case
    if j == 0:
        j0_term = sum([(-1)**omega(M)*cmcd(2*k-2,N//M) for M in divisors(N) if M > 1])
        return main_term - j0_term
    return main_term
        
