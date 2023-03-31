from sage.all import (QQ, divisors, prime_divisors, prod, is_even, is_odd, is_squarefree)
from lmfdb import db

'''
Here we implement the dimenion formulas for cuspidal spaces of paramodular forms,
following the paper [IK] Ibukiyama, Kitiyama - "Dimension formulas of paramodular
forms of squarefree level and comparison with inner twist", J. Math. Soc. Japan, Vol. 69, No. 2 (2017)

Although conjectural at the time, the formulas are now a consequence of the preprint
[RW] Rosner, Weissauer - "Global liftings between inner forms of GSp(4)".

!! TODO - the formulas are missing the forms with other Atkin-Lehner signs (as in Rama-Tornaria) !!
Add that.
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

    Example #2 (total dim., scalar valued, p.666):
    >>> [[paramodular_cusp_dim(k,0,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE 
    [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 3],
     [0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 0, 2, 1, 4, 1, 4, 2, 7],
     [0, 0, 0, 1, 0, 1, 1, 2, 1, 4, 1, 4, 3, 6, 3, 10, 4, 11],
     [0, 0, 1, 1, 1, 2, 2, 4, 4, 6, 5, 9, 8, 13, 12, 18, 16, 25],
     [0, 0, 0, 1, 1, 2, 3, 4, 5, 10, 6, 11, 15, 20, 17, 30, 27, 40],
     [0, 1, 1, 2, 2, 4, 4, 7, 7, 11, 11, 16, 16, 24, 24, 33, 33, 45],
     [0, 1, 1, 2, 4, 6, 7, 12, 15, 21, 23, 32, 38, 52, 55, 72, 81, 103],
     [0, 1, 2, 3, 3, 6, 7, 12, 14, 20, 22, 32, 36, 48, 54, 69, 76, 97],
     [1, 2, 3, 5, 7, 10, 13, 19, 23, 31, 37, 48, 56, 72, 82, 102, 115, 140],
     [0, 1, 2, 3, 5, 10, 12, 19, 26, 36, 42, 58, 70, 92, 105, 132, 152, 189],
     [0, 1, 2, 5, 5, 10, 14, 21, 26, 40, 44, 62, 74, 96, 109, 144, 156, 197]]

    Example #3 (total dim., vector valued j=2, p. 667):
    >>> [[paramodular_cusp_dim(k,2,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE 
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 2, 0, 3],
     [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 4, 4, 6, 4, 9, 7, 12],
     [0, 0, 0, 0, 0, 1, 1, 2, 2, 4, 4, 8, 8, 13, 11, 18, 16, 26],
     [0, 0, 0, 0, 1, 2, 3, 5, 7, 10, 13, 19, 22, 31, 33, 44, 48, 63],
     [0, 0, 0, 1, 2, 4, 6, 11, 15, 20, 25, 38, 46, 59, 65, 87, 96, 121],
     [0, 0, 1, 1, 3, 5, 8, 12, 16, 22, 29, 39, 47, 61, 70, 88, 100, 124],
     [0, 0, 1, 3, 7, 12, 19, 30, 42, 56, 73, 99, 123, 154, 181, 226, 262, 316],
     [0, 1, 2, 4, 9, 14, 21, 31, 43, 57, 75, 97, 119, 151, 178, 217, 254, 304],
     [0, 1, 3, 5, 12, 19, 29, 43, 59, 79, 104, 134, 166, 208, 248, 301, 354, 422],
     [0, 1, 4, 8, 17, 27, 42, 63, 87, 115, 151, 197, 245, 304, 363, 443, 520, 618],
     [0, 1, 4, 8, 17, 29, 44, 65, 89, 121, 157, 205, 253, 318, 377, 461, 538, 646]]

    Example #4 (total dim., vector valued j=4, p. 667):
    >>> [[paramodular_cusp_dim(k,4,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE 
    [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 3, 1, 4, 2, 6],
     [0, 0, 0, 0, 1, 1, 1, 3, 3, 5, 4, 8, 9, 13, 11, 17, 18, 26],
     [0, 0, 0, 1, 0, 2, 3, 5, 5, 10, 8, 16, 17, 23, 24, 35, 33, 49],
     [0, 1, 1, 3, 4, 7, 10, 16, 18, 27, 31, 43, 50, 65, 72, 92, 102, 127],
     [0, 0, 1, 4, 6, 12, 17, 25, 34, 49, 54, 78, 95, 117, 135, 170, 191, 239],
     [0, 1, 1, 5, 7, 12, 18, 28, 34, 49, 59, 79, 95, 120, 138, 172, 196, 237],
     [0, 3, 5, 13, 22, 35, 50, 74, 96, 131, 161, 209, 256, 317, 370, 450, 524, 623],
     [1, 4, 6, 15, 22, 35, 51, 73, 93, 127, 157, 201, 245, 302, 355, 430, 498, 589],
     [1, 5, 9, 20, 31, 49, 71, 101, 131, 176, 220, 280, 342, 420, 497, 598, 696, 820],
     [1, 6, 12, 27, 45, 70, 101, 145, 191, 255, 319, 407, 500, 613, 725, 872, 1020, 1201],
     [0, 5, 11, 28, 42, 71, 104, 148, 194, 264, 326, 422, 515, 632, 750, 907, 1049, 1246]]

    Example #5 (total dim., vector valued j=6, p. 668):
    >>> [[paramodular_cusp_dim(k,6,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE 
    [[0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 1, 3, 2, 5, 3, 7, 4, 9],
     [0, 0, 0, 0, 1, 2, 2, 4, 6, 7, 7, 13, 14, 19, 19, 28, 28, 37],
     [0, 0, 1, 2, 2, 5, 7, 9, 12, 18, 18, 27, 31, 40, 43, 59, 59, 77],
     [0, 0, 3, 4, 7, 13, 17, 24, 33, 43, 52, 69, 81, 102, 118, 145, 162, 197],
     [0, 0, 3, 8, 13, 22, 33, 45, 62, 82, 97, 130, 159, 192, 223, 279, 312, 372],
     [1, 3, 7, 11, 18, 28, 38, 52, 69, 89, 109, 139, 166, 204, 238, 286, 327, 387],
     [0, 4, 12, 23, 42, 63, 91, 127, 170, 218, 274, 347, 422, 513, 605, 725, 840, 982],
     [3, 5, 18, 27, 44, 68, 94, 126, 170, 216, 270, 338, 409, 494, 587, 695, 804, 941],
     [4, 10, 25, 40, 64, 96, 134, 180, 239, 305, 381, 475, 575, 694, 823, 973, 1129, 1316],
     [2, 10, 29, 51, 86, 131, 185, 253, 339, 433, 543, 682, 829, 1001, 1188, 1410, 1639, 1910],
     [2, 10, 31, 55, 88, 137, 195, 263, 351, 455, 565, 710, 863, 1043, 1236, 1472, 1701, 1990]]
    
    '''
    assert is_squarefree(N)
    if is_odd(j):
        return 0
    # assert ((j >= 2) and (k >= 5)) or ((j == 0) and (k >= 3))
    assert k >= 3
    H_part = sum([H[i](k,j,N) for i in range(1,13)])
    I_part = sum([I[i](k,j,N) for i in range(1,11)])
    delta = (j == 0)*(k == 3)
    dim = H_part + I_part + delta
    assert dim >= 0
    return dim

def cmf_label(k,N):
    return str(N)+'.'+str(k) + '.a'

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
    f = db.mf_newspaces.lookup(cmf_label(k,N), ['dim', 'plus_dim'])
    if f['dim'] == 0:
        return 0
    return f['plus_dim'] if is_odd(k//2) else f['dim']-f['plus_dim']

# we could return them together, but we never call them with the same level
def classical_plus_cusp_dim(k,N):
    '''
    Returns the dimension of S_k^{new, +}(N).
    We do it by querying the LMFDB for the total dimension and the
    dimension of the plus subspace

    Example #1 (p. 618):
    >>> [[classical_plus_cusp_dim(2*k-2,N) for k in range(3,12)]
    ... for N in [6,3,2,1]] # doctest: +NORMALIZE_WHITESPACE
    [[1, 1, 1, 1, 2, 1, 2, 2, 2],
     [0, 1, 1, 1, 1, 2, 1, 2, 2],
     [0, 0, 1, 1, 0, 1, 1, 1, 1],
     [0, 0, 0, 0, 1, 0, 1, 0, 1]]
    
    '''
    f = db.mf_newspaces.lookup(cmf_label(k,N), ['dim', 'plus_dim'])
    if f['dim'] == 0:
        return 0
    return f['plus_dim'] if is_even(k//2) else f['dim']-f['plus_dim']

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

    Example #2 (new dim., scalar valued, p. 668):
    >>> [[paramodular_new_cusp_dim(k,0,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 3],
     [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 2, 1, 2, 2, 3],
     [0, 0, 0, 1, 0, 1, 1, 1, 1, 3, 1, 3, 3, 4, 3, 8, 4, 7],
     [0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 5, 8, 8, 11, 12, 16, 16, 21],
     [0, 0, 0, 0, 1, 0, 2, 2, 3, 4, 5, 5, 10, 9, 12, 12, 18, 19],
     [0, 1, 1, 2, 2, 4, 4, 6, 7, 10, 11, 15, 16, 22, 24, 31, 33, 41],
     [0, 1, 0, 1, 3, 3, 5, 7, 9, 12, 16, 18, 25, 29, 35, 40, 51, 57],
     [0, 1, 2, 3, 3, 6, 7, 11, 14, 19, 22, 31, 36, 46, 54, 67, 76, 93],
     [1, 2, 3, 5, 7, 10, 13, 18, 23, 30, 37, 47, 56, 70, 82, 100, 115, 136],
     [0, 0, 1, 1, 3, 4, 7, 10, 15, 19, 25, 32, 43, 50, 63, 73, 91, 106],
     [0, 1, 1, 3, 4, 7, 11, 15, 20, 28, 36, 45, 58, 70, 86, 102, 123, 144]]

    Example #3 (new dim., vector valued j = 2, p. 669):
    >>> [[paramodular_new_cusp_dim(k,2,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 2, 0, 3],
     [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 4, 2, 4, 5, 7, 6],
     [0, 0, 0, 0, 0, 1, 1, 2, 2, 4, 4, 6, 8, 9, 11, 14, 16, 20],
     [0, 0, 0, 0, 1, 2, 3, 5, 7, 10, 13, 17, 22, 27, 33, 40, 48, 57],
     [0, 0, 0, 1, 2, 2, 4, 5, 9, 10, 15, 18, 22, 29, 35, 41, 50, 57],
     [0, 0, 1, 1, 3, 5, 8, 12, 16, 22, 29, 37, 47, 57, 70, 84, 100, 118],
     [0, 0, 1, 3, 5, 8, 13, 18, 26, 34, 45, 57, 71, 88, 107, 128, 152, 178],
     [0, 1, 2, 4, 9, 14, 21, 31, 43, 57, 75, 95, 119, 147, 178, 213, 254, 298],
     [0, 1, 3, 5, 12, 19, 29, 43, 59, 79, 104, 132, 166, 204, 248, 297, 354, 416],
     [0, 1, 2, 6, 11, 17, 26, 37, 53, 69, 91, 115, 143, 178, 215, 257, 306, 358],
     [0, 1, 4, 8, 15, 23, 36, 51, 71, 93, 123, 155, 193, 238, 289, 345, 410, 480]]
    
    Example #4 (new dim., vector valued j = 4, p. 669):
    >>> [[paramodular_new_cusp_dim(k,4,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 3, 1, 4, 2, 6],
     [0, 0, 0, 0, 1, 1, 1, 1, 3, 3, 4, 4, 7, 7, 9, 9, 14, 14],
     [0, 0, 0, 1, 0, 2, 3, 3, 5, 8, 8, 12, 15, 17, 22, 27, 29, 37],
     [0, 1, 1, 3, 4, 7, 10, 14, 18, 25, 31, 39, 48, 59, 70, 84, 98, 115],
     [0, 0, 1, 2, 4, 6, 9, 13, 18, 23, 30, 38, 47, 57, 69, 82, 97, 113],
     [0, 1, 1, 5, 7, 12, 18, 26, 34, 47, 59, 75, 93, 114, 136, 164, 192, 225],
     [0, 1, 3, 7, 12, 19, 28, 40, 54, 71, 91, 115, 142, 173, 208, 248, 292, 341],
     [1, 4, 6, 15, 22, 35, 51, 71, 93, 125, 157, 197, 243, 296, 353, 422, 494, 577],
     [1, 5, 9, 20, 31, 49, 71, 99, 131, 174, 220, 276, 340, 414, 495, 590, 692, 808],
     [1, 4, 10, 17, 29, 44, 63, 87, 117, 151, 193, 241, 296, 359, 431, 510, 600, 699],
     [0, 3, 9, 20, 34, 53, 78, 110, 148, 194, 248, 312, 385, 468, 562, 669, 787, 918]]

    Example #5 (new dim., vector valued j = 6, p. 670):
    >>> [[paramodular_new_cusp_dim(k,6,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 1, 3, 2, 5, 3, 7, 4, 9],
     [0, 0, 0, 0, 1, 0, 2, 2, 4, 3, 5, 7, 10, 9, 13, 14, 20, 19],
     [0, 0, 1, 2, 2, 3, 7, 7, 10, 14, 16, 21, 27, 30, 37, 45, 51, 59],
     [0, 0, 3, 4, 7, 11, 17, 22, 31, 39, 50, 63, 77, 92, 112, 131, 154, 179],
     [0, 0, 1, 4, 7, 12, 15, 23, 30, 40, 51, 62, 77, 94, 111, 133, 154, 180],
     [1, 3, 7, 11, 18, 26, 38, 50, 67, 85, 107, 133, 162, 194, 232, 272, 319, 369],
     [0, 4, 6, 15, 26, 37, 53, 75, 96, 126, 160, 195, 240, 291, 343, 407, 476, 550],
     [3, 5, 18, 27, 44, 66, 94, 124, 168, 212, 268, 332, 405, 484, 581, 681, 796, 923],
     [4, 10, 25, 40, 64, 94, 134, 178, 237, 301, 379, 469, 571, 684, 817, 959, 1121, 1298],
     [0, 4, 15, 29, 48, 75, 105, 145, 193, 249, 315, 390, 477, 575, 686, 810, 945, 1098],
     [2, 10, 23, 43, 70, 105, 147, 201, 265, 341, 429, 530, 647, 779, 926, 1092, 1275, 1478]]
    
    '''
    pcd = paramodular_cusp_dim
    cmcd = classical_minus_cusp_dim
    main_term = sum([(-2)**omega(M)*pcd(k,j,N//M) for M in divisors(N)])
    # still have to deal with the j = 0 case
    if j == 0:
        j0_term = sum([(-1)**omega(M)*cmcd(2*k-2,N//M) for M in divisors(N) if M > 1])
        dim = main_term - j0_term
        assert dim >= 0
        return dim
    return main_term

# Here again we cheat even though there are implemented formulas
def Yoshida_lift_dim_sub(k,j,N,M):
    wt = [j+2, 2*k+j-2]
    lev = [M, N//M]
    lookup = db.mf_newspaces.lookup
    labels = [cmf_label(wt[i],lev[i]) for i in [0,1]]
    dims = [lookup(label, ['dim'])['dim'] for label in labels]
    return prod(dims)
    
def Yoshida_new_lift_dim(k,j,N):
    Ms = [M for M in divisors(N) if is_odd(omega(M))]
    return sum([Yoshida_lift_dim_sub(k,j,N,M) for M in Ms])

def Yoshida_lift_dim(k,j,N):
    '''
    Returns the dimension of the total space of Yoshida lifts, including old forms
    '''
    return sum([Yoshida_new_lift_dim(k,j,M) for M in divisors(N)])

def Saito_Kurokawa_new_lift_dim(k,N):
    # all the rest here are considered old
    return classical_minus_cusp_dim(2*k-2,N)

def Saito_Kurokawa_lift_dim(k,N):
    '''
    Returns the dimension of the total space of SK(=Gritsenko) lifts, including old forms

    Example #1 (dimension of Jacobi cusp forms, p. 666):
    >>> [[Saito_Kurokawa_lift_dim(k,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2],
     [0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 0, 2, 1, 3, 1, 3, 1, 4],
     [0, 0, 0, 1, 0, 1, 1, 2, 1, 3, 1, 3, 2, 4, 2, 5, 2, 5],
     [0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 3, 5, 4, 6, 5, 7, 5, 8],
     [0, 0, 0, 1, 1, 2, 2, 3, 3, 5, 3, 5, 5, 7, 5, 8, 6, 9],
     [0, 1, 1, 2, 2, 3, 3, 5, 4, 6, 5, 7, 6, 9, 7, 10, 8, 11],
     [0, 1, 1, 2, 3, 4, 4, 6, 6, 8, 7, 9, 9, 12, 10, 13, 12, 15],
     [0, 1, 2, 3, 3, 5, 5, 7, 7, 9, 8, 11, 10, 13, 12, 15, 13, 17],
     [1, 2, 3, 4, 5, 6, 7, 9, 9, 11, 11, 13, 13, 16, 15, 18, 17, 20],
     [0, 1, 2, 3, 4, 6, 6, 8, 9, 11, 10, 13, 13, 16, 15, 18, 17, 21],
     [0, 1, 2, 4, 4, 6, 7, 9, 9, 12, 11, 14, 14, 17, 16, 20, 18, 22]]
    
    '''
    return sum([Saito_Kurokawa_new_lift_dim(k,M) for M in divisors(N)])

# This is the total dimension of the (cuspidal) lifts (see (7) in p. 617)
def paramodular_new_lift_dim(k,j,N):
    # Yoshida lifts are not paramodular (see Saha & Schmidt), in fact they are also never of Siegel square-free level
    # Y_dim = Yoshida_new_lift_dim(k,j,N)
    SK_dim = Saito_Kurokawa_new_lift_dim(k,N) if j == 0 else 0
    # delta = (j == 0)*(k == 3)*(is_odd(omega(N)))
    delta = 0
    # return Y_dim + SK_dim + delta
    return SK_dim + delta

def paramodular_non_lift_new_cusp_dim(k,j,N):
    '''
    Returns the dimension of the cuspidal subspace of general type S_{k,j}^{new}_{(G)}(K(N))

    Example #1 (non-lifts, p. 619):
    >>> [[paramodular_non_lift_new_cusp_dim(k,0,N) for k in range(3,12)]
    ... for N in [6,3,2,1]] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 1, 1, 2],
     [0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0]]

    Example #2 (non-lifts, scalar valued, p. 669):
    >>> [[paramodular_non_lift_new_cusp_dim(k,0,N) for k in range(3,21)] for N
    ... in range(1,16) if is_squarefree(N)] # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 1, 5, 2, 4],
     [0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 4, 4, 7, 7, 11, 11, 15],
     [0, 0, 0, 0, 0, 0, 1, 1, 2, 3, 3, 4, 8, 7, 10, 10, 15, 17],
     [0, 0, 0, 0, 0, 1, 1, 2, 3, 5, 6, 9, 10, 15, 17, 23, 25, 32],
     [0, 0, 0, 0, 1, 2, 3, 4, 7, 9, 12, 15, 21, 24, 31, 35, 45, 52],
     [0, 0, 0, 0, 0, 1, 2, 5, 7, 11, 14, 21, 26, 35, 42, 54, 63, 78],
     [0, 0, 0, 1, 2, 4, 6, 10, 14, 20, 26, 35, 43, 56, 67, 84, 98, 118],
     [0, 0, 0, 0, 1, 2, 4, 7, 11, 15, 20, 27, 37, 44, 56, 66, 83, 98],
     [0, 0, 0, 1, 1, 4, 7, 10, 15, 22, 29, 38, 50, 61, 77, 92, 112, 133]]
    
    '''
    dim = paramodular_new_cusp_dim(k,j,N) - paramodular_new_lift_dim(k,j,N)
    assert dim >= 0
    return dim

def smf_dims_paramodular_squarefree(k,j,N):
    space = {}
    space['degree'] =  2 
    space['family'] = 'K' 
    space['level'] = N 
    space['weight'] = [k, j]
    space['char_orbit_index'] = 1
    space['cusp_dim'] = paramodular_cusp_dim(k,j,N)
    # space['cusp_Y_dim'] = Yoshida_lift_dim(k,j,N)
    # Yoshida lifts are not paramodular!!
    space['cusp_Y_dim'] = 0
    space['cusp_P_dim'] = Saito_Kurokawa_lift_dim(k,N) if j == 0 else 0
    space['cusp_G_dim'] = space['cusp_dim'] - space['cusp_Y_dim'] - space['cusp_P_dim']
    space['new_cusp_dim'] = paramodular_new_cusp_dim(k,j,N)
    # space['new_cusp_Y_dim'] = Yoshida_new_lift_dim(k,j,N)
    space['new_cusp_Y_dim'] = 0
    space['new_cusp_P_dim'] = Saito_Kurokawa_new_lift_dim(k,N) if j == 0 else 0
    space['new_cusp_G_dim'] = space['new_cusp_dim'] - space['new_cusp_Y_dim'] - space['new_cusp_P_dim']
    for aut_type in ['Y', 'P', 'G']:
        col = 'cusp_' + aut_type + '_dim'
        space['old_' + col] = space[col] - space['new_' + col]
    space['old_cusp_dim'] = space['cusp_dim'] - space['new_cusp_dim']
    return space

def smf_dims_paramodular_not_squarefree(k,j,N):
    space = {}
    space['degree'] =  2 
    space['family'] = 'K' 
    space['level'] = N 
    space['weight'] = [k, j]
    space['char_orbit_index'] = 1
    space['cusp_dim'] = 'NULL'
    space['cusp_Y_dim'] = 0
    space['cusp_P_dim'] = Saito_Kurokawa_lift_dim(k,N) if j == 0 else 0
    space['cusp_G_dim'] = 'NULL'
    space['new_cusp_dim'] = 'NULL'
    space['new_cusp_Y_dim'] = 0
    space['new_cusp_P_dim'] = Saito_Kurokawa_new_lift_dim(k,N) if j == 0 else 0
    space['new_cusp_G_dim'] = 'NULL'
    for aut_type in ['Y', 'P']:
        col = 'cusp_' + aut_type + '_dim'
        space['old_' + col] = space[col] - space['new_' + col]
    space['old_cusp_dim'] = 'NULL'
    return space

def smf_dims_paramodular(k,j,N):
    if is_squarefree(N):
        return smf_dims_paramodular_squarefree(k,j,N)
    return smf_dims_paramodular_not_squarefree(k,j,N)

def H1_prime(k,j,N):
    C = QQ(2)**(-7) * QQ(3)**(-3) * QQ(5)**(-1)
    D = (j+1)*(k-2)*(j+k-1)*(j+2*k-3)
    return C*D*prod([p**2 - 1 for p in prime_divisors(N)])

def H2_prime(k,j,N):
    C = -QQ(2)**(-7) * QQ(3)**(-2)
    D = (j+k-1)*(k-2)*(-1)**k
    E = 3 if (N == 2) else 0
    return C*D*E

def H3_prime(k,j,N):
    C = QQ(2)**(-5) * QQ(3)**(-1)
    E = 3 if (N == 2) else 0
    return C * C3(k,j) * E

def H4_prime(k,j,N):
    C = QQ(2)**(-3) * QQ(3)**(-3)
    E = 8 if (N == 3) else 0
    return C * C4(k,j) * E

def H5_prime(k,j,N):
    return 0

def H6_1_prime(N):
    a = 1/(QQ(2**5 * 3))
    b = 1/(QQ(2**7 * 3))
    prod_a = prod([p - IK_legendre(-1,p) for p in prime_divisors(N)])
    prod_b = prod([-1 + p*IK_legendre(-1,p) for p in prime_divisors(N)])
    return a*prod_a + b * prod_b

def H6_2_prime(N):
    a = (-1)**omega(N)/(QQ(2**5 * 3))
    b = 1/(QQ(2**7 * 3))
    prod_a = prod([p - IK_legendre(-1,p) for p in prime_divisors(N)])
    prod_b = prod([-1 + p*IK_legendre(-1,p) for p in prime_divisors(N)])
    return a*prod_a - b * prod_b

def H6_prime(k,j,N):
    a = (-1)**(j//2) * (2*k+j-3)
    b = (-1)**(j//2 + k) * (j+1)
    return a * H6_1_prime(N) + b * H6_2_prime(N)

def H7_1_prime(N):
    a = 1/(QQ(2**3 * 3**2))
    b = 1/(QQ(2**3 * 3**3))
    prod_a = prod([p - IK_legendre(-3,p) for p in prime_divisors(N)])
    prod_b = prod([-1 + p*IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a + b * prod_b

def H7_2_prime(N):
    a = (-1)**omega(N)/(QQ(2**3 * 3**2))
    b = 1/(QQ(2**3 * 3**3))
    prod_a = prod([p - IK_legendre(-3,p) for p in prime_divisors(N)])
    prod_b = prod([-1 + p*IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a - b * prod_b 

def H7_prime(k,j,N):
    return H7_1_prime(N)*mux([1,-1,0],j)*(2*k+j-3)+H7_2_prime(N)*mux([0,1,-1],j+2*k)*(j+1)

def H8_prime(k,j,N):
    return 0

def H9_prime(k,j,N):
    E = 3 if (N == 2) else 0
    return -QQ(2)**(-1)*QQ(3)**(-2) * C9(k,j) * E

def H10_prime(k,j,N):
    if len(NN(1,5,N) + NN(4,5,N)) == 0:
        F = 1 if N % 5 == 0 else 2
    else:
        F = 0
    E = (-1)**omega(N)*QQ(2)**(omega(N)-1)*F
    return QQ(5)**(-1) * C10(k,j) * E

def H11_prime(k,j,N):
    F = 1 if len(NN(1,8,N) + NN(7,8,N)) == 0 else 0
    delta = 1 if is_even(N) else 0
    E = (-1)**omega(N)*2**(omega(N)-delta)*F
    return QQ(2)**(-3) * C11(k,j) * E

def H12_1_prime(N):
    a = (-1)**omega(N)/(QQ(2**3 * 3))
    b = 1/(QQ(2**3 * 3))
    prod_a = prod([1 - IK_legendre(3,p) for p in prime_divisors(N)])
    prod_b = prod([IK_legendre(-1,p) - IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a - b * prod_b

def H12_2_prime(N):
    a = (-1)**omega(N)/(QQ(2**3 * 3))
    b = (-1)**omega(N)/(QQ(2**3 * 3))
    prod_a = prod([1 - IK_legendre(3,p) for p in prime_divisors(N)])
    prod_b = prod([IK_legendre(-1,p) - IK_legendre(-3,p) for p in prime_divisors(N)])
    return a*prod_a + b * prod_b 

def H12_prime(k,j,N):
    return H12_1_prime(N)*mux([1,-1,0],j)*(-1)**(j//2+k)+H12_2_prime(N)*mux([0,-1,1],j+2*k)*(-1)**(j//2)

def I1_prime(k,j,N):
    C = QQ(2)**(-3)*QQ(3)**(-1)
    D = j+1
    return C*D*prod([p-1 for p in prime_divisors(N)])

def I2_prime(k,j,N):
    return 0

def I3_prime(k,j,N):
    return 0

def I4_prime(k,j,N):
    return 0

def I5_prime(k,j,N):
    return 0

def I6_prime(k,j,N):
    return 0

def I7_prime(k,j,N):
    return 0

def I8_prime(k,j,N):
    return 0

def I9_prime(k,j,N):
    C = -QQ(2)**(-3)*(-1)**(j//2)
    return C*prod([-1 + IK_legendre(-1,p) for p in prime_divisors(N)])

def I10_prime(k,j,N):
    C = -QQ(2)**(-1)*QQ(3)**(-1)*mux([1,-1,0],j)
    return C*prod([-1 + IK_legendre(-3,p) for p in prime_divisors(N)])

H_prime = { n : eval('H'+str(n)+'_prime') for n in range(1,13)}

I_prime = { n : eval('I'+str(n)+'_prime') for n in range(1,11)}

def make_H_star(h, h_prime):
    def h_star(k,j,N):
        return sum([(-1)**omega(M)*2**omega(M)*h(k,j,N//M) for M in divisors(N)]) - h_prime(k,j,N)
    return h_star

def make_I_star(i, i_prime):
    def i_star(k,j,N):
        return sum([(-1)**omega(M)*2**omega(M)*i(k,j,N//M) for M in divisors(N)]) - (1 if is_even(omega(N)) else 0) * i_prime(k,j,N)
    return i_star

H_star = { i : make_H_star(H[i], H_prime[i]) for i in range(1,13)}

I_star = { i : make_I_star(I[i], I_prime[i]) for i in range(1,11)}

# when omega(N) is odd these should be given by the following (the rest are zero)

def H6_1_star(N):
    a = 1/(QQ(2**4 * 3))
    adm_divs = [M for M in divisors(N) if is_odd(omega(M))]
    s = 0
    for M in adm_divs:
        prod_M = prod([-1 + IK_legendre(-1,p) for p in prime_divisors(M)])
        prod_N_M = prod([p-1 for p in prime_divisors(N//M)])
        s += prod_M*prod_N_M
    return a*s

def H6_2_star(N):
    a = 1/(QQ(2**4 * 3))
    adm_divs = [M for M in divisors(N) if is_even(omega(M))]
    s = 0
    for M in adm_divs:
        prod_M = prod([-1 + IK_legendre(-1,p) for p in prime_divisors(M)])
        prod_N_M = prod([p-1 for p in prime_divisors(N//M)])
        s += prod_M*prod_N_M
    return a*s

def H6_star(k,j,N):
    a = (-1)**(j//2) * (2*k+j-3)
    b = (-1)**(j//2 + k) * (j+1)
    return a * H6_1_star(N) + b * H6_2_star(N)

def H7_1_star(N):
    a = 1/(QQ(2**2 * 3**2))
    adm_divs = [M for M in divisors(N) if is_odd(omega(M))]
    s = 0
    for M in adm_divs:
        prod_M = prod([-1 + IK_legendre(-3,p) for p in prime_divisors(M)])
        prod_N_M = prod([p-1 for p in prime_divisors(N//M)])
        s += prod_M*prod_N_M
    return a*s

def H7_2_star(N):
    a = 1/(QQ(2**2 * 3**2))
    adm_divs = [M for M in divisors(N) if is_even(omega(M))]
    s = 0
    for M in adm_divs:
        prod_M = prod([-1 + IK_legendre(-3,p) for p in prime_divisors(M)])
        prod_N_M = prod([p-1 for p in prime_divisors(N//M)])
        s += prod_M*prod_N_M
    return a*s

def H7_star(k,j,N):
    return H7_1_star(N)*mux([1,-1,0],j)*(2*k+j-3)+H7_2_star(N)*mux([0,1,-1],j+2*k)*(j+1)

def H12_1_star(N):
    a = -1/(QQ(2**2 * 3))
    adm_divs = [M for M in divisors(N) if is_even(omega(M))]
    s = 0
    for M in adm_divs:
        prod_M = prod([-1 + IK_legendre(-1,p) for p in prime_divisors(M)])
        prod_N_M = prod([-1 + IK_legendre(-3,p) for p in prime_divisors(N//M)])
        s += prod_M*prod_N_M
    return a*s

def H12_2_star(N):
    a = -1/(QQ(2**2 * 3))
    adm_divs = [M for M in divisors(N) if is_odd(omega(M))]
    s = 0
    for M in adm_divs:
        prod_M = prod([-1 + IK_legendre(-1,p) for p in prime_divisors(M)])
        prod_N_M = prod([-1 + IK_legendre(-3,p) for p in prime_divisors(N//M)])
        s += prod_M*prod_N_M
    return a*s

def H12_star(k,j,N):
    return H12_1_star(N)*mux([1,-1,0],j)*(-1)**(j//2+k)+H12_2_star(N)*mux([0,1,-1],j+2*k)*(-1)**(j//2)

def I1_star(k,j,N):
    C = QQ(2)**(-3)*QQ(3)**(-1)
    D = j+1
    return C*D*prod([p-1 for p in prime_divisors(N)])

def I3_star(k,j,N):
    return -QQ(2)**(omega(N)-5)*QQ(3)**(-2)*(j+1)*(2*k+j-3)*prod([p-1 for p in prime_divisors(N)])

def I4_star(k,j,N):
    return -QQ(2)**(omega(N)-5)*(-1)**k * prod([-1 + IK_legendre(-1,p) for p in prime_divisors(N)])

def I7_star(k,j,N):
    return -QQ(2)**(omega(N)-1)*QQ(3)**(-2)*mux([1,-1,0],j)*mux([0,1,-1],j+2*k)*prod([IK_legendre(-3,p)-1 for p in prime_divisors(N)])

def I9_star(k,j,N):
    C = -QQ(2)**(-3)*(-1)**(j//2)
    return C*prod([-1 + IK_legendre(-1,p) for p in prime_divisors(N)])

def I10_star(k,j,N):
    C = -QQ(2)**(-1)*QQ(3)**(-1)*mux([1,-1,0],j)
    return C*prod([-1 + IK_legendre(-3,p) for p in prime_divisors(N)])

def new_cmf_dim(N,k):
    d = (k-1) / 12 * prod([p-1 for p in prime_divisors(N)])
    e2 = (-1)**(k/2) / 4 * prod([IK_legendre(-1,p)-1 for p in prime_divisors(N)])
    e3 = mux([-1,0,1],k) / 3 *  prod([IK_legendre(-3,p)-1 for p in prime_divisors(N)])
    e_inf = (-1)**omega(N) if k == 2 else 0
    triv_level = 1/2 if N == 1 else 0
    return d+e2-e3+e_inf-triv_level
    
