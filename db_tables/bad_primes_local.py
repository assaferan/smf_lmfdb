"""
Local data at a prime p dividing the level N of a paramodular newform of weight (k, j),
in the conventions decided for the database (Roberts-Schmidt, as in omf5's
doc/euler-factors-bad-primes.md for the weight (3, 0), and omf5's README for the weights).

omf5's weight (a, b) is the paramodular weight (k, j) = (a - b + 3, 2b); write
w = 2a + 3 = 2k + j - 3 for the motivic weight.  omf5 normalizes T_p by p^a.  For (3, 0),
a = 0 and w = 3:

  p || N  (case 1): omf5 prints c = T_p + p^a, or -(T_p + p^a) when p divides the conductor d;
                    lambda = c - eps p^a, mu = -eps p^a c - p^(2a+1),
                    factor (1 - (c + eps p^(a+1)) x + p^w x^2) (1 + eps p^(a+1) x)
                         = 1 - c x + (p^w + p mu) x^2 + eps p^(w+a+1) x^3.
  p^2 | N (case 2, the radical mod p a line): lambda = T_p, mu = T_{p^2} / 2,
                    factor 1 - lambda x + (p^w + p mu) x^2.

The factors are those of omf5's README (checked there against the functional equation for the
weights (1, 1) and (2, 1)).  lambda and mu are the database's convention: mu is defined in both
cases by the coefficient p^w + p mu of x^2, and lambda = T_p when eps = 1 in case 1.  For a = 0
they are the values of Johnson-Leung and Roberts.

eps is the Atkin-Lehner sign at p.  For a Saito-Kurokawa lift (j = 0) of a newform f of weight
2k - 2 and level N, the factor is (1 - p^(k-2) x)(1 - p^(k-1) x)(1 - a_p(f) x), and lambda, mu
are read off it by the same dictionary.  All values may be elements of a number field.
"""
from sage.all import PolynomialRing, QQ


def local_from_omf5(case, p, eps, c=None, Tp=None, Tp2=None, a=0):
    """(lambda, mu, factor coefficients [1, b1, b2, b3]) from omf5's printed values, for omf5's
    weight (a, b) (only a enters)."""
    w = 2 * a + 3
    if case == 1:
        lam = c - eps * p**a
        # (1 - (c + eps p^(a+1)) x + p^w x^2)(1 + eps p^(a+1) x)
        b1 = -(c + eps * p**(a + 1)) + eps * p**(a + 1)
        b2 = p**w - (c + eps * p**(a + 1)) * eps * p**(a + 1)
        b3 = p**w * eps * p**(a + 1)
        mu = (b2 - p**w) / p
        assert mu == -eps * p**a * c - p**(2 * a + 1)
        return lam, mu, [1, b1, b2, b3]
    if case == 2:
        lam, mu = Tp, Tp2 / 2
        return lam, mu, [1, -lam, p**w + p * mu, 0]
    raise ValueError(case)


def local_from_sk(case, p, eps, ap_f, k=3):
    """(lambda, mu, factor) of the SK lift of weight (k, 0) of f, ap_f = a_p(f) (0 when p^2
    divides the level of f)."""
    a, w = k - 3, 2 * k - 3
    # factor: (1 - p^(k-2) x)(1 - p^(k-1) x)(1 - ap_f x)
    R = PolynomialRing(ap_f.parent() if hasattr(ap_f, 'parent') else QQ, 'x')
    x = R.gen()
    F = (1 - p**(k - 2) * x) * (1 - p**(k - 1) * x) * (1 - ap_f * x)
    co = [F[i] for i in range(4)]
    if case == 1:
        # match (1 - (c + eps p^(a+1)) x + p^w x^2)(1 + eps p^(a+1) x): the coefficient of x is -c
        c = -co[1]
        lam, mu, fac = local_from_omf5(1, p, eps, c=c, a=a)
    else:
        lam = -co[1]
        mu = (co[2] - p**w) / p
        fac = [1, -lam, p**w + p * mu, 0]
    assert [fac[i] - co[i] for i in range(4)] == [0, 0, 0, 0], (case, p, eps, ap_f, fac, co)
    return lam, mu, fac


def prime_power_coefficients(fac, kmax):
    """a_{p^k}, k = 0..kmax: the coefficients of 1 / (1 + b1 x + b2 x^2 + b3 x^3)."""
    a = [fac[0] ** 0]
    for k in range(1, kmax + 1):
        s = 0
        for i in range(1, min(k, 3) + 1):
            s -= fac[i] * a[k - i]
        a.append(s)
    return a
