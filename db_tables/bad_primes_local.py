"""
Local data at a prime p dividing the level N of a paramodular newform of weight (3, 0),
in the conventions decided for the database (Roberts-Schmidt, as in omf5's
doc/euler-factors-bad-primes.md):

  p || N  (case 1): omf5 prints c = T_p + 1, or -(T_p + 1) when p divides the conductor d;
                    lambda = c - eps, mu = -eps c - p,
                    factor (1 - (c + eps p) x + p^3 x^2) (1 + eps p x).
  p^2 | N (case 2, the radical mod p a line): lambda = T_p, mu = T_{p^2} / 2,
                    factor 1 - lambda x + (p^3 + p mu) x^2.

eps is the Atkin-Lehner sign at p.  For a Saito-Kurokawa lift of a newform f of weight 4
and level N, the factor is (1 - p x)(1 - p^2 x)(1 - a_p(f) x), and lambda, mu are read off
it by the same dictionary.  All values may be elements of a number field.
"""
from sage.all import PolynomialRing, QQ


def local_from_omf5(case, p, eps, c=None, Tp=None, Tp2=None):
    """(lambda, mu, factor coefficients [1, b1, b2, b3]) from omf5's printed values."""
    if case == 1:
        lam = c - eps
        mu = -eps * c - p
        # (1 - (c + eps p) x + p^3 x^2)(1 + eps p x)
        b1 = -(c + eps * p) + eps * p
        b2 = p**3 - (c + eps * p) * eps * p
        b3 = p**3 * eps * p
        return lam, mu, [1, b1, b2, b3]
    if case == 2:
        lam, mu = Tp, Tp2 / 2
        return lam, mu, [1, -lam, p**3 + p * mu, 0]
    raise ValueError(case)


def local_from_sk(case, p, eps, ap_f):
    """(lambda, mu, factor) of the SK lift of f, ap_f = a_p(f) (0 when p^2 divides the level of f)."""
    # factor: (1 - p x)(1 - p^2 x)(1 - ap_f x)
    R = PolynomialRing(ap_f.parent() if hasattr(ap_f, 'parent') else QQ, 'x')
    x = R.gen()
    F = (1 - p * x) * (1 - p**2 * x) * (1 - ap_f * x)
    co = [F[i] for i in range(4)]
    if case == 1:
        # match (1 - (c + eps p) x + p^3 x^2)(1 + eps p x): the coefficient of x is -c
        c = -co[1]
        lam, mu, fac = local_from_omf5(1, p, eps, c=c)
    else:
        lam = -co[1]
        mu = (co[2] - p**3) / p
        fac = [1, -lam, p**3 + p * mu, 0]
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
