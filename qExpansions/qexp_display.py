from sage.all import *
# from sage.all import (LaurentPolynomialRing, PolynomialRing, Rationals)
cwd = os.getcwd()
os.chdir('smf_lmfdb/qExpansions')
load('FE_FG_20_0.sage')
os.chdir(cwd)

def make_qexp_display(f, nterms=5):
    '''
    creates a display for the first few terms in the q-expansion from f
    assumes f is given as multivariate Laurent series in 3 variables - q_1, q_2 and q_{12} (in this order)
    '''
    d = f.dict()
    max_d = max([k[0] + k[1] for k in d.keys()])
    key_by_deg = [[k for k in d.keys() if k[0]+k[1] == deg] for deg in range(max_d)]
    s = 0
    i = 0
    while (s < nterms):
        s += len(key_by_deg[i])
        i += 1
    Lq = LaurentPolynomialRing(Rationals(), names=['q12'])
    q12 = Lq.gens()[0]
    PLq = PolynomialRing(Lq, 2, names=['q1','q2'])
    q1, q2 = PLq.gens() 
    terms = []
    for j in range(i):
        exps = sorted(list({(k[0], k[1]) for k in key_by_deg[j]}))
        deg_hom = [sum([d[k]*q12**k[2] for k in key_by_deg[j] if (k[0], k[1]) == exp])*q1**exp[0]*q2**exp[1] for exp in exps]        
        terms += deg_hom
    str_terms = [str(t).replace('q12', 'q_{12}').replace('q1', 'q_{1}').replace('q2', 'q_{2}').replace('*', '') for t in terms]
    return '+'.join(str_terms + ['\\cdots'])

def get_qexp_display_F20G():
    return make_qexp_display(F20G)
