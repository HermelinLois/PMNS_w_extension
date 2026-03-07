from sage.all import PolynomialRing, ZZ, ceil, Integer, exp, random_prime
from pmns_core.parameters.params_gestion import search_base_rho_and_gamma, search_memory_overhead
from pmns_core.parameters.roots_gestion import is_gamma_feasible, search_roots

PR = PolynomialRing(ZZ, "X")

def gen_pol_e(n:int, coefs:list):
    X = PR("X")
    E = X**n
    for coef, deg in coefs:
        E -= coef * X**deg
    return E


def increase_parameters(pol_e, p, k, coefs, _init_coefs, increase_n, phi):
    n = pol_e.degree()
    for coef, deg in coefs:
        temp_coefs = coefs
        current_coef = coef + 1
        
        temp_coefs[deg] = current_coef, temp_coefs[deg][1]
        pol_e = gen_pol_e(n, temp_coefs)
        w = search_memory_overhead(pol_e)
        
        if 2 * w * p**(k/n) < phi:
            return n, temp_coefs
    return increase_n(n,k), _init_coefs


def compute_max_add_coef(n, k, _init_coefs):
    _init_coefs = sorted(_init_coefs, key=lambda x: x[1])
    
    def compute_coef(u):
        if u == 0:
            return 0
        return sum(abs(c[0]) for c in _init_coefs[:-1]) + abs(_init_coefs[-1][0])*compute_coef(u-1)
    
    u = 1
    while k > u*(n-k)+2:
        u += 1
    return compute_coef(u)

def search_minimal_degree(p: int, k: int, return_n, _init_coefs, phi_pow: int):
    pbits = p.nbits()
    
    inner_degree = max(_init_coefs, key=lambda x: x[1])[1]
    n = max(int(pbits * k / phi_pow), inner_degree) + 1
    phi = 2**phi_pow

    while round(2**(k * pbits / n) * ceil(n / exp(1)) * 2 * (compute_max_add_coef(n, k, _init_coefs) * (n - 1) + 1)) >= phi:
        n += 1
        
    return return_n(n, k)


def gen_parameters(p:int, k:int, _init_coefs, phi_pow:int=64, increase_n=lambda x,k: x + 1, return_n=lambda x,k: x):
    
    assert is_gamma_feasible(p, k), "no gamma possibly satisfy the construction"
    assert k > 1, f"extension degree must be at least 2, here {k=}"

    p = Integer(p)
    assert p.nbits() >= phi_pow, f"construction only works if the number of bits in prime (here {p=}) is greater or equal to {phi_pow=}"

    n = search_minimal_degree(p, k, return_n, _init_coefs, phi_pow)
    coefs = _init_coefs
    phi = 2**phi_pow

    parameters_not_found = True
    result = None
    while parameters_not_found:
        pol_e = gen_pol_e(n, coefs)
        print(pol_e)
        roots = search_roots(p, k, pol_e)

        if roots:
            result = search_base_rho_and_gamma(roots, k, p, phi, pol_e)
            parameters_not_found = (result is None)

        if parameters_not_found:
            n, coefs = increase_parameters(pol_e, p, k, coefs, _init_coefs, increase_n, phi)
    
    L, rho, gamma = result
    return {'rho': rho, 'gamma': gamma, 'phi_pow': phi_pow, 'L': L, 'E': pol_e, 'mod': PR(gamma.parent().modulus()), 'p': p, 'k':k}