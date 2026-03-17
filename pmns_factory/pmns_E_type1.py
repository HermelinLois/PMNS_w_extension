from sage.all import PolynomialRing, ZZ, ceil, Integer, GF
from pmns_core.parameters.params_gestion import search_minimal_degree as SMD, search_base_rho_and_gamma, search_memory_overhead, cast_polynomial_to_minimal_representation
from pmns_core.parameters.roots_gestion import is_gamma_feasible, search_roots

PR = PolynomialRing(ZZ, "X")
X = PR("X")
INIT_ALPHA = 1
INIT_BETA = 2

def gen_pol_e(n, k, alpha, beta):
    return X**n - alpha * X**k - beta

def increase_parameters(pol_e, p:int, k:int, phi:int) -> tuple:
    n = pol_e.degree()
    beta = - pol_e[0]
    alpha = -pol_e[k]
    
    n_alpha = alpha + 1
    n_beta = beta + 1

    E_alpha = gen_pol_e(n, k, n_alpha, beta)
    E_beta = gen_pol_e(n, k, alpha, n_beta)
    
    w_alpha = search_memory_overhead(E_alpha)
    w_beta = search_memory_overhead(E_beta)

    p_pow_kn = p**(k/n)
    test_alpha = (2 * w_alpha * p_pow_kn >= phi)
    test_beta = (2 * w_beta * p_pow_kn >= phi)

    if test_alpha and test_beta:
        return INIT_ALPHA, INIT_BETA, n + k
    
    if test_beta:
        return n_alpha, INIT_BETA, n 
    
    return alpha, n_beta, n


def search_minimal_degree(p: int, k: int, phi_pow: int) -> int:
    max_add_coef = INIT_BETA + INIT_ALPHA*(INIT_BETA + INIT_ALPHA)
    n = SMD(p, k, phi_pow, max_add_coef)

    return int(k * ceil(n/k))


def gen_parameters(p:int, k:int, phi_pow:int=64, name:str ="z") -> dict:
    assert is_gamma_feasible(p, k), "no gamma satisfy the construction"
    assert k > 1, f"extension degree must be at least 2, here {k=}"

    p = Integer(p)
    assert p.nbits() >= phi_pow, f"construction only works if the number of bits in prime (here {p=}) is greater or equal to {phi_pow=}"

    n = search_minimal_degree(p, k, phi_pow)
    alpha = INIT_ALPHA
    beta = INIT_BETA
    phi = 2**phi_pow

    K = GF(p**k, name)

    parameters_not_found = True
    result = None
    while parameters_not_found:
        pol_e = gen_pol_e(n, k, alpha, beta)
        roots = search_roots(p, k, pol_e, K)

        if roots:
            result = search_base_rho_and_gamma(roots, k, p, phi, pol_e)
            parameters_not_found = (result is None)

        if parameters_not_found:
            alpha, beta, n = increase_parameters(pol_e, p, k, phi)
    
    L, rho, gamma = result
    return {'rho': rho, 'gamma': gamma, 'phi_pow': phi_pow, 'L': L, 'E': pol_e, 'mod': cast_polynomial_to_minimal_representation(K.modulus(), p), 'p': p, 'k':k}