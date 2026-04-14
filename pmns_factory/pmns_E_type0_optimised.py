# ==================================================
# Objectif is to construct E = X^n - lambda
# with fast search of root by constructing a 
# specific structure of the extension field
# ==================================================

from sage.all import PolynomialRing, ZZ, vector, ceil, Integer, GF, random_prime, factor, randint, gcd, matrix
from core.parameters.params_gestion import search_minimal_degree as SMD, search_base_rho_and_gamma, search_memory_overhead, cast_polynomial_to_minimal_representation
from core.parameters.roots_gestion import is_gamma_feasible, search_roots

# Cantor-Zassenhaus
# Berlekamp
# those algorithms can be used to assertthat polynomials are irreducible

PR = PolynomialRing(ZZ, "X")
X = PR.gen()
INIT_LAMB = 2

def gen_pol_e(n:int, lamb:int):
    return X**n - lamb

def construct_irreducible_polynomial(k, p, phi_pow):
    # check if k = 0 mod(4) and in that case check 
    # if p = 1 mod(4) 
    if k%4 == 0 and p%4 != 1:
        return None

    primes_k = [q for q, _ in factor(k)]

    # fast verification to see if it's impossible to 
    # construct an irreducible polynomial with p and k
    for q in primes_k:
        if (p - 1) % q != 0:
            return None

    # search any element wich isn't a q power of k where q is all prime factor of k
    # generators of Z/pZ are elements wich satisfy those condition 
    exponents = [(p - 1) // q for q in primes_k]
    limit = min(2**(phi_pow), p)
    for epsilon in range(2, limit):
        if all(pow(epsilon, exp, p) != 1 for exp in exponents):
            return X**k - epsilon


def increase_parameters(pol_e, p:int, k:int, phi:int) -> tuple:
    n = pol_e.degree()
    lamb = - pol_e[0]
    n_lamb = lamb + 1

    E = gen_pol_e(n, n_lamb)

    # compute overhead based on: "PMNS for efficient arithmetic and small memory cost" 
    # (J. Robert, P. Véron, F. Dosso, 2022)
    w = search_memory_overhead(E)
    if 2 * w * p**(k/n) >= phi:
        return INIT_LAMB, n + 1
    return n_lamb, n


def search_minimal_degree(p: int, k: int, phi_pow: int) -> int:
    #initialy, max add coef is lamb for external reduction
    # so we serach a minimal n with this lamb
    max_add_coef = lambda n : search_memory_overhead( gen_pol_e(n, INIT_LAMB) )
    n = SMD(p, k, phi_pow, max_add_coef)
    
    # as we search root of polynomial pol_e in extension field, wich ar e intergers of power k, 
    # we have to let the degree n be a multiple of k
    return n


def gen_parameters(p:int, k:int, phi_pow:int=64, name:str="z") -> dict:    
    # first condition ensure there existe element in space suitable for our PMNS
    # second permit to know if we can construct irreducible polynomial
    assert gcd(k, p-1) > 1, "impossible to construct an irreducible polynomial over Z/pZ"

    p = Integer(p)
    assert p.nbits() >= phi_pow, f"construction only works if the number of bits in prime (here {p=}) is greater or equal to {phi_pow=}"

    n = search_minimal_degree(p, k, phi_pow)
    lamb = INIT_LAMB
    phi = 2**phi_pow

    # extension field creation
    # if mod is None, sage will use an irreducible polynomial
    # of a generic form
    mod = construct_irreducible_polynomial(k, p, phi_pow)
    K = GF(p**k, name=name, modulus=mod)

    parameters_not_found = True
    result = None
    iteration = 0
    while parameters_not_found:
        iteration +=1
        # construction of the polynomial pol_e
        pol_e = gen_pol_e(n, lamb)
        roots = search_roots(p, k, pol_e, K)

        if roots:
            # search suitable element to construct a PMNS using found roots
            result = search_base_rho_and_gamma(roots, k, p, phi, pol_e)
            parameters_not_found = (result is None)

        if parameters_not_found:
            lamb, n = increase_parameters(pol_e, p, k, phi)

    L, rho, gamma = result
    return {'rho': rho, 'gamma': gamma, 'phi_pow': phi_pow, 'L': L, 'E': pol_e, 'mod': cast_polynomial_to_minimal_representation(K.modulus(), p), 'p': p, 'k':k, 'it':iteration}
