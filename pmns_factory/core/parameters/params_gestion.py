# ==================================================
# params_generation.py
# Generic functions around parameters generation for
# PMNS construction
# ==================================================

from sage.all import vector, infinity, ZZ, PolynomialRing, ceil, exp, Integer
from core.parameters.matrix_gestion import gen_overflow_matrix, gen_reduce_null_base
from core.operations.reductions.montgomery_reduction import search_m_with_even_degs, search_m_with_odd_deg, search_polynomial_m, search_m_and_n

PR = PolynomialRing(ZZ, "X")

def search_memory_overhead(pol_e) -> int:
    """
    Compute the memory overhead coefficient w used in PMNS bounds.

    Based on:
    "PMNS for efficient arithmetic and small memory cost"
    (J. Robert, P. Véron, F. Dosso, 2022)

    The coefficient bound is given by:

        w = || v1 + v2 * e ||∞

    where:
        e  : overflow matrix associated with polynomial pol_e
        v1 : vector [1,2,...,n]
        v2 : vector [n-1, n-2, ..., 1]

    Args:
        pol_e (polynomial): polynomial used for external reduction in PMNS

    Returns:
        int: approximation of the coefficient growth factor during
             polynomial multiplication in PMNS
    """

    n = pol_e.degree()

    # Overflow matrix
    epsilon = gen_overflow_matrix(pol_e)

    # Construct vectors
    v1 = vector(range(1, n + 1))
    v2 = vector(range(n - 1, 0, -1))

    # Compute bound
    return (v1 + v2 * epsilon).norm(infinity)



def search_base_rho_and_gamma(roots: list, k: int, p: int, phi: int, pol_e):
    """
    Search a suitable PMNS base, rho bound, and gamma root.

    Args:
        roots (list): roots of pol_e in the extension field
        k (int): extension degree
        p (int): prime used to construct the extension field
        phi (int): word size bound
        pol_e (Polynomial): polynomial used for external reduction in PMNS

    Returns:
        tuple:
            base (matrix): LLL-reduced base of null polynomials at gamma
            rho (int): coefficient bound for PMNS representation
            gamma: root of pol_e suitable for PMNS construction

        None if no suitable root is found.
    """

    n = pol_e.degree()
    # coefficient growth factor
    w = search_memory_overhead(pol_e)

    for gamma in roots:
        # Generate reduced lattice base
        base = gen_reduce_null_base(k, p, n, gamma)

        # Norm of the base
        rho = Integer(base.norm(1) - 1)
        
        # Ensure coefficients fit machine word
        if 2 * w * (rho - 1) < phi:
            return base, rho, gamma

    return None



def search_minimal_degree(p: int, k: int, phi_pow: int, init_polynomial: callable) -> int:
    """
    Function that compute a minimal value of n such that we can possibly construct a PMNS

    Args:
        p (Interger): prime used to construction extension field
        k (int): extension degree
        phi_pow (int): word size
        init_polynomial (callable): accepts only polynomial degree and return polynomial with initial 
            coeffients parameters given to construct polynomial

    Returns:
        int: return a degree n minimal for wich we can possibly construct a PMNS
    """
    pbits = p.nbits()
    # knowing that n must verify that (2rho -1)^n > p^k
    # n must be at least represent p over log2(phi) register.
    # So n must be greater that the bit size of p^k, ie, #(p)2 * k 
    n = int(pbits * k / phi_pow) + 1
    phi = 2**phi_pow

    # compute minimal degree n wich can lead to a possible contruction of PMNS
    # here we approximate a value of the laticce G such that rho >= ||G||-1    
    n = int(pbits * k / phi_pow) + 1
    while round( 2 * search_memory_overhead(init_polynomial(n)) * (max(gen_overflow_matrix(init_polynomial(n))._list())/2 * 2**(k * pbits / n) * ceil(n / exp(1)) - 2)) >= phi:
        n += 1
    return n


def cast_polynomial_to_minimal_representation(pol, p):
    """
    Retrun polynomial pol with coefficients in signed samll representation with modulus p.

    Args:
        p (Interger): prime used to construction extension field
        pol(Polynomial): polynomial that need to have reduced coefficients

    Returns:
        Polynomial: return polynomial with coefficients in [-p//2, p//2]
        
    Exemple : (p-1)*X => -X[p]
    """
    coefficients = [int(c) if c <= p//2 else int(c)-p for c in pol]
    return PR(coefficients)