# ==================================================
# matrix_gestion.py
# Generic functions around matrix used for PMNS construction
# and C implementation
# ==================================================

from sage.all import matrix, ZZ, PolynomialRing
from pmns_core.math_utils import square_and_multiply

PR = PolynomialRing(ZZ, "X")

def gen_reduce_null_base(k:int, p:int, n:int, gamma):
    """
    Generate a reduced base of null polynomial for gamma.
    Approach based on: An Alternative Approach for SIDH Arithmetic (C.Bouvier & L.Imbert)
    
    Args:
        p (int): prime used to construct the extension field
        k (int): extension degree of the field
        n (int): degree of pol_e used for polynomial reduction in PMNS
            gamma: element used to construct PMNS (gamma^k is integer)
    
    Returns:
        matrix (LLL) : reduced base where each row evaluates to zero at gamma
    """

    gamma_pow_k = -square_and_multiply(gamma, k)
    base = matrix(ZZ, n, n, 0)

    # fill the diagonal
    for i in range(k):
        base[i, i] = p
        
    for i in range(k, n):
        base[i, i] = 1
        base[i, i-k] = gamma_pow_k

    # LLL reduction
    return base.LLL()


def gen_overflow_matrix(pol_e):
    """
    Create a matrix used to compute the coefficients after reduction modulo pol_e 
    for powers of X greater or equal than the degree of polynomial pol_e.

    Args:
        pol_e (Polynomial): polynomial used for external reduction in PMNS

    Returns:
        matrix (ZZ): matrix representing the reduction of X^(n+i) modulo pol_e
    """
    n = pol_e.degree()
    X = PR("X")

    matrix_coefficients = []
    for i in range(n-1):
        poly_mod = X**(n + i) % pol_e
        # Pad coefficients to length n
        coeffs = list(poly_mod) + [0] * (n - poly_mod.degree()-1)
        matrix_coefficients.append(coeffs)
    
    return matrix(ZZ, matrix_coefficients)