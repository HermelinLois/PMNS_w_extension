# ==================================================
# matrix_generation.py
# Generic functions around matrix used for PMNS construction
# and C implementation
# ==================================================

from sage.all import matrix, ZZ, PolynomialRing
from math_utils import square_and_multiply

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


def gen_overflow_matrix(n: int, pol_e):
    """
    Create a matrix used to compute the coefficients after reduction modulo pol_e 
    for powers of X greater or equal than the degree of polynomial pol_e.

    Args:
        n (int): degree of pol_e
        pol_e (Polynomial): polynomial used for external reduction in PMNS

    Returns:
        matrix (ZZ): matrix representing the reduction of X^(n+i) modulo pol_e
    """
    PR = PolynomialRing(ZZ, "X")
    X = PR("X")

    matrix_coefficients = []
    for i in range(n - 1):
        poly_mod = square_and_multiply(X, n + i, pol_e)
        # Pad coefficients to length n
        coeffs = list(poly_mod) + [0] * (n - len(list(poly_mod)))
        matrix_coefficients.append(coeffs)
    
    return matrix(ZZ, matrix_coefficients)


def gen_external_reduction_matrix(pol_m, pol_e, n: int, phi: int):
    """
    Generate matrices pol_m and N for external reduction in PMNS.

    Args:
        pol_m (Polynomial): polynomial null in the chosen root of pol_e
        pol_e (Polynomial): polynomial used for external reduction
        n (int): degree of pol_e
        phi (int): word size in bits (used for modulo)

    Returns:
        mat_m (matrix): matrix representing pol_m for Montgomery reduction
        mat_n (matrix): matrix representing N = -pol_m^(-1) modulo phi
    """
    PR = PolynomialRing(ZZ, "X")
    X = PR("X")

    matrix_coefficients = []
    for i in range(n):
        # Compute (pol_m * X^i) % pol_e
        poly_mod = (pol_m * X**i) % pol_e
        coeffs = list(poly_mod) + [0] * (n - len(list(poly_mod)))
        matrix_coefficients.append(coeffs)

    mat_m = matrix(ZZ, matrix_coefficients)
    mat_n = -mat_m.inverse() % phi

    return mat_m, mat_n

