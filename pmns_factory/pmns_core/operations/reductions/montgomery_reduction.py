# ==================================================
# montgomery_reduction.py
# functions wich permit to apply Montgomery reduction 
# to a polynomial 
# ==================================================

from sage.all import matrix, ZZ, PolynomialRing
from ...math_utils import square_and_multiply

PR = PolynomialRing(ZZ, "X")

def gen_external_reduction_matrix(pol_m, pol_e, phi: int):
    """
    Generate matrices pol_m and N for external reduction in PMNS.

    Args:
        pol_m (Polynomial): polynomial null in the chosen root of pol_e
        pol_e (Polynomial): polynomial used for external reduction
        phi (int): word size in bits (used for modulo)

    Returns:
        mat_m (matrix): matrix representing pol_m for Montgomery reduction
        mat_n (matrix): matrix representing N = -pol_m^(-1) modulo phi
    """
    n = pol_e.degree()
    X = PR("X")

    matrix_coefficients = []
    for i in range(n):
        # Compute (pol_m * X^i) % pol_e
        poly_mod = (pol_m * square_and_multiply(X,i)) % pol_e
        coeffs = list(poly_mod) + [0] * (n - len(list(poly_mod)))
        matrix_coefficients.append(coeffs)

    mat_m = matrix(ZZ, matrix_coefficients)
    mat_n = -mat_m.inverse() % phi

    return mat_m, mat_n


def montgomery_reduction(pol_p, pol_m, pol_n, pol_e, phi:int):
    """
    Reduction of a polynomial with Montgomery reduction

    Args:
        pol_p (Polynomial | matrix): elment wich need to be reduced 
        pol_m (Polynomial | matrix): element wich represent M such that M(gamma)=0
        pol_n (Polynomial | matrix): element wich represent N = -M^-1
        pol_e (Polynomial): polynomial for external reduction
        phi (int): represent bit word size use by the architecture

    Returns:
        Polynomial : reduction of pol_p which still represent the same element
    """
    Q = ((pol_p * pol_n) % pol_e) % phi
    T = (Q * pol_m) % pol_e 
    return (pol_p + T) // phi