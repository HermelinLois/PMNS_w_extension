# ==================================================
# montgomery_reduction.py
# functions wich permit to apply Montgomery reduction 
# to a polynomial 
# ==================================================

from sage.all import matrix, ZZ, PolynomialRing, Integer, vector
from ...math_utils import square_and_multiply

PR = PolynomialRing(ZZ, "X")

def gen_mn_reduction_matrix(M, E, phi: int):
    """
    Generate matrices M and N for external reduction in PMNS.

    Args:
        M (Polynomial): polynomial null in the chosen root of E
        E (Polynomial): polynomial used for external reduction
        phi (int): word size in bits (used for modulo)

    Returns:
        mat_m (matrix): matrix representing M for Montgomery reduction
        mat_n (matrix): matrix representing N = -M^(-1) modulo phi
    """
    n = E.degree()
    X = PR("X")

    matrix_coefficients = []
    for i in range(n):
        # Compute (M * X^i) % E
        poly_mod = (M * square_and_multiply(X,i)) % E
        coeffs = list(poly_mod) + [0] * (n - len(list(poly_mod)))
        matrix_coefficients.append(coeffs)

    mat_m = matrix(ZZ, matrix_coefficients)
    mat_n = -mat_m.inverse() % phi

    return mat_m, mat_n


def montgomery_reduction(pol_p, M, N, E, gamma, phi:int =2**64):
    """
    Reduction of a polynomial with Montgomery reduction

    Args:
        pol_p (Polynomial | matrix): elment wich need to be reduced 
        M (Polynomial | matrix): element wich represent M such that M(gamma)=0
        N (Polynomial | matrix): element wich represent N = -M^-1
        E (Polynomial): polynomial for external reduction
        phi (int, opt): represent bit word size use by the architecture. usualy a power of 2 and is 2^64 by default
        gamma (extension field element): root of E

    Returns:
        Polynomial : reduction of pol_p which still represent the same element
    """
    phi = Integer(phi)

    Q = ((pol_p * N) % E) % phi
    Q = PR([elmt - phi * (elmt>2**(phi.nbits() - 2)) for elmt in Q])

    T = (Q * M) % E 
    reduction = (pol_p + T) // phi

    assert reduction(gamma) * phi == pol_p(gamma), f"Error occurring during reduction. Please check parameters"
    
    return reduction


def fast_montgomery_reduction(pol_p, E, L, L_inv, phi:int = 2**64):
    n = L.nrows()
    phi = Integer(phi)

    v_p = vector(ZZ, pol_p.list() + [0] * (n - len(pol_p.list())))
    q_raw = (v_p * L_inv) % phi
    q_centered = vector(ZZ, [(c % phi) - (phi if (c % phi) > phi // 2 else 0) for c in q_raw])
    
    T = PR((q_centered * L).list())
        
    return (pol_p + T)//phi