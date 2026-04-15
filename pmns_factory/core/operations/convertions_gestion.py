# ==================================================
# conversions_gestion.py
# function wich permit to convert element to polynomial
# or/and pmns representation using Montgomery or Babai
# reduction
# ==================================================

from sage.all import vector, matrix, PolynomialRing, ZZ, log
from core.operations.reductions.montgomery_reduction import fast_montgomery_reduction

PR = PolynomialRing(ZZ, "X")

def gen_transition_matrix(gamma, k: int) -> matrix:
    """
    Generate transposition matrix to express elements in gamma basis

    Args:
        gamma (extension field element): root of the external reduction
        k (int): degree of the extension

    Returns:
        matrix: transposition matrix from canonical basis to gamma basis
    """
    mat = matrix([(gamma**i)._vector_() for i in range(k)])
    return matrix(ZZ, mat.inverse())


def convert_element_to_polynomial(element, gamma, transition_matrix: matrix):
    """
    Represent an extension field element as a polynomial in gamma.

    Args:
        element (extension field element): element to represent
        transition_matrix (matrix): matrix from canonical to gamma basis

    Returns:
        Polynomial: polynomial representing the element
    """
    gamma_decomposition = element._vector_() * transition_matrix

    polynomial_of_element = PR(list(gamma_decomposition))
    
    assert polynomial_of_element(gamma) == element, f"error in the construction, polynomial doesn't represent {element=}"

    return polynomial_of_element


def convert_element_to_pmns_montgomery(element, transition_matrix, pmns):
    """
    Convert an extension field element to PMNS using Montgomery reduction.

    Args:
        element (extension field element): element to convert
        transition_matrix (matrix):  matrix from canonical to gamma basis
        pmns: must include phi_pow, rho, gamma, M, N, E

    Returns:
        Polynomial: PMNS representation
    """
    phi_pow, rho, gamma = pmns['phi_pow'], pmns['rho'], pmns['gamma']
    E = pmns['E']
    L, L_inv = pmns['L'], pmns['L_inv']
    k = pmns['k']

    # retrieve parameters from given elements
    n = E.degree()
    phi = 2**phi_pow
    nb_iteration = n

    alpha = element * phi**nb_iteration
    V = convert_element_to_polynomial(alpha, gamma, transition_matrix)   
    
    for i in range(nb_iteration):
        V = fast_montgomery_reduction(V, E, L, L_inv, phi)

    assert V(gamma) == element, f"polynomial doesn't represent {element=}\n{V(gamma)}"
    assert all(abs(c) < rho for c in V), f"{rho=} too low for {element=}\n{V =}"

    return V