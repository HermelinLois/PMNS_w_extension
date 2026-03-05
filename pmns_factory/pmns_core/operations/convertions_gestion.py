# ==================================================
# conversions_gestion.py
# function wich permit to convert element to polynomial
# or/and pmns representation using Montgomery or Babai
# reduction
# ==================================================

from sage.all import vector, matrix, PolynomialRing, ZZ
from reductions.montgomery_reduction import montgomery_reduction

def gen_gamma_base(gamma, k: int) -> matrix:
    """
    Generate base of extension field using powers of gamma.

    Args:
        gamma (extension field element): root of the external reduction
        k (int): degree of the extension

    Returns:
        matrix: each column is gamma^i for i=0..k-1
    """
    gamma_base = [gamma**i for i in range(k)]
    return matrix([g._vector_() for g in gamma_base]).transpose()


def convert_element_to_polynomial(element, gamma_base: matrix):
    """
    Represent an extension field element as a polynomial in gamma.

    Args:
        element (extension field element): element to represent
        gamma_base (matrix): matrix of powers of gamma (columns = powers)

    Returns:
        Polynomial: polynomial representing the element
    """
    target_values = vector(element._vector_())
    gamma_decomposition = gamma_base.solve_right(target_values)

    PR = PolynomialRing(ZZ, "X")
    polynomial_of_element = PR(list(gamma_decomposition))

    # gamma_base[1] corresponds to gamma^1
    gamma = gamma_base.column(1)
    assert polynomial_of_element(gamma) == element, \
        f"error in the construction, polynomial doesn't represent {element=}"

    return polynomial_of_element


def convert_element_to_pmns_montgomery(element, gamma_base, **kwargs):
    """
    Convert an extension field element to PMNS using Montgomery reduction.

    Args:
        element (extension field element): element to convert
        gamma_base (matrix): powers of chosen gamma
        kwargs: must include phi, rho, gamma, M, N, E

    Returns:
        Polynomial: PMNS representation
    """
    phi, rho, gamma = kwargs['phi'], kwargs['rho'], kwargs['gamma']
    M, N, E = kwargs['M'], kwargs['N'], kwargs['E']
    n = E.degree()

    alpha = element * phi**n
    V = convert_element_to_polynomial(alpha, gamma_base)

    for _ in range(n):
        V = montgomery_reduction(V, M, N, E, phi)

    assert V(gamma) == element, f"polynomial doesn't represent {element=}"
    assert all(abs(c) < rho for c in V), f"{rho=} too low for {element=}"

    return V


#TODO implement babai convertion to pmns
def convert_element_to_pmns_babai(element, **kwargs):
    raise Warning("method not implemented yet")