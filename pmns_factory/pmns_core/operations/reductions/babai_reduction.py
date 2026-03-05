# ==================================================
# babai_reduction.py
# Functions to apply Babai lattice reduction for PMNS.
# Provides three variants:
#   > nearest plane reduction
#   > rounding reduction with unlimited precision
#   > rounding reduction with limited precision
# ==================================================

from sage.all import vector, ZZ, floor, matrix, ceil
from ...parameters.matrix_gestion import gen_overflow_matrix

def babai_nearest_plane_reduction(l_base, pol_p):
    """
    Reduce a polynomial using Babai's nearest plane algorithm.
    Code adapted from N. Meloni's Sage implementation.

    Args:
        pol_p (Polynomial): polynomial to reduce.
        l_base (matrix): lattice basis of null polynomials for the PMNS.

    Returns:
        Polynomial: reduced polynomial representing the same element.
    """

    n = l_base.nrows()

    # Gram-Schmidt orthogonalization of the lattice basis
    gram, _ = l_base.gram_schmidt()

    # normalized Gram vectors
    G = [gram[i] / (gram[i].norm()**2) for i in range(n)]

    s = vector(pol_p.list() + [0]*(n - pol_p.degree() - 1))

    # reduction
    for idx in range(n-1, -1, -1):
        coef = round(s.dot_product(G[idx]))
        s -= coef * l_base[idx]

    return ZZ["X"](list(s))


def babai_rounding_unlimited_reduction(l_base, pol_p):
    """
    Reduce a polynomial using Babai rounding with unlimited precision.
    Code adapted from N. Meloni's Sage implementation.

    Args:
        pol_p (Polynomial): polynomial to reduce.
        l_base (matrix): lattice basis of null polynomials for the PMNS.

    Returns:
        Polynomial: reduced polynomial representing the same element.
    """

    n = l_base.nrows()

    s = vector(pol_p.list() + [0]*(n - pol_p.degree() - 1))

    l_inv = l_base.inverse()

    coefs = vector(map(round, s * l_inv))

    return ZZ["X"](list(s - coefs * l_base))


def babai_rounding_limited_reduction(pol_p, h1, h2, l_base, l_inv_babai):
    """
    Reduce a polynomial using Babai rounding with limited precision.
    Code adapted from N. Meloni's Sage implementation.

    This version is designed for architectures with fixed precision
    (e.g., C implementations).

    Args:
        pol_p (Polynomial): polynomial to reduce.
        h1 (int): scaling parameter for Babai rounding.
        h2 (int): secondary scaling parameter for coefficient truncation.
        l_base (matrix): lattice basis of null polynomials for the PMNS.
        l_inv_babai (matrix): precomputed scaled inverse of the lattice basis.

    Returns:
        Polynomial: reduced polynomial representing the same element.
    """

    n = l_base.nrows()

    resize_vect = pol_p.list() + [0]*(n - len(pol_p.list()))

    v = vector([floor(x / 2**h2) for x in resize_vect])

    s = v * l_inv_babai
    s = vector([floor(x / 2**(h1 - h2)) for x in s])

    return ZZ["X"](list(vector(resize_vect) - s * l_base))


def gen_params_for_babai(l_base, phi:int, rho:int, pol_e):
    """
    Generate parameters for babi limited implementation.

    Args:
        l_base (matrix): lattice base of null polynomial in gamma
        phi (int): represent bit word size use by the architecture
        rho (int): limit for absolute coefficients values
        pol_e (Polynomial): polynomial use in external reduction

    Returns:
        h1 (int): scaling parameter for Babai rounding.
        h2 (int): secondary scaling parameter for coefficient truncation.
        l_inv_babai (matrix): inver in ratianal field of the lattice base l_base
    """
    l_inv = l_base.inverse()

    # h1 must verify that round(2**h1 * B^-1) <= phi/2
    # so 2**h1 * B^-1 <= phi/2 - 0.5
    # so h1 <= log2((phi/2 - 0.5) / B^-1)
    # as this condition is true for all element and knowing that the limit is at its lower when B^-1 = max(|B^-1|)

    maximum_value = max(abs(c) for c in l_inv)
    h1 = ceil((phi/2 - 0.5) / (maximum_value)).nbits()
    
    l_inv_babai = matrix([[round(2**h1 * x) for x in vect] for vect in l_inv])

    # h2 must verify that low(v / 2**h2) <= phi/2
    # so v / 2**h2 <= phi/2 <==> h2 > log2(2*v/phi)
    # this contion must be true for all element and thus for max(v) giving us our lower limit of h2
    # but element are under rho so the product (without reduction) have coefficients under rho**2
    # with reduction, overhead * rho**2 is added at maximum to element
    # so maximum element are less than rho**2(overhead + 1)
        
    epsilon = gen_overflow_matrix(pol_e)
    maximum_value = rho**2 * (1 + epsilon.norm(1))
    h2 = int(2*maximum_value / phi).bit_length()
    
    return h1, h2, l_inv_babai