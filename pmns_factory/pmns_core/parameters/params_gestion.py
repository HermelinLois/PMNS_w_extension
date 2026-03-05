# ==================================================
# params_generation.py
# Generic functions around parameters geenration for
# PMNS construction
# ==================================================

from sage.all import vector, infinity, ceil, exp
from pmns_factory.pmns_core.parameters.matrix_gestion import gen_overflow_matrix, gen_reduce_null_base

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
        phi (int): word size bound (typically 2^word_size)
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
        base_norm = int(base.norm(1))

        # Coefficient bound estimation
        rho = 2 * w * base_norm

        # Ensure coefficients fit machine word
        if rho < phi:
            return base, rho, gamma

    return None