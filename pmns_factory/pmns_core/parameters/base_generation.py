# ==================================================
# base_generation.py
# Generic function for PMNS base generation
# ==================================================

from sage.all import matrix, ZZ

def gen_reduce_null_base(k:int, p:int, n:int, gamma) -> matrix:
    """
    Generate a reduced base of null polynomial for gamma.
    Approach based on: An Alternative Approach for SIDH Arithmetic (C.Bouvier & L.Imbert)
    
    Args:
        p (int): prime used to construct the extension field
        k (int): extension degree of the field
        n (int): degree of E used for polynomial reduction in PMNS
            gamma: element used to construct PMNS (gamma^k is integer)
    
    Returns:
        matrix (LLL) : reduced base where each row evaluates to zero at gamma
    """

    gamma_pow_k = -gamma**k
    base = matrix(ZZ, n, n, 0)

    # fill the diagonal
    for i in range(k):
        base[i, i] = p
    for i in range(k, n):
        base[i, i] = 1
        base[i, i-k] = gamma_pow_k

    # LLL reduction
    return base.LLL()

