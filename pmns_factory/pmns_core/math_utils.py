# ==================================================
# math_utils.py
# Generic math functions for PMNS parameter generation
# ==================================================

from sage.all import Integer, gcd



def square_and_multiply(base, exponent):
    """
    Compute base raised to the power of exponent efficiently.

    Args:
        base (int, Field element or extension field element) : The base to be exponentiated
        exponent (int): the exponent

    Returns:
        result (same type as base) : value of base ** exponent
    """
    
    result = 1
    while exponent:
        if exponent & 1:
            result *= base
        base = base ** 2
        exponent >>= 1
    return result


def is_gamma_feasible(p, k) -> bool:
    """
    Check if a suitable gamma possibly exist for PMNS construction
    
    Args:
        p (int, Integer): prime use to construct field
        k (int, Integer): extension degree use to field

    Returns:
        bool: return True if there possibly exist gamma such that
            gamma isn't an interger and gamma^k is an interger
    """
    
    # Ensure p and k are Sage Integer to handle large primes
    p, k = Integer(p), Integer(k)
    value = (square_and_multiply(p, k) - 1) // (p-1)
    
    return gcd(k, value) > 1
    