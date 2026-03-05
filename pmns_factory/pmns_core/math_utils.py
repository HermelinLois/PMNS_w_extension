# ==================================================
# math_utils.py
# Generic math functions for PMNS parameter generation
# ==================================================

from sage.all import Integer, gcd

def square_and_multiply(base, exponent, mod=None):
    """
    Compute base raised to the power of exponent efficiently with modulus.

    Args:
        base : The base to be exponentiated
        exponent : the exponent
        mad : the modulus

    Returns:
        result (same type as base) : value of base ** exponent
    """
    result = 1
    while exponent:
        if exponent & 1:
            result *= base
        base = base ** 2
        exponent >>= 1
    return result if mod is None else result % mod
