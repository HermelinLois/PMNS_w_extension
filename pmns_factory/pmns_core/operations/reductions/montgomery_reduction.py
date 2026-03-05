# ==================================================
# montgomery_reduction.py
# function wich permit to apply Montgomery reduction 
# to a polynomial
# ==================================================

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