# ==================================================
# params_generation.py
# Generic functions around parameters generation for
# PMNS construction
# ==================================================

from sage.all import vector, infinity, ZZ, matrix, GF, gcd, PolynomialRing, xgcd, ceil, exp
from pmns_core.parameters.matrix_gestion import gen_overflow_matrix, gen_reduce_null_base
from pmns_core.math_utils import square_and_multiply

PR = PolynomialRing(ZZ, "X")
PR2 = PolynomialRing(GF(2), "X")

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
        phi (int): word size bound
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



def search_m_with_even_degs(base, pol_e):
    """
    Optimized search for an invertible polynomial modulo pol_e
    by looking for elements whose degree-0 coefficient is odd.

    Args:
        base (matrix): matrix of null polynomials evaluated at gamma
        pol_e (Polynomial): polynomial used for external PMNS reduction

    Returns:
        Polynomial: polynomial invertible modulo pol_e and null over gamma
    """
    
    # loop over rows to find a polynomial with odd constant term
    for row in base:
        if row[0] & 1:  # check if constant term is odd
            return PR(list(row)) % pol_e
    return None


def search_m_with_odd_deg(k: int, p: int, gamma, pol_e):
    """
    General search for an invertible polynomial modulo pol_e.

    Args:
        k (int): extension degree
        p (int): prime used to construct the extension field
        gamma: root of the external reduction polynomial
        pol_e (Polynomial): external reduction polynomial

    Returns:
        Polynomial: polynomial invertible modulo pol_e
    """
    n = pol_e.degree()
    base = matrix(ZZ, n, n, 0)
    
    # precompute element to fill the matrix
    gamma_pk = -int((gamma**k))
    gamma_pk_mod2 = gamma_pk & 1
    
    # use the same process as for classical PLNS representation
    # by manipulating base coefficients and adding odd coefficient if gamma_pk is odd
    coef = gamma_pk + p * gamma_pk_mod2

    for i in range(k):
        base[i, i] = p

    for i in range(k, n):
        base[i, i] = 1
        base[i, i - k] = coef

    reduced_base = base.LLL()
    
    # Note : with this structure, there always exists an element invertible in this base if E is irrecdutible.
    # we evaluate polynomial in an extension field of characteristic p. Our construction therefore
    # doesn't change the polynomial value when evaluated in the extension field.
    
    # Note: This function is a general search for the inverse but can be improved if all coefficients 
    # are even, by applying the search of the function 'search_m_with_even_deg'
    
    reference_polynomial = PR2(pol_e)

    # search for a linear combination that is invertible mod pol_e
    for linear_combination in range(1, 2**n):
        # do polynomial combination using bit representation
        result = sum(reduced_base[i] for i in range(n) if (linear_combination >> i) & 1)
        
        # create the polynomial and cast it to GF(2)
        polynomial = PR(list(result)) % pol_e
        bin_polynomial = PR2(polynomial)

        # check that the polynomial is invertible and not a constant polynomial
        if bin_polynomial != 1 and gcd(bin_polynomial, reference_polynomial) == 1:
            return polynomial
    return None


def search_polynomial_m(base, k:int, p:int, gamma, pol_e):
    """
    Function which searches for an invertible element modulo pol_e with optimization
    if the degree-0 coefficient of pol_e is even.

    Args:
        base (matrix): matrix of null polynomials when evaluated over gamma
        k (int): extension degree
        p (int): prime used to construct the extension field
        gamma: root of the external reduction polynomial
        pol_e (Polynomial): external reduction polynomial
        
    Returns:
        Polynomial: polynomial invertible modulo pol_e
    """
    no_optimisation = any(coef&1 for coef in pol_e)
    if no_optimisation:
        return search_m_with_odd_deg(k, p, gamma, pol_e)
    return search_m_with_even_degs(base, pol_e)


def search_m_and_n(k: int, p: int, gamma, base, pol_e, phi: int=64):
    """
    Function that retrieves a polynomial M invertible modulo pol_e and N = -M^(-1) mod phi.

    Args:
        k (int): extension degree
        p (int): prime used to construct the extension field
        gamma (extension field element): root of E suitable for PMNS construction
        base (matrix): reduced base of null polynomial over gamma
        phi_pow (int, Optional): word size bound. Equal to 64 by default
        pol_e (Polynomial): polynomial used for external reduction in PMNS

    Returns:
        Polynomial: M, a polynomial null over gamma and invertible modulo pol_e
        Polynomial: N, a polynomial such that N = -M^(-1) mod phi
    """
    # retrieve an invertible polynomial M modulo pol_e
    M = search_polynomial_m(base, k, p, gamma, pol_e)
    
    assert M(gamma) == 0, "problem occuring with the base. Polynomial M isn't inverssible"

    # with xgcd, we get d, u, v such that M*u + pol_e*v = d
    # modulo pol_e, we have M^(-1) = u * d^(-1)
    d, u, _ = xgcd(M, pol_e)

    # ensure d is integer and invertible modulo phi
    d = int(d)
    d_inv = pow(d, -1, phi)

    # compute -M^(-1) modulo phi
    N = PR((-d_inv*u) % phi)

    return M, N



def search_minimal_degree(p: int, k: int, phi_pow: int, max_add_coef: int) -> int:
    """
    Function that compute a minimal value of n such that we can possibly construct a PMNS

    Args:
        p (Interger): prime used to construction extension field
        k (int): extension degree
        phi_pow (int): word size
        max_add_coef (int): minimal value add to coefficient after internal reduction with initial parameters
        
        exemple : E = X^n - alpha X^k - beta => max_add_coef can be approximate by  beta + alpha(alpha + beta)

    Returns:
        int: return a degree n minimal for wich we can possibly construct a PMNS
    """
    pbits = p.nbits()
    n = int(pbits * k / phi_pow) + 1
    phi = 2**phi_pow

    # compute minimal degree n wich can lead to a possible contruction of PMNS
    # here :
    #   > p**(k/n) is due to Minkowski theorem saying that element are upper bound by p^(k/n)
    #       to avoid problem with high bit prime, we use an approximation of p by a power of 2
    #   > ceil(n / exp(1)) is due to stirling solving approximation of the hyperspheric condition 
    #       over the volume of the lattice (V.Pascale, N.Méloni, F.Palma)
    #   > factor 2 is due to Babai rounding: even with error, factor 2 allows stability for Babai algorithm
    #   > (max_add_coef*(n - 1) + 1) represent the value of overleaping coefficient after an external reduction 
    #       of the product of two theoretical PMNS polynomial
    
    while round(2**(k * pbits / n) * ceil(n / exp(1)) * 2 * (max_add_coef * (n - 1) + 1)) >= phi:
        n += 1
        
    return n

def cast_polynomial_to_minimal_representation(pol, p):
    """
    Retrun polynomial pol with coefficients in signed samll representation with modulus p.

    Args:
        p (Interger): prime used to construction extension field
        pol(Polynomial): polynomial that need to have reduced coefficients

    Returns:
        Polynomial: return polynomial with coefficients in [-p//2, p//2]
        
    Exemple : (p-1)*X => -X[p]
    """
    coefficients = [int(c) if c <= p//2 else int(c)-p for c in pol]
    return PR(coefficients)