# ==================================================
# conversions_gestion.py
# function wich permit to convert element to polynomial
# or/and pmns representation using Montgomery or Babai
# reduction
# ==================================================

from sage.all import vector, matrix, PolynomialRing, ZZ, log, RR, ceil, Integer
from pmns_factory.core.operations.reductions.montgomery_reduction import fast_montgomery_reduction

PR = PolynomialRing(ZZ, "X")

def gen_transition_matrix(gamma, k: int) -> matrix:
    """
    Generate transition matrix to express elements in gamma basis

    Args:
        gamma (extension field element): root of the external reduction
        k (int): degree of the extension

    Returns:
        matrix: transition matrix from canonical basis to gamma basis
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




def montgomery_pseudo_fast_conversion(element, container):
    k = container.get('k')
    E = container.get('E')
    L, L_inv =container.get('L_origin'), container.get('L_inv_origin')
    theta_pow = container.get('theta_pow')
    phi = 2**container.get('phi_pow')
    mask = (1<<theta_pow) - 1 
    rho, gamma = container.get('rho'), container.get('gamma')

    n_pols = container.get('n_pol')
    ipols = container.get('int_pols')[0]
    zpols = container.get('z_pols')
    
    polynomial = PR(0)
    for deg in range(k):
        current_element = int(element._vector_()[deg])
        current_pol = PR(0)

        for pol_coeffs in ipols:
            part = current_element & mask
            current_pol += part * PR(pol_coeffs)

            current_element >>= theta_pow

        polynomial += current_pol * PR(zpols[deg]) % E

    for i in range(container.get('n_red_pseudo')):
        polynomial = fast_montgomery_reduction(polynomial, L, L_inv, phi)
    
    assert polynomial(gamma) == element, f"polynomial doesn't represent {element=}\n{polynomial(gamma)}"
    assert all(abs(c) < rho for c in polynomial), f"{rho=} too low for {element=}\n{polynomial =}"
    
    return polynomial



def montgomery_fast_conversion(element, container):
    k = container.get('k')
    L, L_inv =container.get('L_origin'), container.get('L_inv_origin')
    theta_pow = container.get('theta_pow')
    pols = container.get('fast_pols')
    phi = 2**container.get('phi_pow')
    mask = (1<<theta_pow) - 1 
    rho, gamma = container.get('rho'), container.get('gamma')
    
    polynomial = PR(0)
    for deg in range(k):
        current_element = int(element._vector_()[deg])
        
        for pol_coeffs in pols[deg]:
            part = current_element & mask
            polynomial += part * PR(pol_coeffs)
            
            current_element >>= theta_pow
    polynomial = fast_montgomery_reduction(polynomial, L, L_inv, phi)
    
    assert polynomial(gamma) == element, f"polynomial doesn't represent {element=}\n{polynomial(gamma)}"
    assert all(abs(c) < rho for c in polynomial), f"{rho=} too low for {element=}\n{polynomial =}"
    
    return polynomial
    

def montgomery_exact_conversion(element, container, add_red=0):
    """
    Convert an extension field element to PMNS using Montgomery reduction.

    Args:
        element (extension field element): element to convert
        container(PMNSContainer): container used to contain generated pmns

    Returns:
        Polynomial: PMNS representation
    """
    
    phi_pow = container.get('phi_pow')
    gamma = container.get('gamma')
    rho = container.get('rho')
    L = container.get('L_origin')
    L_inv = container.get('L_inv_origin')
    transition_matrix = container.get('T_mat_origin')
    
    phi = 2**phi_pow
    nb_iteration = container.get('n_red_exact') + add_red
    alpha = element * phi**nb_iteration

    V = convert_element_to_polynomial(alpha, gamma, transition_matrix)

    for _ in range(nb_iteration):
        V = fast_montgomery_reduction(V, L, L_inv, phi)

    assert V(gamma) == element, f"polynomial doesn't represent {element=}\n{V(gamma)}"
    assert all(abs(c) < rho for c in V), f"{rho=} too low for {element=}\n{V =}"

    return V

def compute_nb_internal_reductions(num, phi, rho, sublattice):
    augment = phi / (phi - 1)
    norm1 = sublattice.norm(1)
    denom = rho - (norm1 / 2) * augment
    
    rr_phi = RR(phi)
    rr_denom = RR(denom)
    rr_num = RR(num)

    log_phi = log(rr_phi)
    log_x = log(rr_num)
    log_denom = log(rr_denom)
    result = ceil((log_x - log_denom) / log_phi)
    return result
    
    
    
def compute_conversion_tables(container, psi, nb_red, npols, over_field=True):
    n, k = container.get('n'), container.get('k')
    phi_red = (2**container.get('phi_pow'))**nb_red
    mask = 1 << psi
    z = container.get('gamma').parent().gen()
    
    num_deg = k if over_field else 1

    return [[montgomery_exact_conversion((mask**i) * phi_red * (z**deg), container, add_red=1).list() for i in range(npols)] for deg in range(num_deg)]
