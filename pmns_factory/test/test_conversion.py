from sage.all import *
from pathlib import Path
import sys

PMNS_FACTORY_DIR = Path(__file__).resolve().parents[1]
if str(PMNS_FACTORY_DIR) not in sys.path:
    sys.path.append(str(PMNS_FACTORY_DIR))

from core.operations.convertions_gestion import gen_transition_matrix, convert_element_to_polynomial
from core.parameters.roots_gestion import *
from core.parameters.params_gestion import search_m_and_n
from core.operations.reductions.montgomery_reduction import montgomery_reduction
from core.operations.convertions_gestion import convert_element_to_pmns_montgomery
import pmns_E_type0_optimised as otype0
import pmns_E_type0 as type0


m = 128
p = random_prime(2**m, lbound=2**(m-1))
k = 3
pmns = otype0.gen_parameters(p, k)

E = pmns['E']
L = pmns['L']
gamma = pmns['gamma']
n = E.degree()
delta = gamma**k
mu = max(vector(ZZ, gamma._vector_()), key=abs)
K = gamma.parent()
M,N = search_m_and_n(k, p, gamma, pmns['L'], E)
pmns.update({'M': M, 'N': N})

PR = PolynomialRing(ZZ,"X")
X = PR.gen()

print("p = ",p)
print("E = ",E)
print("rho = ", pmns['rho'])



# fast convertion implementation 
def gen_parameters(transition_matrix, pmns):
    m = pmns['p'].nbits()
    n = pmns['E'].degree()
    w = pmns['rho'].nbits()
    k = pmns['k']

    init_coef = n//k
    coefs_augment = floor((m+1)/(2 * w * init_coef)) + 1
    nb_coefs = coefs_augment * init_coef
    near_2pow = ceil((m + 1)/(init_coef * coefs_augment))

    pow2_pmns = []
    for i in range(nb_coefs):
        pow2_i = convert_element_to_pmns_montgomery(K([2**(near_2pow * i)]), transition_matrix, **pmns)
        pow2_pmns.append(pow2_i)

    print("m = ", m)
    print("M = ", coefs_augment)
    print("n = ", n//k)

    return {'base': pow2_pmns, 'near_2pow': near_2pow, 'augment': coefs_augment, 'nb_coefs': nb_coefs, 'init': init_coef}


def gen_convertion_to_gamma_pols(transition_matrix, pmns, mu):
    p = pmns['p']
    k = pmns['k']
    E = pmns['E']
    phi = 2**pmns['phi_pow']

    mu_inv = K(mu)**(-1)
    pmns_mu_inv = convert_element_to_pmns_montgomery(mu_inv, transition_matrix, **pmns)

    conversion_pols = []
    for i in range(k):
        current_pol = (pmns_mu_inv**i * X**i * phi) % p
        current_pol = montgomery_reduction(current_pol % E, M, N, E, gamma)
        conversion_pols.append(current_pol)
    return conversion_pols


def convert_to_pol(base:list, nb_coefs:list):
    return sum(c*base[i] for i, c in enumerate(nb_coefs))


def fast_decompose_bits(element, n, base_pow):
    current_element = Integer(element)
    result = [0]*n
    mask = (1<<base_pow) - 1
    for i in range(n):
        digit = mask & current_element
        result[i] = digit
        current_element >>= base_pow
    return result


def fast_int_conversion_to_pmns(element, pow2_base, phi_deg, pmns):
    phi = 2**pmns['phi_pow']
    p = pmns['p']

    decomp = fast_decompose_bits((element * phi**phi_deg) % p, pow2_base['nb_coefs'], pow2_base['near_2pow'])
    P = convert_to_pol(pow2_base['base'] , decomp)

    return P


def fast_field_conversion_to_pmns(element, pow2_base, conversion_pols, pmns):
    E = pmns['E']
    M = pmns['M']
    N = pmns['N']
    gamma = pmns['gamma']
    rho = pmns['rho']
    elements = element._vector_()
    
    phi_deg = 1 if element.polynomial().degree() <= 0 else 4

    P = 0
    for i, e in enumerate(elements):
        P += (fast_int_conversion_to_pmns(e, pow2_base, phi_deg, pmns) * conversion_pols[i]) % E %p

    for _ in range(phi_deg):
        P = montgomery_reduction(P, M, N, E, gamma)

    assert all(abs(c) < rho for c in P)
    assert P(gamma) == element

    return P



transition_matrix = gen_transition_matrix(pmns['gamma'], k)
pow2_base = gen_parameters(transition_matrix, pmns)
conversion_pols = gen_convertion_to_gamma_pols(transition_matrix, pmns, mu)

P = fast_field_conversion_to_pmns(K([mu, 1]), pow2_base, conversion_pols, pmns)
print(P)