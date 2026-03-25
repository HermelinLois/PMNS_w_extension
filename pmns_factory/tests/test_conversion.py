from sage.all import *
from pathlib import Path
import sys
from time import time

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

PR = PolynomialRing(ZZ,"X")
X = PR.gen()

m = 1024
p = random_prime(2**m, lbound=2**(m-1))
k = 2
pmns = otype0.gen_parameters(p, k)


E = pmns['E']
L = pmns['L']
gamma = pmns['gamma']
rho = pmns['rho']
n = E.degree()

K = gamma.parent()
delta = gamma**k
mu = max(vector(ZZ, gamma._vector_()), key=abs)

M, N = search_m_and_n(k, p, gamma, L, E)
pmns.update({'M': M, 'N': N})

element = K([randint(0,p-1) for _ in range(k)])

print("p = ",p)
print("E = ",E)
print("gamma = ", gamma)
print("rho = ", rho)
print("element = ", element)

transition_matrix = gen_transition_matrix(pmns['gamma'], k)

# common functions 
def convert_to_pol(base:list, nb_coefs:list):
    return sum(c*base[i] for i, c in enumerate(nb_coefs))


def fast_int_conversion_to_pmns(element, params_2pow):
    decomp = fast_decompose_bits(element, params_2pow)
    P = convert_to_pol(params_2pow['pow2_pmns'] , decomp)
    return P


def fast_decompose_bits(element, params_2pow):
    pow2_bits = params_2pow['near_2pow']
    nb_elements = params_2pow['nb_elements']

    current_element = Integer(element)
    result = []
    mask = (1<<pow2_bits) - 1

    for i in range(nb_elements):
        digit = mask & current_element
        result.append(digit)
        current_element >>= pow2_bits

    return result


def gen_2pow_parameters(transition_matrix, pmns, base_int):
    m = base_int.nbits()
    n = pmns['E'].degree()
    t = pmns['rho'].nbits()
    k = pmns['k']

    coef = n//k
    augment = ceil((m+1)/((t-1) * coef))
    nb_elements = augment * coef
    near_2pow = ceil((m + 1)/nb_elements)

    pow2_pmns = []
    for i in range(nb_elements):
        pow2_i = convert_element_to_pmns_montgomery(K([2**(near_2pow * i)]), transition_matrix, **pmns)
        pow2_pmns.append(pow2_i)

    print("m = ", m)
    print("t = ", t)
    print("M = ", augment)
    print("under rho = ", floor((m+1)/((t-1) * coef)) + 1)
    print("~n = ", n//k)

    return {'near_2pow': near_2pow, 'pow2_pmns': pow2_pmns, 'nb_elements': nb_elements}

#############################################################################################
#############################################################################################
#############################################################################################

head = "<====> CLASSICAL CONVERSION <====>"
print("\n"+head)

# this function as a slowdown of his capacity due to equality check over gamma and check of
# coefficients under rho
start = time()
P = convert_element_to_pmns_montgomery(element, transition_matrix, **pmns)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element

print(f"conversion of element is : {P = } ({stop - start})")
print("=" * len(head))


#############################################################################################
#############################################################################################
#############################################################################################

head = "<====> FIRST CONVERSION (no transp. matrix) <====>"
print("\n" + head)

def gen_convertion_to_gamma_pols(transition_matrix, pmns, mu):
    p = pmns['p']
    k = pmns['k']
    E = pmns['E']
    M = pmns['M']
    N = pmns['N']
    gamma = pmns['gamma']
    phi = 2**pmns['phi_pow']

    mu_inv = K(mu)**(-1)
    pmns_mu_inv = convert_element_to_pmns_montgomery(mu_inv, transition_matrix, **pmns)

    conversion_pols = []
    for i in range(k):
        current_pol = pmns_mu_inv**i * X**i * phi % E
        current_pol = montgomery_reduction(current_pol, M, N, E, gamma)
        conversion_pols.append(current_pol)
    return conversion_pols


def fast_field_conversion_to_pmns(element, params_2pow, conversion_pols, pmns):
    E = pmns['E']
    M = pmns['M']
    N = pmns['N']
    gamma = pmns['gamma']
    rho = pmns['rho']
    phi = 2**pmns['phi_pow']

    nb_internal_reduction = 0 if element.polynomial().degree() <= 0 else 3
    elements = element._vector_() * phi**nb_internal_reduction

    P = 0
    for idx, elmt in enumerate(elements):
        P += fast_int_conversion_to_pmns(elmt, params_2pow) * conversion_pols[idx] % E

    for _ in range(nb_internal_reduction):
        P = montgomery_reduction(P, M, N, E, gamma)
    return P


params_2pow = gen_2pow_parameters(transition_matrix, pmns, p)
conversion_pols = gen_convertion_to_gamma_pols(transition_matrix, pmns, mu)

start = time()
P = fast_field_conversion_to_pmns(element, params_2pow, conversion_pols, pmns)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element

print(f"\nconversion of element is : {P = } ({stop - start})")
print("=" * len(head))

#############################################################################################
#############################################################################################
#############################################################################################

head = "<====> SECOND CONVERSION (transp. matrix) <====>"
print("\n"+head)

def fast_field_conversion_to_pmns(element, transition_matrix, pmns, params_2pow):
    E = pmns['E']
    M = pmns['M']
    N = pmns['N']
    rho =pmns['rho']
    gamma = pmns['gamma']
    phi = 2**pmns['phi_pow']

    nb_internal_reduction = 2
    element_in_gamma = convert_element_to_polynomial(element * phi**nb_internal_reduction, gamma, transition_matrix)
    
    P = 0
    for idx, e in enumerate(element_in_gamma):
        P += fast_int_conversion_to_pmns(e, params_2pow) * X**idx % E
    
    for i in range(nb_internal_reduction):
        P = montgomery_reduction(P, M, N, E, gamma)

    return P

params_2pow = gen_2pow_parameters(transition_matrix, pmns, p)
conversion_pols = gen_convertion_to_gamma_pols(transition_matrix, pmns, mu)

start = time()
P = fast_field_conversion_to_pmns(element, transition_matrix, pmns, params_2pow)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element

print(f"\nconversion of element is : {P = } ({stop - start})")
print("=" * len(head))
