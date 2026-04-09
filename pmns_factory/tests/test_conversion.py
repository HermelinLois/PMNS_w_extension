from sage.all import *
from pathlib import Path
import sys
from time import time

PMNS_FACTORY_DIR = Path(__file__).resolve().parents[1]
if str(PMNS_FACTORY_DIR) not in sys.path:
    sys.path.append(str(PMNS_FACTORY_DIR))

from core.operations.convertions_gestion import gen_transition_matrix, convert_element_to_polynomial
from core.parameters.roots_gestion import *
from core.parameters.params_gestion import search_memory_overhead
from core.operations.reductions.montgomery_reduction import fast_montgomery_reduction
from core.operations.convertions_gestion import convert_element_to_pmns_montgomery
import pmns_E_type0_optimised as otype0
import pmns_E_type0 as type0

PR = PolynomialRing(ZZ,"X")
X = PR.gen()

m = None
p = Integer(101744864283287535450564907935948129100692627616929520537342055149603313146187) #random_prime(2**m, lbound=2**(m-1))
k = 13
pmns = otype0.gen_parameters(p, k)

E = pmns['E']
L = pmns['L']
gamma = pmns['gamma']
rho = pmns['rho']
n = E.degree()
phi = 2**pmns['phi_pow']

K = gamma.parent()
delta = gamma**k
pmns.update({'L_inv': (-(L.inverse()%phi))})

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


def int_conversion_to_pmns(element, params_2pow):
    decomp = decompose_bits(element, params_2pow)
    P = convert_to_pol(params_2pow['pow2_pmns'] , decomp)
    return P


def decompose_bits(element, params_2pow):
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
    zeta = ceil((m+1)/((t-1) * coef))
    nb_elements = zeta * coef
    near_2pow = ceil((m + 1)/nb_elements)

    pow2_pmns = []
    for i in range(nb_elements):
        pow2_i = convert_element_to_pmns_montgomery(K([2**(near_2pow * i)]), transition_matrix, pmns)
        pow2_pmns.append(pow2_i)

    """    print("m = ", m)
    print("t = ", t)
    print("zeta = ", zeta)
    print("no added reduction ? ", zeta* n//k *((k-1)*w+1), rho, zeta* n//k *((k-1)*w+1) < rho)
    print("~n = ", n//k)"""

    return {'near_2pow': near_2pow, 'pow2_pmns': pow2_pmns, 'nb_elements': nb_elements}


#############################################################################################
#############################################################################################
#############################################################################################

head = "<====> CLASSICAL CONVERSION <====>"
print("\n"+head)

# this function is underperfoming due to equality check over gamma and check of
# coefficients under rhorandom_prime(2**m, lbound=2**(m-1)) # 

start = time()
P = convert_element_to_pmns_montgomery(element, transition_matrix, pmns)
stop = time()

print(f"conversion of element is : {P = } ({stop - start})")
print("=" * len(head))




#############################################################################################
#############################################################################################
#############################################################################################

head = "<====> PSEUDO-FAST CONVERSION (no transp. matrix) <====>"
print("\n" + head)

def gen_convertion_to_gamma_pols(pmns):
    E = pmns['E']
    L = pmns['L']
    L_inv = pmns['L_inv']
    gamma = pmns['gamma']
    p = pmns['p']
    k = pmns['k'] 
    phi = 2**pmns['phi_pow']
    
    conversion_pols = []
    transition_matrix = gen_transition_matrix(gamma, k)
    
    for i in range(k):
        z_pow_i = K([0]*i + [1])
        
        C_i = convert_element_to_pmns_montgomery(z_pow_i, transition_matrix, pmns)
        conversion_pols.append(C_i)
    
    for i, C_i in enumerate(conversion_pols):
        assert C_i(gamma) == pow(K([0,1]), i, p), f"C_{i}(gamma) = {C_i(gamma)} != z^{i}"
    
    return conversion_pols

def field_conversion_to_pmns_v1(element, params_2pow, conversion_pols, pmns):
    E = pmns['E']
    L = pmns['L']
    L_inv = pmns['L_inv']
    rho = pmns['rho']
    phi = 2**pmns['phi_pow']
    n = E.degree()
    k = pmns['k']
    w = search_memory_overhead(E)
    
    
    m_plus_1 = pmns['p'].nbits()  
    t = rho.nbits()
    zeta = ceil((m_plus_1) / ((t-1) * (n//k)))
    bound_P = zeta * n * w * rho**3
    nb_internal_reduction = ceil(log(zeta * n * w * rho**3/(rho - L.norm(1)/2 * phi/(phi-1)), phi))

    
    element = element * phi**nb_internal_reduction
    P = PR(0)
    for i, x_i in enumerate(element._vector_()):
        Qi = int_conversion_to_pmns(x_i, params_2pow)        
        P += Qi * conversion_pols[i] % E
    
    for _ in range(nb_internal_reduction):
        P = fast_montgomery_reduction(P, E, L, L_inv, phi)
    
    return P

params_2pow = gen_2pow_parameters(transition_matrix, pmns, p)
conversion_pols = gen_convertion_to_gamma_pols(pmns)

start = time()
P = field_conversion_to_pmns_v1(element, params_2pow, conversion_pols, pmns)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element

print(f"\nconversion of element is : {P = } ({stop - start})")
print("=" * len(head))





#############################################################################################
#############################################################################################
#############################################################################################

head = "<====> PSEUDO-FAST CONVERSION (transp. matrix) <====>"
print("\n"+head)

def field_conversion_to_pmns_v2(element, transition_matrix, pmns, params_2pow):
    E = pmns['E']
    L = pmns['L']
    L_inv = pmns['L_inv']
    rho =pmns['rho']
    gamma = pmns['gamma']
    phi = 2**pmns['phi_pow']

    m_plus_1 = pmns['p'].nbits()  
    t = rho.nbits()
    zeta = ceil((m_plus_1) / ((t-1) * (n//k)))
    w = search_memory_overhead(E)

    nb_internal_reduction = ceil(log(zeta* n//k *((k-1)*w + 1) * rho**2/(rho - L.norm(1)/2 * (phi/(phi-1))), phi))
    element_in_gamma = convert_element_to_polynomial(element * phi**nb_internal_reduction, gamma, transition_matrix)
    
    P = 0
    for idx, e in enumerate(element_in_gamma):
        P += int_conversion_to_pmns(e, params_2pow) * X**idx % E
    
    for i in range(nb_internal_reduction):
        P = fast_montgomery_reduction(P, E, L, L_inv, phi)

    return P

params_2pow = gen_2pow_parameters(transition_matrix, pmns, p)

start = time()
P = field_conversion_to_pmns_v2(element, transition_matrix, pmns, params_2pow)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element

print(f"\nconversion of element is : {P = } ({stop - start})")
print("=" * len(head))

#############################################################################################
#############################################################################################
#############################################################################################
head = "<====> FAST CONVERSION (transp. matrix)<====>"
print("\n"+head)

def decompose_bits(element, nb_elements, near_2pow):
    current_element = Integer(element)
    result = []
    mask = (1 << near_2pow) - 1
    for _ in range(nb_elements):
        result.append(current_element & mask)
        current_element >>= near_2pow
    return result

def gen_2pow_gamma_i_parameters(transition_matrix, pmns, base_int):
    m = base_int.nbits()
    n = pmns['E'].degree()
    t = pmns['rho'].nbits()
    k = pmns['k']
    K = pmns['gamma'].parent()
    gamma = pmns['gamma']

    coef_div = n // k
    zeta = ceil((m + 1) / ((t - 1) * coef_div))
    nb_elements = zeta * coef_div 
    near_2pow = ceil((m + 1) / nb_elements)

    pow2_pmns_gammas = [[] for _ in range(k)]
    
    for j in range(k):
        gamma_j = gamma**j
        for i in range(nb_elements):
            val_to_convert = K(2**(near_2pow * i) * gamma_j)
            pow2_val = convert_element_to_pmns_montgomery(val_to_convert, transition_matrix, pmns)
            pow2_pmns_gammas[j].append(pow2_val)

    return {
        'pow2_pmns': pow2_pmns_gammas, 
        'nb_elements': nb_elements,
        'near_2pow': near_2pow
    }


def field_conversion_to_pmns_v3(element, transition_matrix, pmns, params_2pow_gamma_i):
    E = pmns['E']
    L = pmns['L']
    L_inv = pmns['L_inv']
    gamma = pmns['gamma']
    phi = 2**pmns['phi_pow']
    rho = pmns['rho']
    
    nb_red = 1 
    element_in_gamma = convert_element_to_polynomial(element * (phi**nb_red), gamma, transition_matrix)
    
    P = PR(0)
    for idx, e in enumerate(element_in_gamma):
        decomp = decompose_bits(e, params_2pow_gamma_i['nb_elements'], params_2pow_gamma_i['near_2pow'])
        table_gamma_idx = params_2pow_gamma_i['pow2_pmns'][idx]

        P += sum(c * table_gamma_idx[i] for i, c in enumerate(decomp))

    return fast_montgomery_reduction(P, E, L, L_inv, phi)



pow2_pmns_gammas = gen_2pow_gamma_i_parameters(transition_matrix, pmns, p)
start = time()
P = field_conversion_to_pmns_v3(element, transition_matrix, pmns, pow2_pmns_gammas)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element, f'{P(gamma)} \n{element}'

print(f"\nconversion of element is : {P = } ({stop - start})")
print("=" * len(head))

head = "<====> FAST CONVERSION (no transp. matrix)<====>"
print("\n"+head)

def gen_2pow_gamma_i_parameters(transition_matrix, pmns, base_int):
    m = base_int.nbits()
    n = pmns['E'].degree()
    t = pmns['rho'].nbits()
    k = pmns['k']
    K = pmns['gamma'].parent()
    gamma = pmns['gamma']

    coef_div = n // k
    zeta = ceil((m + 1) / ((t - 1) * coef_div))
    nb_elements = zeta * coef_div 
    near_2pow = ceil((m + 1) / nb_elements)

    pow2_pmns_gammas = [[] for _ in range(k)]
    
    for j in range(k):
        z_j = K([0]*j + [1])
        for i in range(nb_elements):
            val_to_convert = K(2**(near_2pow * i) * z_j)
            pow2_val = convert_element_to_pmns_montgomery(val_to_convert, transition_matrix, pmns)
            pow2_pmns_gammas[j].append(pow2_val)

    return {
        'pow2_pmns': pow2_pmns_gammas, 
        'nb_elements': nb_elements,
        'near_2pow': near_2pow
    }


def field_conversion_to_pmns_v3(element, transition_matrix, pmns, params_2pow_gamma_i):
    E = pmns['E']
    L = pmns['L']
    L_inv = pmns['L_inv']
    gamma = pmns['gamma']
    phi = 2**pmns['phi_pow']
    element *= phi
    
    P = PR(0)
    for idx, e in enumerate(element._vector_()):
        decomp = decompose_bits(e, params_2pow_gamma_i['nb_elements'], params_2pow_gamma_i['near_2pow'])
        table_gamma_idx = params_2pow_gamma_i['pow2_pmns'][idx]

        P += sum(c * table_gamma_idx[i] for i, c in enumerate(decomp))

    return fast_montgomery_reduction(P, E, L, L_inv, phi)



pow2_pmns_gammas = gen_2pow_gamma_i_parameters(transition_matrix, pmns, p)
start = time()
P = field_conversion_to_pmns_v3(element, transition_matrix, pmns, pow2_pmns_gammas)
stop = time()

assert all(abs(c) < rho for c in P)
assert P(gamma) == element, f'{P(gamma)} \n{element}'

print(f"\nconversion of element is : {P = } ({stop - start})")
print("=" * len(head))

