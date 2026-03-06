# ==================================================
# test_generation.py
# File to test a PMNS construction with specific E type
# ==================================================

from sage.all import random_prime, randint, GF
from pmns_E_type0 import gen_parameters
from pmns_core.parameters.params_gestion import search_m_and_n
from pmns_core.operations.convertions_gestion import gen_gamma_base, convert_element_to_pmns_montgomery
from pmns_core.operations.reductions.montgomery_reduction import montgomery_reduction

if __name__ == "__main__":
    m = 64
    p = 11038636823290202869#random_prime(2**m - 1, lbound=2**(m-1))
    k = 2


    print("=> start PMNS parameters generation <=")
    parameters = gen_parameters(p, k)
    print("=> generation achieved <=")

    rho, gamma, phi_pow, L, E = parameters['rho'], parameters['gamma'], parameters['phi_pow'], parameters['L'], parameters['E']
    phi = 2**phi_pow
    M, N = search_m_and_n(k, p, gamma, L, E, phi)
    parameters.update({'M': M, 'N': N})
    
    print(f"\n{p = }\n{k = }\n{E = }\n{M = }\n{N = }\n{gamma = }\n{phi_pow = }\n")

    K = gamma.parent()
    a = K([randint(0, p-1) for _ in range(k)])
    b = K([randint(0, p-1) for _ in range(k)])

    gamma_base = gen_gamma_base(gamma, k) 

    A = convert_element_to_pmns_montgomery(a, gamma_base, **parameters)
    B = convert_element_to_pmns_montgomery(b, gamma_base, **parameters)
    print(f"{a = } is represented by:\n{A}\n({A(gamma) == a})\n")
    print(f"{b = } is represented by:\n{B}\n({B(gamma) == b})\n")
    
    C = A * B % E
    rC = montgomery_reduction(C, M, N, E, phi)
    print(f"product before reduction is:\n{C}\n")
    print(f"product after reduction is:\n{rC}\ncoefficients are under {rho = }? ({all(abs(c)<rho for c in rC)})")
