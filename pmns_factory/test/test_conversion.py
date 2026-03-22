from sage.all import *
from pathlib import Path
import sys

# Allow running this file directly with: sage test_conversion.py
PMNS_FACTORY_DIR = Path(__file__).resolve().parents[1]
if str(PMNS_FACTORY_DIR) not in sys.path:
    sys.path.append(str(PMNS_FACTORY_DIR))

from core.operations.convertions_gestion import gen_transition_matrix, convert_element_to_polynomial
from core.parameters.roots_gestion import *
from core.parameters.params_gestion import search_m_and_n
import pmns_E_type0_optimised as otype0
import pmns_E_type0 as type0


m = 128
p = random_prime(2**m, lbound=2**(m-1))
k = 3
pmns = otype0.gen_parameters(p, k)

E = pmns['E']
gamma = pmns['gamma']
#print("gamma : ", gamma)
K = gamma.parent()
n = E.degree()

PR = PolynomialRing(K,"X")
Ek = PR(E)

transition_matrix = gen_transition_matrix(gamma, k)
params_head = "<=== PARAMETERS ===>"
#print(params_head)
#print("p = ",p)
#print("k = ", k)
#print("gamma = ",gamma)
#print("E = ", E)
#print("<" + "="*(len(params_head) - 2) + ">\n")

method1_head = "<=== METHOD MATRIX PARAMETERS ===>"
#print(method1_head)
#print("transition matrix = ")
#print(transition_matrix)
#print("<" + "="*(len(method1_head) - 2) + ">\n")

method2_head = "<=== METHOD INT PMNS PARAMETERS ===>"
#print(method2_head)

def decompose_centered(target, base, p):
    coeffs = []
    K = GF(p)
    current = target
    for v in base:
        c, current = divmod(current, v)
        coeffs.append(K(c).lift_centered())

    assert current == 0

    return coeffs


def gen_lattice(gamma, k, p, n,coef=1):
    nb_coefs = coef * (int(n / k) + 1)
    M = matrix(ZZ, nb_coefs, nb_coefs, 0)
    M[0,0] = p
    
    delta = pow(gamma,k)
    current = 1
    for i in range(1, nb_coefs):
        current *= delta
        M[i, i] = 1
        M[i, 0] = -current
    return M.LLL()


def decompose_babai(target, L):
    n = L.nrows()
    s = vector([target] + [0]*(n - 1))
    l_inv = L.inverse()
    coefs = vector(map(round, s * l_inv))

    return list(s - coefs * L)


def convert_list_to_pol(coefs, X, k):
    return sum(c * X**(k * deg) for deg, c in enumerate(coefs))

def gen_convertion_pol_basis(gamma, p, k, E):
    PR = PolynomialRing(ZZ, "X")
    X = PR.gen()
    n = E.degree()

    mu = max(vector(ZZ, gamma._vector_()), key=abs)
    mu_inv = Integer(pow(mu,-1, p))

    base1 = gen_lattice(gamma, k, p, n, coef=1)
    base2 = gen_lattice(gamma, k, p, n, coef=4)
    
    c1 = PR(list(base1[1]))
    c2 = PR(E*X**(1))
    print(list(c2))
    
    print(c1)
    print(c2)
    print()
    
    
    decomposition1 = decompose_babai(mu_inv, base1)
    decomposition2 = decompose_babai(mu_inv, base2)
    
    P1 = convert_list_to_pol(decomposition1, X, k)
    P2 = convert_list_to_pol(decomposition2, X, k) 
    assert P1(gamma) == mu_inv
    assert P2(gamma) == mu_inv
    assert (P2%c2)(gamma == mu_inv)

    nP2 = P2%c2
    print(P1)
    print(P2)
    print(nP2)
    print(ZZ(nP2[n]*gamma**k+ nP2[0]) -p)
    print(mu_inv)
    return P1

""" assert pol_mu_inv(gamma) == mu_inv

    pol_base = [0] * k
    for deg in range(k):
        pol_base[deg] = (pol_mu_inv**deg % E) % p

    return pol_base
"""
def convert_to_polynomial():
    pass

print('E = ', E)
print('rho=', pmns['rho'])

r = gen_convertion_pol_basis(gamma, p, k, E)