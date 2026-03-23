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


m = 2048
p = random_prime(2**m, lbound=2**(m-1))
k = 3
pmns = otype0.gen_parameters(p, k)

E = pmns['E']

gamma = pmns['gamma']
#print("gamma : ", gamma)
K = gamma.parent()
n = E.degree()

PR = PolynomialRing(K,"X")
X = PR.gen()
Ek = PR(E)
nE = X**9 - X**3 - 8
transition_matrix = gen_transition_matrix(gamma, k)

def decompose_centered(target, base, p):
    coeffs = []
    K = GF(p)
    current = target
    for v in base:
        c, current = divmod(current, v)
        coeffs.append(K(c).lift_centered())

    assert current == 0

    return coeffs


def gen_lattice(gamma, k, p, n):
    nb_coefs = n / k
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

def gen_convertion_pol_basis(gamma, p, k, E, L):
    PR = PolynomialRing(ZZ, "X")
    X = PR.gen()
    n = E.degree()

    mu = max(vector(ZZ, gamma._vector_()), key=abs)
    mu_inv = Integer(pow(mu,-1, p))

    base1 = gen_lattice(gamma, k, p, n)
    
    decomposition1 = decompose_babai(mu_inv, base1)
    P1 = convert_list_to_pol(decomposition1, X, k)

    print(P1)

def convert_to_polynomial():
    pass
L = pmns['L']
print('E = ', E)
print('rho=', pmns['rho'])

r = gen_convertion_pol_basis(gamma, p, k, E, L)