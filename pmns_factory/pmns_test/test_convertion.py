from sage.all import *
from pmns_core.operations.convertions_gestion import gen_transition_matrix, convert_element_to_polynomial
from pmns_core.parameters.roots_gestion import *
from pmns_core.parameters.params_gestion import search_m_and_n
import pmns_E_type0_optimised as otype0
import pmns_E_type0 as type0


m = 512
p = random_prime(2**m, lbound=2**(m-1))
k = 5
pmns = otype0.gen_parameters(p, k)

E = pmns['E']
gamma = pmns['gamma']
#print("gamma : ", gamma)
K = gamma.parent()
n = E.degree()

PR = PolynomialRing(K,"X")
Ek = PR(E)

#print(E)

"""rts = select_roots(Ek.roots(), p, k)

for root in rts :
    if is_root_pow_in_base_field(root, k, p) and is_root_free(root, k, p):
        rk = root**k
        #print([Integer(rk**i).nbits() for i in range((n-1)//k + 1)])
        

M,N = search_m_and_n(k, p, pmns['gamma'], pmns['L'], pmns['E'])
pmns.update({"M": M, "N": N})"""



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


def gen_lattice(gamma, k, p, n):
    nb_coefs = int((n-1) / k) + 1
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
    s = vector(target + [0]*(n - len(target)))
    l_inv = L.inverse()
    coefs = vector(map(round, s * l_inv))
    
    return [(c,i) for i,c in enumerate(list(s - coefs * L))]



def convert_list_to_pol(coefs, X, k):
    return sum(c*X**(k*deg) for c, deg in coefs)

def gen_convertion_pol_basis(gamma, rho, p, k, E):
    PR = PolynomialRing(ZZ, "X")
    X = PR.gen()
    n = E.degree()

    mu = max(vector(ZZ, gamma._vector_()), key=abs)
    mu_inv = Integer(pow(mu,-1, p))

    base = gen_lattice(gamma, k, p, n)
    decomposition = decompose_babai([mu_inv], base)
    print(decomposition)
    for i in range(k):
        decomposition = decompose_babai([c for c,_ in decomposition], base)
        print(decomposition, convert_list_to_pol(decomposition, X, k)(gamma)==mu_inv)

    """pol_mu_inv = convert_list_to_pol(decomposition, X, k)

    assert pol_mu_inv(gamma) == mu_inv
    print(rho)

    pol_base = [0] * k
    for deg in range(k):
        pol_base[deg] = pol_mu_inv**deg % E

    return pol_base
"""
def convert_to_polynomial():
    pass


r = gen_convertion_pol_basis(gamma, pmns['rho'], p, k, E)
print(p)
print(r)