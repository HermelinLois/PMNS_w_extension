from sage.all import randint, ceil, log
from jinja2 import Environment, FileSystemLoader
from pmns_factory.core.operations.reductions.montgomery_reduction import search_m_and_n, fast_montgomery_reduction
from pmns_factory.core.operations.convertions_gestion import convert_element_to_pmns_montgomery, gen_transition_matrix, convert_element_to_polynomial
import pmns_generator.writers.format_to_c.int_to_c as fint
from pathlib import Path
import inspect

CURRENT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = CURRENT_DIR / "templates" / "test_values_templates"


def write_reduction_test(output_dir:str , n_test:int, reduction_method: callable,  pmns_params:dict) -> None:
    # Create test_values subdirectory
    test_dir = Path(output_dir) / "test_values"
    test_dir.mkdir(exist_ok=True)
    
    # we use montgomery reduction to represent element in pmns_params
    k, p, phi_pow = pmns_params['k'], pmns_params['p'], pmns_params['phi_pow']
    gamma = pmns_params['gamma']
    L, E = pmns_params['L'], pmns_params['E']

    phi = 2**phi_pow

    # new elements generation to allow pmns_params representation using Montgomery
    constructed_parameters = pmns_params.keys()

    if 'M' not in constructed_parameters or 'N' not in constructed_parameters :
        M, N = search_m_and_n(k, p, gamma, L, E, phi)
        pmns_params.update({'M': M, 'N':N})

    if 'L_inv' not in constructed_parameters :
        pmns_params.update({'L_inv': -(L.inverse()%phi)})

    K = gamma.parent()
    transition_matrix = gen_transition_matrix(gamma, k)
    
    # choose element for reduction method
    sig = inspect.signature(reduction_method)
    usefull_args = {k: v for k, v in pmns_params.items() if k in sig.parameters}
  
    # create random elements that will be turned into pmns_params representation
    # save those elements and compute reduce product of element
    polynomials_a = []
    polynomials_b = []
    polynomials_reduced = []
    for _ in range(n_test):
        A = K.random_element()
        B = K.random_element()
        
        Pa = convert_element_to_pmns_montgomery(A, transition_matrix, pmns_params)
        Pb = convert_element_to_pmns_montgomery(B, transition_matrix, pmns_params)
        
        P = (Pa * Pb) % E
        R = reduction_method(P, **usefull_args)
        
        polynomials_a.append(Pa.list())
        polynomials_b.append(Pb.list())
        polynomials_reduced.append(R.list())

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("test_reduction_template.j2")
    
    params = {'n': E.degree(), 'n_test': n_test, 'polA_str': fint.format_matrix(polynomials_a), 'polB_str': fint.format_matrix(polynomials_b), 'reds_str': fint.format_matrix(polynomials_reduced)}
    rendered_params = template.render(params)
    
    output_path = test_dir / "test_reduction.h"
    output_path.write_text(rendered_params)


def write_conversion_test(output_dir:str , n_test:int,  pmns_params:dict):
    # Create test_values subdirectory
    test_dir = Path(output_dir) / "test_values"
    test_dir.mkdir(exist_ok=True)
    
    k, p, = pmns_params['k'], pmns_params['p']
    gamma = pmns_params['gamma']
    n = pmns_params['E'].degree()
    rho = pmns_params['rho']
    phi = 2**pmns_params['phi_pow']
    L = pmns_params['L']

    K = gamma.parent()

    elements = []
    for _ in range(n_test):
        element = K.random_element()   
        elements.append([int(c) for c in element._vector_()])

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("test_conversion_template.j2")

    nb_limbs = ceil(p.nbits()/pmns_params['phi_pow'])
    params = {'k': k, 'n_test': n_test, 'elements_decompose': fint.format_matrix_BigInt(elements, nb_limbs), 'nb_limbs': nb_limbs, 'nb_intern_red_calssical': n//k, "phi_pow":pmns_params['phi_pow'], 'rho':rho}
    
    rendered_params = template.render(params)
    output_path = test_dir / "test_conversion.h"
    output_path.write_text(rendered_params)

def write_test(output_dir:str , n_test:int, reduction_method: callable,  pmns_params:dict):
    write_reduction_test(output_dir, n_test, reduction_method,  pmns_params)
    write_conversion_test(output_dir, n_test,  pmns_params)