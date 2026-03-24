from sage.all import randint
from jinja2 import Environment, FileSystemLoader
from pmns_factory.core.parameters.params_gestion import search_m_and_n
from pmns_factory.core.operations.convertions_gestion import convert_element_to_pmns_montgomery, gen_transition_matrix
import pmns_generator.writers.format_to_c.int_to_c as fint
from pathlib import Path
import inspect

CURRENT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = CURRENT_DIR / "templates/tests_templates"


def write_reduction_test(output_dir:str , n_test:int, reduction_method: callable,  **pmns_params:dict) -> None:
    # we use montgomery reduction to represent element in PMNS
    k, p, phi_pow = pmns_params['k'], pmns_params['p'], pmns_params['phi_pow']
    gamma = pmns_params['gamma']
    L, E = pmns_params['L'], pmns_params['E']

    phi = 2**phi_pow

    # new elements generation to allow pmns representation using Montgomery
    constructed_parameters = pmns_params.keys()
    if 'M' not in constructed_parameters or 'N' not in constructed_parameters:
        M, N = search_m_and_n(k, p, gamma, L, E, phi)
        pmns_params.update({'M': M, 'N':N})

    K = gamma.parent()
    transition_matrix = gen_transition_matrix(gamma, k)
    
    # choose element for reduction method
    sig = inspect.signature(reduction_method)
    usefull_args = {k: v for k, v in pmns_params.items() if k in sig.parameters}
  
    # create random elements that will be turned into PMNS representation
    # save those elements and compute reduce product of element
    polynomials_a = []
    polynomials_b = []
    polynomials_reduced = []
    for _ in range(n_test):
        A = K([randint(0, p) for _ in range(k)])
        B = K([randint(0, p) for _ in range(k)])
        
        Pa = convert_element_to_pmns_montgomery(A, transition_matrix, **pmns_params)
        Pb = convert_element_to_pmns_montgomery(B, transition_matrix, **pmns_params)
        
        P = (Pa * Pb) % E
        R = reduction_method(P, **usefull_args)
        
        polynomials_a.append(Pa.list())
        polynomials_b.append(Pb.list())
        polynomials_reduced.append(R.list())

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("test_reduction_template.j2")
    
    params = {'n': E.degree(), 'n_test': n_test, 'polA_str': fint.format_matrix(polynomials_a), 'polB_str': fint.format_matrix(polynomials_b), 'reds_str': fint.format_matrix(polynomials_reduced)}
    rendered_params = template.render(params)
    
    output_path = Path(output_dir) / "test_reduction.h"
    output_path.write_text(rendered_params)


def write_conversion_test(output_dir:str , n_test:int,  **pmns_params):
    k, p, = pmns_params['k'], pmns_params['p']
    gamma = pmns_params['gamma']

    K = gamma.parent()
    fname = "INT128"

    elements = []
    for _ in range(n_test):
        element = [randint(0, p) for _ in range(k)]
        elements.append(element)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("test_conversion_template.j2")
    
    params = {'k': k, 'n_test': n_test,'fname': fname, 'elements_str': fint.format_matrix128(elements, fname=fname)}
    rendered_params = template.render(params)
    
    output_path = Path(output_dir) / "test_conversion.h"
    output_path.write_text(rendered_params)

def write_test(output_dir:str , n_test:int, reduction_method: callable,  **pmns_params:dict):
    write_reduction_test(output_dir, n_test, reduction_method,  **pmns_params)
    write_conversion_test(output_dir, n_test,  **pmns_params)