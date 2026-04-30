from jinja2 import Environment, FileSystemLoader
from pmns_factory.core.operations.convertions_gestion import gen_transition_matrix, montgomery_fast_conversion, montgomery_exact_conversion, montgomery_pseudo_fast_conversion
from pmns_factory.core.operations.reductions.babai_reduction import babai_rounding_limited_reduction
from pmns_factory.core.operations.reductions.montgomery_reduction import fast_montgomery_reduction
from pmns_generator.writers.format.container import PMNSContainer
from .format import format_element as format
from random import choice
from sage.all import PolynomialRing, ZZ

def write_conversions_values(env, output_dir:str, n_test:int, container:PMNSContainer) -> list: 
    gamma = container.get('gamma')
    K = gamma.parent()

    elements = []
    conversions_exact = []
    conversions_fast = []
    conversions_pseudo_fast = []

    for _ in range(n_test):
        element = K.random_element()   
        elements.append([int(c) for c in element._vector_()])

        exact_poly = montgomery_exact_conversion(element, container)
        conversions_exact.append(exact_poly)

        pseudo_fast_poly = montgomery_pseudo_fast_conversion(element, container)
        conversions_pseudo_fast.append(pseudo_fast_poly)

        fast_poly = montgomery_fast_conversion(element, container)
        conversions_fast.append(fast_poly)
        
    tests_params = {'elements_mpn': format.format_matrix_to_mpn(elements, container.get('n_limbs')), 
                    'conversions_exact': format.format_matrix_to_int64(conversions_exact),
                    'conversions_fast': format.format_matrix_to_int64(conversions_fast),
                    'conversions_pseudo_fast': format.format_matrix_to_int64(conversions_pseudo_fast)}
        
    template = env.get_template("conversions_values_template.j2")
    rendered_params = template.render(tests_params)
    
    output_path = output_dir / "conversions_values.h"
    output_path.write_text(rendered_params)
    
    return conversions_exact



def write_reduction_values(env, output_dir, n_test, container, convs_pool) -> None:  
    E = container.get('E')
    PR = PolynomialRing(ZZ,"X")
    L = container.get('L_origin')
    L_inv = container.get('L_inv_origin')
    
    polynomials_a = []
    polynomials_b = []
    montgomery_reductions = []
    babai_reductions = []
    for _ in range(n_test):
        
        pol_A = choice(convs_pool)
        pol_B = choice(convs_pool)
            
        prod = (PR(pol_A) * PR(pol_B)) % E    
        polynomials_a.append(pol_A)
        polynomials_b.append(pol_B)
        
        montgomery_red = fast_montgomery_reduction(prod, L, L_inv)
        babai_red = babai_rounding_limited_reduction(prod, container)
            
        montgomery_reductions.append(montgomery_red)
        babai_reductions.append(babai_red)
            
    tests_params = {'polA': format.format_matrix_to_int64(polynomials_a), 
                    'polB': format.format_matrix_to_int64(polynomials_b), 
                    'montgomery_red': format.format_matrix_to_int64(montgomery_reductions), 
                    'babai_red': format.format_matrix_to_int64(babai_reductions)}
    
    template = env.get_template("reductions_values_template.j2")
    rendered_params = template.render(tests_params)
    
    output_path = output_dir / "reductions_values.h"
    output_path.write_text(rendered_params)


def write_values(templates_dir, output_dir:str, n_test:int,  container:dict):
    VALUES_DIR =  templates_dir / "values"
    VALUES_OUTPUT = output_dir / "tests"
    VALUES_OUTPUT.mkdir(exist_ok=True)
    env = Environment(loader=FileSystemLoader(str(VALUES_DIR)))

    convs_pool = write_conversions_values(env, VALUES_OUTPUT, n_test, container)
    write_reduction_values(env, VALUES_OUTPUT, n_test, container, convs_pool)