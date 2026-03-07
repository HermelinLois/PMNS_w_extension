from sage.all import randint
from jinja2 import Environment, FileSystemLoader
from pmns_factory.pmns_core.parameters.params_gestion import search_m_and_n
from pmns_factory.pmns_core.operations.convertions_gestion import convert_element_to_pmns_montgomery, gen_gamma_base
import pmns_generator.writers.format_to_c.int_to_c as fint
from pathlib import Path
import sys

CURRENT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = CURRENT_DIR / "templates"
    
def write_test(output_dir:str , n_test:int, reduction_method: callable,  **pmns_params:dict) -> None:
    # we use montgomery reduction to represent element in PMNS
    k, p, phi_pow = pmns_params['k'], pmns_params['p'], pmns_params['phi_pow']
    gamma = pmns_params['gamma']
    L,E = pmns_params['L'], pmns_params['E']
    
    # we search M and N to represent element in PMNS
    M, N = search_m_and_n(k, p, gamma, L, E, phi_pow)
    
    # we extract data gamma
    K = gamma.parent()
    gamma_base = gen_gamma_base(gamma, k)
    
    # we set dict wich will be used in the conversion process
    convert_params =  {'M': M, 'N': N, **pmns_params}
    
    # name use to recompose decompose int128 in C
    fname = "INT128"
    
    # create and stock polynomial product
    Prods = []
    Reds = []
    for _ in range(n_test):
        A = K([randint(0, p) for _ in range(k)])
        B = K([randint(0, p) for _ in range(k)])
        
        Pa = convert_element_to_pmns_montgomery(A, gamma_base, **convert_params)
        Pb = convert_element_to_pmns_montgomery(B, gamma_base, **convert_params)
        
        P = (Pa * Pb) % E
        R = reduction_method(P, M, N, E, 2**phi_pow)
        
        Prods.append(P.list())
        Reds.append(R.list())
        
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("test_template.j2")
    
    params = {'n': E.degree(), 'n_test': n_test,'fname': fname, 'prods': fint.format_matrix128(Prods,fname), 'reds': fint.format_matrix(Reds)}
    rendered_params = template.render(params)
    
    output_path = Path(output_dir) / "test.h"
    output_path.write_text(rendered_params)
