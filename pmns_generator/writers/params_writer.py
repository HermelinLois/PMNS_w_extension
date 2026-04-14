# ==================================
# params_writer.py
# Functions to write parameters in C 
# and in txt format
# ==================================

from sage.all import matrix
from jinja2 import Environment, FileSystemLoader
from .format_to_c.int_to_c import format_matrix
from pathlib import Path
import sys

CURRENT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = CURRENT_DIR / "templates"
ROOT_DIR = CURRENT_DIR.parent.parent
ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)
    
from config import METHOD_MONTGOMERY, METHOD_BABAI
from pmns_factory.core.parameters.matrix_gestion import gen_overflow_matrix

def compute_additional_params(method, params:dict) -> None:
    E = params['E']
    L = params['L']
    n = E.degree()
    phi_pow = params['phi_pow']

    external_reduction_matrix = gen_overflow_matrix(E)
    zero_row = matrix(external_reduction_matrix.base_ring(), 1, external_reduction_matrix.ncols())

    params.update({'ext_red_mat_str': format_matrix(external_reduction_matrix.stack(zero_row))})

    if method == METHOD_MONTGOMERY:
        from pmns_factory.core.operations.reductions.montgomery_reduction import gen_mn_reduction_matrix
        from pmns_factory.core.operations.reductions.montgomery_reduction import search_m_and_n
        
        phi = 2**phi_pow
        k = params['k']
        p = params['p']
        gamma = params['gamma']
        
        M, N = search_m_and_n(k, p, gamma, L, E, phi)
        mat_m, mat_n = gen_mn_reduction_matrix(M, E, phi)

        params['phi'] = phi
        params.update({'n': n, 'mat_m_str': format_matrix(mat_m), 'mat_n_str': format_matrix(mat_n), 'M': M, 'N': N})
        return
        
    if method == METHOD_BABAI:
        from pmns_factory.core.operations.reductions.babai_reduction import gen_params_for_babai
        
        rho = params['rho']

        h1, h2, l_inv_babai = gen_params_for_babai(L, phi_pow, rho, E)

        params.update({'n': n, 'h1': h1, 'h2': h2, 'L_inv_babai': l_inv_babai, 'L_inv_babai_str': format_matrix(l_inv_babai), 'L_str': format_matrix(L)})
        return


def write_c_params(env, output_dir:str, method, pmns_params:dict) -> None:
    template = env.get_template("params_c_template.j2")

    is_method_montgomery = (method == METHOD_MONTGOMERY)

    compute_additional_params(method, pmns_params)
    used_params = {"is_method_montgomery" : is_method_montgomery, **pmns_params}
    
    rendered_params = template.render(used_params)
    
    output_path = Path(output_dir) / "param.h"
    output_path.write_text(rendered_params)


def write_params(output_dir, method, pmns_params):
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR / "add_parameters_templates")))
    write_c_params(env, output_dir, method, pmns_params)