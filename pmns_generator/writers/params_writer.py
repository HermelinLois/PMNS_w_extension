# ==================================
# params_writer.py
# Functions to write parameters in C 
# and in txt format
# ==================================

from sage.all import PolynomialRing,ZZ
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

PR = PolynomialRing(ZZ,"X")

def compute_additional_params(method, params:dict) -> None:
    E = params['E']
    L = params['L']
    n = E.degree()
    phi_pow = params['phi_pow']
    
    if method == METHOD_MONTGOMERY:
        from pmns_factory.pmns_core.operations.reductions.montgomery_reduction import gen_external_reduction_matrix
        from pmns_factory.pmns_core.parameters.params_gestion import search_polynomial_m
        
        phi = 2**phi_pow
        k = params['k']
        p = params['p']
        gamma = params['gamma']
        
        M = search_polynomial_m(L, k, p, gamma, E)
        mat_m, mat_n = gen_external_reduction_matrix(M, E, phi)
        return {'n': n, 'mat_m': format_matrix(mat_m), 'mat_n': format_matrix(mat_n)}
        
    if method == METHOD_BABAI:
        from pmns_factory.pmns_core.operations.reductions.babai_reduction import gen_params_for_babai
        
        rho = params['rho']
        h1, h2, l_inv_babai = gen_params_for_babai(L, phi_pow, rho, E)
        return {'n': n, 'h1': h1, 'h2': h2, 'L_inv_babai': format_matrix(l_inv_babai), 'L': format_matrix(L)}


def write_c_params(output_dir:str, method, pmns_params:dict) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("params_c_template.j2")

    is_method_montgomery = (method == METHOD_MONTGOMERY)
    params = {"is_method_montgomery" : is_method_montgomery, **compute_additional_params(method, pmns_params)}
    
    rendered_params = template.render(params)
    
    output_path = Path(output_dir) / "param.h"
    output_path.write_text(rendered_params)


def write_txt_params(output_dir:str, params:dict) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("params_txt_template.j2")
    
    mod = PR(params['mod'])
    params['mod'] = mod
    rendered_params = template.render(params)
    
    output_path = Path(output_dir) / "pmns_parameters.txt"
    output_path.write_text(rendered_params)


def write_params(output_dir, method, pmns_params):
    write_txt_params(output_dir, pmns_params)
    write_c_params(output_dir, method, pmns_params)