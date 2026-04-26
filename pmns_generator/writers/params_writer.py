from pmns_factory.core.operations.reductions.babai_reduction import gen_params_for_babai
from pmns_factory.core.operations.reductions.montgomery_reduction import gen_mn_reduction_matrix, search_m_and_n
from pmns_factory.core.operations.convertions_gestion import gen_transition_matrix, compute_nb_internal_reductions, montgomery_exact_conversion
from pmns_factory.core.parameters.matrix_gestion import gen_overflow_matrix
from pmns_generator.writers.format.container import PMNSContainer
from .format import format_element as format
from jinja2 import Environment, FileSystemLoader
from sage.all import ceil


def write_pmns_params(env, output_dir, n_test, container): 
    params = {'n': container.get('n'), 
              'rho': container.get('rho'),
              'n_limbs': container.get('n_limbs'),
              'k': container.get('k'),
              'n_test' : n_test,
              'gamma': container.get('gamma', True),
              'p': container.get('p', True)}
    
    template = env.get_template("pmns_params_template.j2")
    rendered_params = template.render(params)
    
    output_path = output_dir / "pmns_params.h"
    output_path.write_text(rendered_params)



def write_reductions_params(env, output_dir, container):    
    template = env.get_template("reductions_params_template.j2")
    
    params = {"h1": container.get("h1"),
              "h2": container.get('h2'),
              "ext_red": container.get('ext_red_mat', True),
              "M_mat": container.get('M_mat', True),
              "N_mat": container.get('N_mat', True),
              'L': container.get('L', True),
              'L_inv': container.get('L_inv', True),
              'L_inv_babai': container.get('L_inv_babai', True)}
    
    rendered_params = template.render(params)
    
    output_path = output_dir / "reductions_params.h"
    output_path.write_text(rendered_params)



def write_conversions_params(env, output_dir, container):    
    template = env.get_template("conversions_params_template.j2")
    
    params = {"n_red_exact": container.get('n_red_exact'),
              "n_red_pseudo": container.get('n_red_pseudo'),
              "n_red_fast": container.get('n_red_fast'),
              "theta_pow": container.get('theta_pow'),
              "n_pol": container.get('n_pol'),
              "T_mat": container.get('T_mat', True),
              "int_pols": container.get('int_pols', True),
              "z_pols": container.get('z_pols', True),
              "fast_pols": container.get('fast_pols', True)}

    
    rendered_params = template.render(params)
    output_path = output_dir / "conversions_params.h"
    output_path.write_text(rendered_params)
    
    

def write_params(templates_dir, output_dir:str, n_test, container:PMNSContainer):
    PARAMS_DIR =  templates_dir / "params"
    PARAMS_OUTPUT = output_dir / "params"
    
    PARAMS_OUTPUT.mkdir(exist_ok=True)
    env = Environment(loader=FileSystemLoader(str(PARAMS_DIR)))
        
    write_pmns_params(env, PARAMS_OUTPUT, n_test, container)
    write_reductions_params(env, PARAMS_OUTPUT, container)
    
    container.add_conversions_parameters()
    
    write_conversions_params(env, PARAMS_OUTPUT, container)
    