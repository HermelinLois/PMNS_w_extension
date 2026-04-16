# ==================================
# code_writer.py
# Function to write C reduction code 
# ==================================

from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from pmns_factory.core.operations.convertions_gestion import gen_transition_matrix
import pmns_generator.writers.format_to_c.int_to_c as fint

CURRENT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = CURRENT_DIR / "templates" / "code_templates"

def write_reduction_code(output_dir, config):
    # Create code subdirectory
    code_dir = Path(output_dir) / "code"
    code_dir.mkdir(exist_ok=True)
    
    # acces to templates files to writes specific file
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("reduction_templates/general_reduction_template.j2")

    params = {
        "reduction_file": f"reduction_templates/{config['template_file']}",
        "method_name": config['c_method_name']
    }

    rendered_code = template.render(params)
    
    output_path = code_dir / "reduction.c"
    output_path.write_text(rendered_code)


def write_conversion_code(output_dir, pmns_params):
    code_dir = Path(output_dir) / "code"
    code_dir.mkdir(exist_ok=True)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("conversion_templates/conversion..j2")

    p = pmns_params['p']
    k = pmns_params['k']
    phi_pow = pmns_params['phi_pow']
    nb_limbs = (p.nbits() + phi_pow - 1) // phi_pow

    transition_matrix = gen_transition_matrix(pmns_params['gamma'], k)
    params = {
        'p_decompose': fint.format_element_to_BigInt(p, nb_limbs),
        'transition_matrix_decompose': fint.format_matrix_BigInt(transition_matrix, nb_limbs)
    }

    rendered_code = template.render(params)

    output_path = code_dir / "conversion.c"
    output_path.write_text(rendered_code)


def write_c_code(output_dir, config, pmns_params):
    write_reduction_code(output_dir, config)
    write_conversion_code(output_dir, pmns_params)