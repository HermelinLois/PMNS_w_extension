# ==================================
# code_writer.py
# Function to write C reduction code 
# ==================================

from jinja2 import Environment, FileSystemLoader
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = CURRENT_DIR / "templates"
# a donner par orchestrator :UTPU_DIR = CURRENT_DIR.resolve().parents[2]

def write_c_code(output_dir, config):
    
    # acces to templates files to writes specific file
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("general_template.j2")
    
    params = {
        "reduction_file": config['template_file'],
        "method_name": config['c_method_name']
    }

    rendered_code = template.render(params)
    
    output_path = Path(output_dir) / "reduction.c"
    output_path.write_text(rendered_code)
