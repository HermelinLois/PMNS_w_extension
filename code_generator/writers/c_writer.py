from jinja2 import Environment, FileSystemLoader
from pathlib import Path

CURRENT_DIR = Path(__file__)
TEMPLATES_DIR = CURRENT_DIR / "templates"
# a donner par orchestrator :UTPU_DIR = CURRENT_DIR.resolve().parents[2]

def write_c_code(root_path, config):
    # acces to templates files to writes specific file
    env = Environment(loader=FileSystemLoader("templates"))
    template = env.get_template("general_template.j2")
    
    params = {
        "reduction_file": config['template_file'],
        "method_name": config['c_method_name']
    }

    rendered_code = template.render(params)
    Path("test_pmns.c").write_text(rendered_code)
