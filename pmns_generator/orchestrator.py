from sage.all import random_prime
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from writers.format.container import PMNSContainer
import sys
import argparse

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent

ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)

from writers import params_writer, values_writer
from config import PMNS_CONFIG

OUTPUT_DIR_NAME = "pmns_exec"


def write_pmns_config(output_dir: Path, container: dict) -> None:
    templates_dir = CURRENT_DIR / "writers" / "templates"
    env = Environment(loader=FileSystemLoader(str(templates_dir)))
    template = env.get_template("pmns_config_template.j2")


    params = {'mod': container.get('mod'),
              'p': container.get('p'),
              'k': container .get('k'),
              'E': container.get('E'), 
              'gamma': container.get('gamma'),
              'rho': container.get('rho'),
              'phi_pow': container.get('phi_pow')}
    
    rendered_params = template.render(params)

    output_path = output_dir / "config_pmns"
    output_path.write_text(rendered_params)    


def get_container(args):
    if args.load:
        try:
            return PMNSContainer.load(args.nbits, args.k, args.Etype), True
        except (FileNotFoundError, ValueError) as e:
            width = 50
            head = f"<{ '=' * (width-2) }>"

            print(f"\033[93m{head}\033[0m") 
            print(f"\033[93mExtraction failed : {e}\033[0m".center(width))
            print("\033[93m--- Generating new parameters ---\033[0m".center(width))
            print(f"\033[93m{head}\033[0m")

    p = random_prime(2**args.nbits - 1, lbound=2**(args.nbits-1))
    
    pmns_gen = PMNS_CONFIG[args.Etype]
    pmns_params = pmns_gen.gen_parameters(p, args.k)

    return PMNSContainer(args.Etype, pmns_params), False

def get_args():
    parser = argparse.ArgumentParser("PMNS parameters")
    parser.add_argument("-ntests", type=int, default=100)
    parser.add_argument("-nbits", type=int, default=128)
    parser.add_argument("-k", type=int, default=2)
    parser.add_argument("-Etype", type=int, default=0)
    parser.add_argument("--load", action="store_true", help="load pmns in saves files")
    return parser.parse_args()

def main():
    args = get_args()
    container, was_loaded = get_container(args)

    output_path = Path(OUTPUT_DIR_NAME)
    OUTPUT_DIR = output_path if output_path.is_absolute() else (ROOT_DIR / output_path)
    TEMPLATES_DIR = CURRENT_DIR / "writers" / "templates"
    
    write_pmns_config(OUTPUT_DIR, container)

    params_writer.write_params(TEMPLATES_DIR, OUTPUT_DIR, args.ntests, container)
    values_writer.write_values(TEMPLATES_DIR, OUTPUT_DIR, args.ntests, container)
    
    if not was_loaded :
        container.save()

if __name__ == "__main__":
    main()