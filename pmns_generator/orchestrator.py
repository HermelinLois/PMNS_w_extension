from writers import code_writer, params_writer, tests_writer
from sage.all import random_prime
from pathlib import Path
import sys
import argparse


CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent

ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)

from config import PMNS_CONFIG, REDUCTION_CONFIG

def write_pmns_data(n_test:int, m:int, k:int, Etype:int, method:int, name:str) -> None:
    assert Etype in PMNS_CONFIG.keys()
    assert method in REDUCTION_CONFIG.keys()
    
    OUTPUT_DIR = ROOT_DIR / name
    OUTPUT_DIR.mkdir(exist_ok=True)

    p = random_prime(2**m - 1, lbound=2**(m-1))

    PMNS = PMNS_CONFIG[Etype]
    config = REDUCTION_CONFIG[method]
    
    pmns_params = PMNS.gen_parameters(p, k)
    
    code_writer.write_c_code(OUTPUT_DIR, config)
    params_writer.write_params(OUTPUT_DIR, method, pmns_params)
    tests_writer.write_test(OUTPUT_DIR, n_test, config['py_func'],  **pmns_params)
    

def write_all():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ntest", type=int, default=1000)
    parser.add_argument("-nbits", type=int, default=128)
    parser.add_argument("-k", type=int, default=2)
    parser.add_argument("-Etype", type=int, default=0)
    parser.add_argument("-method", type=int, default=0)
    parser.add_argument("-name", type=str, default="generated_code")

    args = parser.parse_args()

    write_pmns_data(args.ntest, args.nbits, args.k, args.Etype, args.method, args.name)

if __name__ == "__main__":
    write_all()