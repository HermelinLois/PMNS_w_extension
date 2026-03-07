from writers import code_writer, params_writer, test_writer
from sage.all import random_prime
from pathlib import Path
import sys


CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent
OUTPUT_DIR = ROOT_DIR / "generated_code"

ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)

from config import PMNS_CONFIG, REDUCTION_CONFIG, E_TYPE0, METHOD_MONTGOMERY, METHOD_BABAI


def write_pmns_data(n_test:int, m:int, k:int, Etype:int, method:int) -> None:
    assert Etype in PMNS_CONFIG.keys()
    assert method in REDUCTION_CONFIG.keys()

    OUTPUT_DIR.mkdir(exist_ok=True)

    p = random_prime(2**m - 1, lbound=2**(m-1))

    PMNS = PMNS_CONFIG[Etype]
    config = REDUCTION_CONFIG[method]
    
    pmns_params = PMNS.gen_parameters(p, k)
    
    code_writer.write_c_code(OUTPUT_DIR, config)
    params_writer.write_params(OUTPUT_DIR, method, pmns_params)
    test_writer.write_test(OUTPUT_DIR, n_test, config['py_func'],  **pmns_params)
    
write_pmns_data(100, 128, 2, E_TYPE0, METHOD_MONTGOMERY)