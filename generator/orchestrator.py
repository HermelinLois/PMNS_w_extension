from writers import c_writer

from sage.all import *
from pathlib import Path
import sys

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent
OUTPUT_DIR = ROOT_DIR / "generated_code"

ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)

from config import *

def gen_parameters_for_method(method:int, **kwargs) -> None:
    pass


def write_pmns_data(n_test:int, m:int, k:int, Etype:int, method:int) -> None:
    assert Etype in PMNS_CONFIG.keys()
    assert method in REDUCTION_CONFIG.keys()

    OUTPUT_DIR.mkdir(exist_ok=True)

    p = random_prime(2**m - 1, lbound=2**(m-1))

    PMNS = PMNS_CONFIG[Etype]
    config = REDUCTION_CONFIG[method]
    
    params = PMNS.gen_parameters(p, k)
    print(params)
    
    c_writer.write_c_code(OUTPUT_DIR, config)
    
write_pmns_data(0, 64, 2, E_TYPE0, METHOD_MONTGOMERY)