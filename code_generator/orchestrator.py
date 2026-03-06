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

def gen_parameters_for_method(method:int, **kwargs) -> None:
    pass


def write_pmns_data(n_test:int, m:int, k:int, Etype:int, method:int) -> None:
    assert k > 1, "extension degree must be at least 2"
    assert Etype in ETYPES_CONFIG.keys()
    assert method in REDUCTION_CONFIG.keys()

    OUTPUT_DIR.mkdir(exist_ok=True)

    p = random_prime(2**m - 1, lbound=2**(m-1))

    config = ETYPES_CONFIG[Etype]
    params = Econfig.gen_parameters(p,k)  
