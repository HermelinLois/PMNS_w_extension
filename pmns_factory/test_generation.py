# ==================================================
# test_generation.py
# File to test a PMNS construction with specific E type
# ==================================================

from sage.all import random_prime, randint, GF, gcd
from pmns_core.parameters.params_gestion import search_m_and_n
from pmns_core.operations.convertions_gestion import gen_gamma_base, convert_element_to_pmns_montgomery
from pmns_core.operations.reductions.montgomery_reduction import montgomery_reduction
import pmns_E_type0 as type0
import pmns_E_type1 as type1
from tqdm import tqdm

import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import signal
from cysignals.signals import AlarmInterrupt
import AlarmInterrupt
class MyTimeoutError(Exception):
    pass

def alarm_handler(signum, frame):
    raise MyTimeoutError("over time limit")

GOOD = "good"
BAD = "bad"

def single_test_pmns(args):
    p, k, pmns_module, timeout = args

    signal.alarm(timeout)
    try:
        pmns, round_count = pmns_module.gen_parameters(p, k)
        return {"p": p, "type" : pmns_module.__name__, "status": GOOD,  "norm" : sum(abs(c) for c in pmns['E'].list()), "count": round_count}

    except AlarmInterrupt:
        return {"p": p, "type" : pmns_module.__name__, "status": BAD, "error": "TIMEOUT"}

    except Exception as e:
        return {"p": p, "type" : pmns_module.__name__, "status": BAD, "error": str(e)}

    finally:
        signal.alarm(0)
    


if __name__ == "__main__":
    timeout = 60
    n_test = 20
    k = 2

    tasks = []
    for m in [64, 128, 256]:
        for _ in range(n_test):
            p = random_prime(2**m - 1, lbound=2**(m-1))
            tasks.append((p, k, type0, timeout))
            tasks.append((p, k, type1, timeout))
    
    results = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        results = list(tqdm(executor.map(single_test_pmns, tasks), total=len(tasks), desc="Tests PMNS"))

    s_results = sorted(results, key=lambda x: (x["p"], x["type"]))
    with open(f"tests_result{k}.txt", "a") as f:
        for idx, e in enumerate(s_results):
            status = e["status"]
            f.write("=" * 20 + "\n")

            f.write(f"PRIME : {e['p']}\nTYPE : {e['type']}\nSTATUS : {e['status']}\n")
            if status == GOOD:
                f.write(f"NORM : {e['norm']}\nROUND_COUNT : {e['count']}\n")
            else:
                f.write(f"ERROR : {e['error']}\n")

            f.write("=" * 20 + "\n")
            if idx&1:
                f.write("\n")