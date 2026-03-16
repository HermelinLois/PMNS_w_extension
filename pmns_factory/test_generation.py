# ==================================================
# test_generation.py
# File to test a PMNS construction with specific E type
# ==================================================

from sage.all import random_prime
import pmns_E_type0
import pmns_E_type1
from tqdm import tqdm

import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import signal
from cysignals.signals import AlarmInterrupt


class MyTimeoutError(Exception):
    pass

def alarm_handler(signum, frame):
    raise MyTimeoutError("over time limit")

GOOD = "good"
BAD = "bad"
RANGE_TEST = [256, 512]
TIMEOUT = 120
N_TEST = 2000

def single_test_pmns(args):
    p, k, pmns_module, TIMEOUT = args

    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(TIMEOUT)
    try:
        pmns, round_count = pmns_module.gen_parameters(p, k)
        return {"p": p, "type" : pmns_module.__name__, "status": GOOD,  "norm" : sum(abs(c) for c in pmns['E'].list()), "count": round_count, "det": pmns['L'].det()}

    except (AlarmInterrupt, MyTimeoutError):
        return {"p": p, "type" : pmns_module.__name__, "status": BAD, "error": "TIMEOUT"}

    except Exception as e:
        return {"p": p, "type" : pmns_module.__name__, "status": BAD, "error": str(e)}

    finally:
        signal.alarm(0)
    
def write_all_data(f, results):
    for idx, e in enumerate(results):
        status = e["status"]
        f.write("=" * 20 + "\n")
        f.write(f"PRIME : {e['p']}\nTYPE : {e['type']}\nSTATUS : {e['status']}\n")
        f.write(f"NORM : {e['norm']}\nROUND_COUNT : {e['count']}\n") if status == GOOD else f.write(f"ERROR : {e['error']}\n")
        f.write("=" * 20 + "\n")
            
        if idx&1:
            f.write("\n")
            
def write_resume_data(f, results):
    STATUS = 0
    ROUND = 1
    NORM = 2
    datas = {
        pmns_E_type0.__name__: [[0, 0, 0] for _ in range(len(RANGE_TEST))],
        pmns_E_type1.__name__: [[0, 0, 0] for _ in range(len(RANGE_TEST))],
    }
    
    for e in results:
        data = datas[e['type']]
        m = e['p'].nbits()
        idx = RANGE_TEST.index(m)
        
        status = e['status']
        if status == GOOD:
            data[idx][STATUS] += 1
            data[idx][ROUND] += e['count']
            data[idx][NORM] += e['norm']
    

    f.write("#{TIMEOUT = }")
    W = (10, 10, len(f"{N_TEST}/{N_TEST}"), 10, 10, 10)
    header = f"{'SIZE':>{W[0]}} | {'TYPE':>{W[1]}} | {'GOOD':>{W[2]}} | {'ROUND_COUNT':>{W[3]}} | {'NORM':>{W[4]}} | {'DET':>{W[5]}}"
    sep_thick = '=' * len(header)
    sep_thin  = '-' * len(header)

    f.write(sep_thick + "\n")
    f.write(header + "\n")
    f.write(sep_thick + "\n")

    for idx in range(len(RANGE_TEST)):
        size = RANGE_TEST[idx]
        type0 = datas[pmns_E_type0.__name__][idx]
        type1 = datas[pmns_E_type1.__name__][idx]
        avg = lambda d, i: float(d[i]) / d[STATUS] if d[STATUS] else 0.0
        f.write(f"{size:>{W[0]}} | {'TYPE0':<{W[1]}} | {f'{type0[STATUS]}/{N_TEST}':>{W[2]}} | {avg(type0, ROUND):>{W[3]}.2f} | {avg(type0, NORM):>{W[4]}.2f} | {avg(type0, DET)}|\n")
        f.write(f"{size:>{W[0]}} | {'TYPE1':<{W[1]}} | {f'{type1[STATUS]}/{N_TEST}':>{W[2]}} | {avg(type1, ROUND):>{W[3]}.2f} | {avg(type1, NORM):>{W[4]}.2f} | {avg(type0, DET)}\n")
        f.write(sep_thin + "\n")
        

def gen_prime(m):
    return random_prime(2**m - 1, lbound=2**(m-1))

if __name__ == "__main__":
    k = 2

    primes_sizes = []
    for m in RANGE_TEST:
        primes_sizes.extend([m] * N_TEST)

    primes = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor :
        primes = list(tqdm(executor.map(gen_prime, primes_sizes), total=len(RANGE_TEST) * N_TEST, desc="Prime generation"))
    
    tasks = []
    for p in primes:
        tasks.append((p, k, pmns_E_type0, TIMEOUT))
        tasks.append((p, k, pmns_E_type1, TIMEOUT))

    results = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        results = list(tqdm(executor.map(single_test_pmns, tasks), total=len(tasks), desc="PMNS generation"))

    s_results = sorted(results, key=lambda x: (x["p"], x["type"]))
    
    with open(f"tests_result{k}.txt", "a") as f:
        #write_all_data(f, s_results)
        write_resume_data(f, results)