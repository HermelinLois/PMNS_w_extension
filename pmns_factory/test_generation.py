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
RANGE_TEST = [64, 96, 128, 192, 256, 384, 512, 767, 1024]
TIMEOUT = 60
N_TEST = 1000

def single_test_pmns(args):
    p, k, pmns_module, TIMEOUT = args

    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(TIMEOUT)
    try:
        pmns, round_count = pmns_module.gen_parameters(p, k)
        return {"p": p, "type" : pmns_module.__name__, "status": GOOD,  "norm" : sum(abs(c) for c in pmns['E'].list()), "count": round_count}

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
    #average round pour la taille et average norme pour la taille
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
            
    W = (6, 6, len(f"{N_TEST}/{N_TEST}"), 13, 8)
    header = f"{'SIZE':>{W[0]}} | {'TYPE':<{W[1]}} | {'GOOD':>{W[2]}} | {'ROUND_COUNT':>{W[3]}} | {'NORM':>{W[4]}}"
    sep_thick = '=' * len(header)
    sep_thin  = '-' * len(header)

    f.write(sep_thick + "\n")
    f.write(header + "\n")
    f.write(sep_thick + "\n")

    for idx in range(len(RANGE_TEST)):
        size = RANGE_TEST[idx]
        d0 = datas[pmns_E_type0.__name__][idx]
        d1 = datas[pmns_E_type1.__name__][idx]
        avg = lambda d, i: float(d[i]) / d[STATUS] if d[STATUS] else 0.0
        f.write(f"{size:>{W[0]}} | {'TYPE0':<{W[1]}} | {f'{d0[STATUS]}/{N_TEST}':>{W[2]}} | {avg(d0, ROUND):>{W[3]}.2f} | {avg(d0, NORM):>{W[4]}.2f}\n")
        f.write(f"{size:>{W[0]}} | {'TYPE1':<{W[1]}} | {f'{d1[STATUS]}/{N_TEST}':>{W[2]}} | {avg(d1, ROUND):>{W[3]}.2f} | {avg(d1, NORM):>{W[4]}.2f}\n")
        f.write(sep_thin + "\n")
        
        

if __name__ == "__main__":
    k = 2

    tasks = []
    print("start generation of test parameters")
    for m in RANGE_TEST :
        for _ in range(N_TEST):
            p = random_prime(2**m - 1, lbound=2**(m-1))
            tasks.append((p, k, pmns_E_type0, TIMEOUT))
            tasks.append((p, k, pmns_E_type1, TIMEOUT))
    print("end of the generation")
    
    results = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        results = list(tqdm(executor.map(single_test_pmns, tasks), total=len(tasks), desc="Tests PMNS"))

    s_results = sorted(results, key=lambda x: (x["p"], x["type"]))
    
    with open(f"tests_result{k}.txt", "a") as f:
        #write_all_data(f, s_results)
        write_resume_data(f, results)