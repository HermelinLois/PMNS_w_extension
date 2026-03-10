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

import signal
from cysignals.signals import AlarmInterrupt

class MyTimeoutError(Exception):
    pass

def alarm_handler(signum, frame):
    raise MyTimeoutError("over time limit")

signal.signal(signal.SIGALRM, alarm_handler)




if __name__ == "__main__":
    time_out = 10
    n_test = 10
    k = 5

    for m in [64, 128, 256]:
        for _ in range(n_test):
            p = random_prime(2**m - 1, lbound=2**(m-1))

            print(gcd(p-1, (p**k - 1)/gcd(k, p**k-1))>1)
            type0_is_ok, type1_is_ok = True, True
            count0, count1 = -1, -1 
            error_type0, error_type1 = None, None

            signal.alarm(time_out)
            try:
                _, count1 = type1.gen_parameters(p, k)
            except (MyTimeoutError, Exception) as e:
                type1_is_ok = False
                error_type1 = "Timeout" if "time limit" in str(e) else str(e)
            finally:
                signal.alarm(0)
                

            signal.alarm(time_out)
            try:
                _, count0 = type0.gen_parameters(p, k)
            except (MyTimeoutError, Exception) as e:
                type0_is_ok = False
                error_type0 = "Timeout" if "time limit" in str(e) else str(e)
            finally:
                signal.alarm(0)
            
            with open(f"recap_{k}.txt", 'a') as f:
                f.write(f"{m}\n{p}\n{k}\n({error_type0}, {error_type1})\n({count0}, {count1})\n\n")