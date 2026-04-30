# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include "../codes/conversions_interface.h"
# include "conversions_values.h"
# include "../codes/measurement_utils.c"


void print_pol(int64_t *P){
    for (int idx=0; idx<DEGREE; idx++)
        printf("%lld  ", (long long)P[idx]);
    printf("\n");
}



void check_validity(int64_t P[DEGREE], int64_t C[DEGREE], const char* method){
    for (int i=0; i<DEGREE; i++){
        if(C[i] != P[i]){
            printf("/!\\ Error has been found /!\\ \n");
            printf("Generate with Python :\n");
            print_pol(C);
            printf("Generate with C (with %s) :\n", method);
            print_pol(P);
            exit(1);
        }
    }
}



static inline void reset_polynomial(int64_t polynomial[DEGREE]){
    for (int i=0; i<DEGREE; i++) 
        polynomial[i]=0;
}



void test_equality(){
    int64_t polynomial[DEGREE];

    for (int idx=0; idx<N_TESTS; idx++){
        convert_element_to_pmns_exact(polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_validity(polynomial, CONVERTED_ELEMENTS_EXACT[idx], "Exact");

        reset_polynomial(polynomial);

        convert_element_to_pmns_pseudo_fast(polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_validity(polynomial, CONVERTED_ELEMENTS_PSEUDO_FAST[idx], "Pseudo-Fast");

        reset_polynomial(polynomial);

        convert_element_to_pmns_fast(polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_validity(polynomial, CONVERTED_ELEMENTS_FAST[idx], "Fast");
    }
    printf("Montgomery conversions seems to work with given parameters\n");
}



void rand_int(mp_limb_t out[EXTENSION_DEGREE][N_LIMBS]){
	mpz_t p_mpz, rand_val;
    mpz_init(p_mpz);
    mpz_init(rand_val);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));


	mpz_import(p_mpz, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);

    for (int i=0; i<EXTENSION_DEGREE; i++){
        mpz_urandomm(rand_val, state, p_mpz);
        mpn_copyi(out[i], mpz_limbs_read(rand_val), N_LIMBS);
    }

    gmp_randclear(state);
    mpz_clear(p_mpz);
    mpz_clear(rand_val);
}



void do_bench(void (*to_pmns)(int64_t pmns[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]), char* method_name){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

	mp_limb_t a[EXTENSION_DEGREE][N_LIMBS];
    int64_t polynomial[DEGREE];
	
	for(int i=0;i<N_BENCH_TESTS;i++){
		rand_int(a);
		to_pmns(polynomial, a);
	}
	
	for(int i=0;i<N_BENCH_SAMPLES;i++){
		rand_int(a);

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
            reset_polynomial(polynomial);

			t1 = cpucyclesStart();
			to_pmns(polynomial, a);
			t2 = cpucyclesStop();

			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2 + diff_t + 1;
			} else { diff_t = t2-t1; }

			if(timermin > diff_t) 
                timermin = diff_t;
			else 
                if(timermax < diff_t) 
                    timermax = diff_t;
			cycles[j]=diff_t;
		}
        
		meanTimermin += timermin;
		meanTimermax += timermax;

		statTimer = quartiles(cycles,N_BENCH_TESTS);
		
		medianTimer += statTimer[1];
		free(statTimer);
	}
	
	free(cycles);

    printf("====================================================\n");
    printf("|%*s%*s|\n", 25 + (int)(strlen(method_name)/2), method_name, 25 - (int)(strlen(method_name)/2), "");
    printf("====================================================\n");
    printf("| %-20s : %-26llu|\n", "Minimum cycles", (unsigned long long)(meanTimermin/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Maximum cycles", (unsigned long long)(meanTimermax/N_BENCH_SAMPLES));
    printf("====================================================\n");
}

void test_speed(){
    do_bench(convert_element_to_pmns_exact, "Exact");
    do_bench(convert_element_to_pmns_pseudo_fast, "Pseudo-Fast");
    do_bench(convert_element_to_pmns_fast, "Fast");
}

int main(){
    test_equality();
    test_speed();
    
    return 0;
}