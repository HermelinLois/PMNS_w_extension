# include <stdio.h>
# include <stdlib.h>

# include "../codes/reductions_interface.h"
# include "reductions_values.h"

void product(__int128 out[DEGREE], int64_t PolA[DEGREE], int64_t PolB[DEGREE]){
    const int PROD_SIZE = 2*DEGREE;
    __int128 P[PROD_SIZE];

    for (int deg=0; deg<PROD_SIZE; deg++) P[deg]=0;

    for (int degA=0; degA<DEGREE; degA++){
        __int128 ai = PolA[degA];
        for (int degB=0; degB<DEGREE; degB++)
            P[degA + degB] += ai * PolB[degB];
    }

    for (int i=0; i<DEGREE; i++) out[i] = P[i];

    ext_red_w_matrix(__int128, out, P + DEGREE, EXT_MAT);
}


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


int main(){
    __int128 polynomial[DEGREE];
    int64_t out[DEGREE];

    for (int idx=0; idx<N_TEST; idx++){
        for (int i=0; i<DEGREE; i++) out[i] = 0;
        product(polynomial, POL_A[idx], POL_B[idx]);

        reduction_montgomery_int128(out, polynomial, MAT_M, MAT_N);
        check_validity(out, MONTGOMERY_PROD_RED[idx], "Montgomery");

        reduction_babai_int128(out, polynomial, L, L_INV_BABAI);
        check_validity(out, BABAI_PROD_RED[idx], "Babai");
    }
    printf("Montgomery and Babai reductions seems to work with given parameters\n/!\\ this only check if reductions operations give the same result compare to python but doesn't check if coefficients are under rho /!\\ \n");
    return 0;
}