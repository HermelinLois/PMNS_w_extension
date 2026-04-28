# include <stdio.h>
# include <stdlib.h>

# include "../codes/conversions_interface.h"
# include "conversions_values.h"


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

int main(){
    int64_t polynomial[DEGREE];

    for (int idx=0; idx<N_TEST; idx++){
        convert_element_to_pmns_exact(polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_validity(polynomial, CONVERTED_ELEMENTS[idx], "Classical");

        reset_polynomial(polynomial);

        //convert_element_to_pmns_pseudo_fast(polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        //check_validity(polynomial, CONVERTED_ELEMENTS[idx], "Pseudo-Fast");

        reset_polynomial(polynomial);

        convert_element_to_pmns_fast(polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_validity(polynomial, CONVERTED_ELEMENTS[idx], "Fast");
    }
    printf("Montgomery conversions seems to work with given parameters\n");
    return 0;
}