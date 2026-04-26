# include "reductions_interface.h"
# include "../params/conversions_params.h"
# include <stdio.h>
# include <inttypes.h>















// afficher un mpn signé (N_LIMBS limbs)
static inline void print_mpn(const mp_limb_t x[N_LIMBS]) {
    mpz_t z;
    mpz_init(z);
    mpz_import(z, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, x);
    // sign-extend : si le bit haut est à 1, soustraire 2^(N_LIMBS*64)
    if ((x[N_LIMBS-1] >> (GMP_NUMB_BITS-1)) & 1) {
        mpz_t mod;
        mpz_init_set_ui(mod, 1);
        mpz_mul_2exp(mod, mod, (mp_bitcnt_t)(N_LIMBS * GMP_NUMB_BITS));
        mpz_sub(z, z, mod);
        mpz_clear(mod);
    }
    gmp_printf("%Zd", z);
    mpz_clear(z);
}





















// init mp_limb_element with element data and multiply the values by phi^power mod p
static inline void init_element_times_phipow(mp_limb_t out[DEGREE][N_LIMBS], const mp_limb_t data[EXTENSION_DEGREE][N_LIMBS], unsigned int power) {
    int total_limbs = N_LIMBS + power;
    
    mp_limb_t temp_el[total_limbs];
    mp_limb_t quotient[power + 1];
    mp_limb_t remainder[N_LIMBS];

    for (int deg = 0; deg < EXTENSION_DEGREE; deg++) {
        mpn_zero(temp_el, power);
        for(int i = 0; i < N_LIMBS; i++) temp_el[power + i] = data[deg][i];


        // apply modulus p
        mpn_tdiv_qr(quotient, remainder, 0, temp_el, total_limbs, P, N_LIMBS);
        mpn_copyi(out[deg], remainder, N_LIMBS);
    }

    for (int deg = EXTENSION_DEGREE; deg < DEGREE; deg++)
        mpn_zero(out[deg], N_LIMBS);
}


// apply change basis operation to represent element in basis gamma
static inline void element_change_basis(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t extension_element[DEGREE][N_LIMBS]){
    int ACC_SIZE = 2 * N_LIMBS + 3;

    mp_limb_t quotient[ACC_SIZE - N_LIMBS + 1];
    mp_limb_t remainder[N_LIMBS];
    mp_limb_t acc[ACC_SIZE];
    mp_limb_t prod[ACC_SIZE];

    for (int j = 0; j < EXTENSION_DEGREE; j++){
        mpn_zero(acc, ACC_SIZE);
        mpn_zero(prod, ACC_SIZE);

        for (int i = 0; i < EXTENSION_DEGREE; i++){
            mpn_zero(prod, ACC_SIZE);
            mpn_mul_n(prod, extension_element[i], TRANSITION_MATRIX[i][j], N_LIMBS);
            (void) mpn_add_n(acc, acc, prod, ACC_SIZE);
        }

        // Reduce accumulated result modulo P
        mpn_tdiv_qr(quotient, remainder, 0, acc, ACC_SIZE, P, N_LIMBS);
        mpn_copyi(out[j], remainder, N_LIMBS);
    }

    for (int j = EXTENSION_DEGREE; j < DEGREE; j++)
        mpn_zero(out[j], N_LIMBS);
}




static inline void addmul_polynomial_int64(__int128 out[DEGREE], const int64_t polynomial[DEGREE], int64_t coef){
    if (coef == 0) return;

    for (int deg=0; deg<DEGREE; deg++)
        out[deg] += (__int128)polynomial[deg] * coef;
}


static inline uint64_t apply_mask(mp_limb_t element_data[N_LIMBS], uint64_t mask_size, size_t *mask_start_pos) {
    uint64_t start = *mask_start_pos;
    uint64_t limb_idx = start / GMP_NUMB_BITS;

    const uint64_t mask = (1ULL << mask_size) - 1ULL;
    uint64_t shift = start % GMP_NUMB_BITS;
    uint64_t chunk = 0;

    // apply mask which is lower than phi_po
    if (limb_idx < N_LIMBS) {
        chunk = (uint64_t)(element_data[limb_idx] >> shift);

        // apply a second mask on higher register in case the mask apply on 2 register at the same time
        if (mask_size > (GMP_NUMB_BITS - shift) && (limb_idx + 1 < N_LIMBS)) {
            chunk |= (uint64_t)(element_data[limb_idx + 1] << (GMP_NUMB_BITS - shift));
        }
    }
    // increase starting point using given mask_size
    *mask_start_pos += mask_size;
    return chunk & mask;
}


void convert_element_to_pmns_fast(int64_t out[DEGREE], mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    __int128 polynomial[DEGREE] = {0};    

    for (int deg=0; deg<EXTENSION_DEGREE; deg++){
        size_t mask_pos = 0;

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_polynomial_int64(polynomial, PMNS_THETA_FAST[deg][i], part);
        }
    }
    reduction_montgomery_int128(out, polynomial, L, L_INV);
}   



// convert element using exact conversion
void convert_element_to_pmns_classical(int64_t out[DEGREE], mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    mp_limb_t vector[DEGREE][N_LIMBS];
    mp_limb_t vector2[DEGREE][N_LIMBS];

    init_element_times_phipow(vector, element_data, N_INT_RED_CLASSICAL);
    element_change_basis(vector, vector);
    
    int sign[DEGREE] = {0};
    for (int it = 0; it < N_INT_RED_CLASSICAL; it++)
        reduction_montgomery_mpn_test(vector, vector, sign, L, L_INV, it);

    for (int i=0; i<DEGREE; i++){
        out[i] = sign[i] ? -(int64_t) vector[i][0]:  (int64_t) vector[i][0];
    }

    // init_element_times_phipow(vector2, element_data, N_INT_RED_CLASSICAL);
    // element_change_basis(vector2, vector2);
    // for (int it = 0; it < N_INT_RED_CLASSICAL; it++)
    //     reduction_montgomery_mpn(vector2, vector2, L, L_INV, it);

    // for (int i=0; i<DEGREE; i++){
    //     out[i] = (int64_t)vector2[i][0];
    // }
}


// void convert_element_to_pmns_pseudo_fast(int64_t out[DEGREE], mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
//     mp_limb_t vector[DEGREE][N_INT_RED_PSEUDO_FAST + 1];

//     for (int deg=0; deg<EXTENSION_DEGREE; deg++){
//         size_t mask_pos = 0;
//         for (int i=0; i<N_POL; i++)
//             uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
//     }
// }
