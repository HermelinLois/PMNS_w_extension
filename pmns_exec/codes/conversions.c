# include "reductions_interface.h"
# include "../params/conversions_params.h"
# include <stdio.h>

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



// compute the product of the given element with the transition matrix and apply modulus p to the result
static inline void element_change_basis(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t extension_element[DEGREE][N_LIMBS]){
    int PROD_SIZE = 2 * N_LIMBS;
    int ACC_SIZE  = 2 * N_LIMBS + 1;

    mp_limb_t quotient[ACC_SIZE - N_LIMBS + 1];
    mp_limb_t remainder[N_LIMBS];
    mp_limb_t acc[ACC_SIZE];
    mp_limb_t prod[ACC_SIZE];
    
    mp_limb_t tmp[EXTENSION_DEGREE][N_LIMBS];
    for (int i=0; i<EXTENSION_DEGREE; i++)
        mpn_copyi(tmp[i], extension_element[i], N_LIMBS);

    for (int j=0; j<EXTENSION_DEGREE; j++){
        mpn_zero(acc, ACC_SIZE);

        for (int i=0; i<EXTENSION_DEGREE; i++){
            mpn_zero(prod, ACC_SIZE);
            mpn_mul_n(prod, tmp[i], TRANSITION_MATRIX[i][j], N_LIMBS);
            mp_limb_t carry = mpn_add_n(acc, acc, prod, PROD_SIZE);
            acc[PROD_SIZE] += carry;
        }

        mpn_tdiv_qr(quotient, remainder, 0, acc, ACC_SIZE, P, N_LIMBS);
        mpn_copyi(out[j], remainder, N_LIMBS);
    }

    for (int j=EXTENSION_DEGREE; j<DEGREE; j++)
        mpn_zero(out[j], N_LIMBS);
}



// compute the product of an int64_t polynomial with an int64_t coefficient and add the result to an int128 polynomial
static inline void addmul_pol64_int64(__int128 out[DEGREE], const int64_t polynomial[DEGREE], int64_t coef){
    if (coef == 0) return;

    for (int deg=0; deg<DEGREE; deg++)
        out[deg] += (__int128)polynomial[deg] * coef;
}


static void negate_mpn(mp_limb_t *res, const mp_limb_t *src, int n) {
    mpn_com(res, src, n);
    mpn_add_1(res, res, n, 1);
}

// compute the product of an int128 polynomial with an int64_t coefficient and add the result to an mp_limb_t polynomial
// use complement to 2 representation for negative values
static inline void addmul_pol128_pol64_c2(int n_limbs, mp_limb_t out[DEGREE][n_limbs], __int128 pol128[DEGREE], const int64_t pol64[DEGREE]){
    int TMP_LIMBS = n_limbs + 1;

    mp_limb_t tmp_pol[DEGREE][TMP_LIMBS];
    for (int k=0; k<DEGREE; k++) mpn_zero(tmp_pol[k], TMP_LIMBS);

    for (int i=0; i<DEGREE; i++) {
        int c128_sign = (pol128[i] < 0);
        unsigned __int128 abs_c128 = c128_sign ? (unsigned __int128)(-pol128[i]) : (unsigned __int128)(pol128[i]);
        mp_limb_t mpn_c128[2] = { (mp_limb_t)(unsigned __int128)abs_c128, (mp_limb_t)((unsigned __int128)abs_c128 >> GMP_NUMB_BITS) };

        for (int j=0; j<DEGREE; j++) {
            int c64_sign = (pol64[j] < 0);
            uint64_t abs_c64 = c64_sign ? (uint64_t)(-pol64[j]) : (uint64_t)(pol64[j]);

            mp_limb_t tmp_prod[TMP_LIMBS];
            tmp_prod[2] = mpn_mul_1(tmp_prod, mpn_c128, 2, abs_c64);
            int prod_sign = c128_sign ^ c64_sign;

            if (prod_sign) negate_mpn(tmp_prod, tmp_prod, TMP_LIMBS);

            if (i + j < DEGREE) {
                mpn_add_n(tmp_pol[i+j], tmp_pol[i+j], tmp_prod, TMP_LIMBS);

            } else {
                int k = i + j - DEGREE;
                for (int l=0; l<DEGREE; l++) {
                    int64_t coeff_mat = EXT_MAT[k][l];
                    if (coeff_mat == 0) continue;

                    int mat_sign = (coeff_mat < 0);

                    uint64_t abs_coeff_mat = mat_sign ? (uint64_t)(-coeff_mat) : (uint64_t)coeff_mat;           

                    mp_limb_t red_prod[TMP_LIMBS + 1];
                    red_prod[TMP_LIMBS] = mpn_mul_1(red_prod, tmp_prod, TMP_LIMBS, abs_coeff_mat);

                    if (mat_sign) {
                        negate_mpn(red_prod, red_prod, TMP_LIMBS);
                    }
                    mpn_add_n(tmp_pol[l], tmp_pol[l], red_prod, TMP_LIMBS);
                }
            }
        }
    }

    for (int i = 0; i < DEGREE; i++) {
        mpn_add_n(out[i], out[i], tmp_pol[i], n_limbs);
    }
}


// apply a mask of size mask_size on the element data starting from the position given by mask_start_pos and 
// return the result as a uint64_t, then increase the mask_start_pos by mask_size
static inline uint64_t apply_mask(const mp_limb_t element_data[N_LIMBS], uint64_t mask_size, size_t *mask_start_pos) {
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



// convert an element of the extension field given to its representation in PMNS using the classical method
void convert_element_to_pmns_exact(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    mp_limb_t vector[DEGREE][N_LIMBS];
    
    init_element_times_phipow(vector, element_data, N_INT_RED_CLASSICAL);
    element_change_basis(vector, vector);

    for (int it = 0; it < N_INT_RED_CLASSICAL; it++)
        reduction_montgomery_mpn(N_LIMBS, vector, vector, L, L_INV, it);


    for (int i=0; i<DEGREE; i++)
        out[i] = (int64_t)vector[i][0];
}


// convert an element of the extension field given to its representation in PMNS using a fast method
void convert_element_to_pmns_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    __int128 polynomial[DEGREE] = {0};    

    for (int deg=0; deg<EXTENSION_DEGREE; deg++){
        size_t mask_pos = 0;

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64(polynomial, PMNS_THETA_FAST[deg][i], part);
        }
    }
    reduction_montgomery_int128(out, polynomial, L, L_INV);
}



// convert an element of the extension field given to its representation in PMNS using a pseudo-fast method
void convert_element_to_pmns_pseudo_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    int POL_LIMBS = N_INT_RED_PSEUDO_FAST + 1;
    mp_limb_t vector[DEGREE][POL_LIMBS];

    for (int i=0; i<DEGREE; i++) mpn_zero(vector[i], POL_LIMBS);

    for (int deg=0; deg<EXTENSION_DEGREE; deg++){
        size_t mask_pos = 0;
	    __int128 partial_polynomial[DEGREE] = {0};

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64(partial_polynomial, PMNS_THETA_PSEUDO_FAST[0][i], part);
        }
        addmul_pol128_pol64_c2(POL_LIMBS, vector, partial_polynomial, PMNS_FIELD_ROOTS[deg]);
    }

    for (int it = 0; it < N_INT_RED_PSEUDO_FAST; it++)
        reduction_montgomery_mpn(POL_LIMBS, vector, vector, L, L_INV, it+1);

    for (int i=0; i<DEGREE; i++)
        out[i] = (int64_t)vector[i][0];
}