#ifndef REDUCTION_API
#define REDUCTION_API

#include <stdint.h>
#include "../params/pmns_params.h"
#include "../params/reductions_params.h"

#define ext_red_w_matrix(RETURN_TYPE, RESULT, POLYNOMIAL, MATRIX){  \
	for (int j = 0; j < DEGREE; j++) {                              \
	    RETURN_TYPE pi = (RETURN_TYPE)(POLYNOMIAL)[j];              \
		for (int i = 0; i < DEGREE; i++)                            \
			(RESULT)[i] += pi * (MATRIX)[j][i];                     \
     }                                                              \
} 

void reduction_montgomery_mpn(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t P[DEGREE][N_LIMBS], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE], int it);

void reduction_montgomery_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);

void reduction_babai_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);

#endif