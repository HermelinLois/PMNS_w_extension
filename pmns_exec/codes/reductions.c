# include "reductions_interface.h"
# include <stdio.h>
# include <gmp.h>

#define DEFINE_SHIFT(NAME, OUT_TYPE, IN_TYPE)                                                   \
static inline void NAME(OUT_TYPE out[DEGREE], const IN_TYPE polynomial[DEGREE], int n_shift) {  \
    for (int i = 0; i < DEGREE; i++) {                                                          \
        out[i] = (OUT_TYPE)(polynomial[i] >> n_shift);                                          \
    }                                                                                           \
}

#define DEFINE_VECT_MAT_PROD(NAME, IN_TYPE)                                                                         \
static inline void NAME(__int128 out[DEGREE], const IN_TYPE vector[DEGREE], const int64_t matrix[DEGREE][DEGREE]) { \
    for (int i = 0; i < DEGREE; i++) {                                                                              \
        __int128 coeff = 0;                                                                                         \
        for (int j = 0; j < DEGREE; j++)                                                                            \
            coeff += (__int128)vector[j] * (__int128)matrix[j][i];                                                  \
        out[i] = coeff;                                                                                             \
    }                                                                                                               \
}

DEFINE_SHIFT(shift64, int64_t, __int128)
DEFINE_SHIFT(shift128, __int128, __int128)

DEFINE_VECT_MAT_PROD(vm_prod64, int64_t)
DEFINE_VECT_MAT_PROD(vm_prod128, __int128)

void reduction_babai_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
    __int128 S[DEGREE];
    __int128 SL[DEGREE];
    int64_t PH2[DEGREE];

    shift64(PH2, polynomial, H2);
    vm_prod64(S, PH2, sublattice_inv);
    shift128(S, S, (H1 - H2));
    vm_prod128(SL, S, sublattice);

    for (int i = 0; i < DEGREE; i++)
        out[i] = polynomial[i] - (int64_t)SL[i];
}


void reduction_montgomery_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
    int64_t Q[DEGREE] = {0};
    __int128 T[DEGREE] = {0};

    prod_pol_mat(uint64_t, Q, polynomial, sublattice_inv);

    for (int deg = 0; deg < DEGREE; deg++)
        Q[deg] = (uint64_t)Q[deg];

    prod_pol_mat(__int128, T, Q, sublattice);

    for (int deg = 0; deg < DEGREE; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}



void reduction_montgomery_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t P[DEGREE][n_limbs], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE], int it){
    int SIZE = n_limbs + 1;

    int64_t Q[DEGREE];
    mp_limb_t T[DEGREE][SIZE];
    mp_limb_t tmp[SIZE];

    for (int j = 0; j < DEGREE; j++) {
        uint64_t acc = 0;
        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice_inv[i][j];
            if (mat_coeff == 0) continue;

            unsigned __int128 P_part = P[i][0] | ((unsigned __int128)P[i][1] << GMP_NUMB_BITS);
            uint64_t p_low = (uint64_t)P_part;

            acc += p_low * (uint64_t)mat_coeff;
        }
        Q[j] = (int64_t)acc;
    }


    for (int j = 0; j < DEGREE; j++) {
        mpn_zero(T[j], SIZE);

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t Q_coeff   = Q[i];
            if (mat_coeff == 0 || Q_coeff == 0) continue;
            
            int mat_sign = (mat_coeff < 0);
            int Q_sign = (Q_coeff < 0);
            unsigned __int128 mat_abs = mat_sign ? -mat_coeff : mat_coeff;
            unsigned __int128 Q_abs = Q_sign ? -Q_coeff : Q_coeff;
            unsigned __int128 prod = mat_abs * Q_abs;
            
            mp_limb_t prod_limbs[SIZE];
            mpn_zero(prod_limbs, SIZE);
            
            prod_limbs[0] = (mp_limb_t)prod;
            prod_limbs[1] = (mp_limb_t)(prod >> GMP_NUMB_BITS);
            int sign_prod = mat_sign ^ Q_sign;

            if (sign_prod) {
                mpn_com(prod_limbs, prod_limbs, SIZE);
                mpn_add_1(prod_limbs, prod_limbs, SIZE, 1);
            }
            mpn_add_n(T[j], T[j], prod_limbs, SIZE);
        }
    }


    for (int deg = 0; deg < DEGREE; deg++) {
        mpn_zero(tmp, SIZE);

        mp_limb_t *P_val = P[deg];
        int sP = (it > 0) ? ((P_val[n_limbs - 1] >> (GMP_NUMB_BITS - 1)) & 1) : 0;

        mp_limb_t coeff_extended[SIZE];
        mpn_copyi(coeff_extended, P_val, n_limbs);
        coeff_extended[n_limbs] = sP ? ~0UL : 0UL;
        
        mpn_add_n(tmp, coeff_extended, T[deg], SIZE);
        mpn_copyi(out[deg], tmp + 1, n_limbs);
    }
}