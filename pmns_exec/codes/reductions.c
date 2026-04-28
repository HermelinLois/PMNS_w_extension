# include "reductions_interface.h"
# include <stdio.h>
# include <gmp.h>

//static mpz_t PHI;
//static int PHI_ready = 0;
//static const mp_limb_t phi_limbs[2] = {0ULL, 1ULL};


//static inline void ensure_PHI() {
//    if (PHI_ready)
//        return;
//
//    mpz_init(PHI);
//    mpz_import(PHI, 2, -1, sizeof(mp_limb_t), 0, 0, phi_limbs);
//    PHI_ready = 1;
//}


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

    ext_red_w_matrix(uint64_t, Q, polynomial, sublattice_inv);

    for (int deg = 0; deg < DEGREE; deg++)
        Q[deg] = (uint64_t)Q[deg];

    ext_red_w_matrix(__int128, T, Q, sublattice);

    for (int deg = 0; deg < DEGREE; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}



void reduction_montgomery_mpn(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t P[DEGREE][N_LIMBS], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE], int it){
    int SIZE = N_LIMBS + 1;

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
            
            mp_limb_t prod_limbs[2] = {(mp_limb_t)prod, (mp_limb_t)(prod >> GMP_NUMB_BITS)};
            int sign_prod = mat_sign ^ Q_sign;

            if (sign_prod) {
                mp_limb_t borrow = mpn_sub_n(T[j], T[j], prod_limbs, 2);
                mpn_sub_1(T[j] + 2, T[j] + 2, SIZE - 2, borrow);
            } else {
                mp_limb_t carry = mpn_add_n(T[j], T[j], prod_limbs, 2);
                mpn_add_1(T[j] + 2, T[j] + 2, SIZE - 2, carry);
            }
        }
    }


    for (int deg = 0; deg < DEGREE; deg++) {
        mpn_zero(tmp, SIZE);

        mp_limb_t *P_val = P[deg];
        int sP = (it > 0) ? ((P_val[N_LIMBS - 1] >> (GMP_NUMB_BITS - 1)) & 1) : 0; 

        mp_limb_t coeff_extended[SIZE];
        mpn_copyi(coeff_extended, P_val, N_LIMBS);
        coeff_extended[N_LIMBS] = sP ? ~0UL : 0UL;
        
        mpn_add_n(tmp, coeff_extended, T[deg], SIZE);
        mpn_copyi(out[deg], tmp + 1, N_LIMBS);

        int is_negative = (tmp[SIZE - 1] >> (GMP_NUMB_BITS - 1)) & 1;
        if (is_negative && tmp[0] != 0) {
            mpn_sub_1(out[deg], out[deg], N_LIMBS, 1);
        }
    }
}


// void reduction_montgomery_mpz(mpz_t out[DEGREE], mpz_t polynomial[DEGREE],  const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
//     mpz_t Q[DEGREE], T[DEGREE];

//     for (int deg=0; deg<DEGREE; deg++){
//         mpz_init_set_ui(Q[deg], 0);
//         mpz_init_set_ui(T[deg], 0);
//     }

//     ensure_PHI();

//     // compute (polynomial * sublattice_inv) % phi
//     for (int j=0; j<DEGREE; j++){
//         for (int i=0; i<DEGREE; i++){
//             int64_t coeff = sublattice_inv[i][j];
//             if (coeff == 0) continue;
            
//             int sign = (coeff < 0);
//             uint64_t abs_coef = (sign ? (uint64_t) -coeff : (uint64_t) coeff);

//             (sign ? mpz_submul_ui : mpz_addmul_ui)(Q[j], polynomial[i], abs_coef);
//         }
//         mpz_fdiv_r_2exp(Q[j], Q[j], PHI_POW);
        
//         if (mpz_tstbit(Q[j], PHI_POW - 1))
//             mpz_sub(Q[j], Q[j], PHI);
//     }
    
//     // compute T = Q * sublattice            
//     for (int i=0; i<DEGREE; i++){                           
//         for (int j=0; j<DEGREE; j++){
//             int64_t coeff = sublattice[i][j];
//             if (coeff == 0) continue;
            
//             int sign = (coeff < 0);
//             uint64_t abs_coef = (sign ? (uint64_t) -coeff : (uint64_t) coeff);

//             (sign ? mpz_submul_ui : mpz_addmul_ui)(T[j], Q[i], abs_coef);
//         }
//         mpz_clear(Q[i]);
//     }

//     for (int deg=0; deg<DEGREE; deg++){
//         mpz_add(out[deg], polynomial[deg], T[deg]);
//         mpz_fdiv_q_2exp(out[deg], out[deg], PHI_POW);
//         mpz_clear(T[deg]);
//     }        
// }







void reduction_montgomery_mpn_test(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t abs_polynomial[DEGREE][N_LIMBS], int coeffs_sign[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]){
    int ACC_SIZE = N_LIMBS + 1;
    mp_limb_t acc[ACC_SIZE];
    int64_t Q[DEGREE];
    mp_limb_t T[DEGREE][ACC_SIZE];
    int T_sign[DEGREE];
    mp_limb_t tmp[ACC_SIZE];

    for (int j = 0; j < DEGREE; j++) {
        uint64_t acc64 = 0;
        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice_inv[i][j];
            if (mat_coeff == 0) continue;
            unsigned __int128 p_abs = (unsigned __int128)abs_polynomial[i][1] << GMP_NUMB_BITS | abs_polynomial[i][0];

            unsigned __int128 p_c2 = coeffs_sign[i] ? -p_abs : p_abs;
            uint64_t p_low = (uint64_t)p_c2;
            acc64 += p_low * (uint64_t)mat_coeff;
        }
        Q[j] = (int64_t)acc64;
    }

    for (int j = 0; j < DEGREE; j++) {
        mpn_zero(T[j], ACC_SIZE);
        int sign = 0;

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t Q_coeff   = Q[i];
            if (mat_coeff == 0 || Q_coeff == 0) continue;

            int Q_sign   = (Q_coeff < 0);
            int mat_sign = (mat_coeff < 0);
            uint64_t mat_abs = mat_sign ? -(uint64_t)mat_coeff : (uint64_t)mat_coeff;
            uint64_t Q_abs   = Q_sign   ? -(uint64_t)Q_coeff   : (uint64_t)Q_coeff;

            int sign_prod = mat_sign ^ Q_sign;
            __int128 prod = (__int128)mat_abs * Q_abs;

            mpn_zero(tmp, ACC_SIZE);
            tmp[0] = (mp_limb_t)prod;
            tmp[1] = (mp_limb_t)(prod >> GMP_NUMB_BITS);

            if (sign == sign_prod) {
                mpn_add_n(T[j], T[j], tmp, ACC_SIZE);
            } else {
                if (mpn_cmp(T[j], tmp, ACC_SIZE) >= 0) {
                    mpn_sub_n(T[j], T[j], tmp, ACC_SIZE);
                } else {
                    mpn_sub_n(T[j], tmp, T[j], ACC_SIZE);
                    sign = sign_prod;
                }
            }
        }
        T_sign[j] = sign;
    }

    for (int deg = 0; deg < DEGREE; deg++) {
        mpn_zero(acc, ACC_SIZE);
        int sP   = coeffs_sign[deg];
        int sT   = T_sign[deg];
        mp_limb_t *P_val = abs_polynomial[deg];
        mp_limb_t *T_val = T[deg];
        int sign;

        if (sP == sT) {
            mp_limb_t carry = mpn_add_n(acc, P_val, T_val, N_LIMBS);
            acc[N_LIMBS] = carry;
            sign = sP;
        } else {
            if (mpn_cmp(P_val, T_val, N_LIMBS) >= 0) {
                mpn_sub_n(acc, P_val, T_val, N_LIMBS);
                sign = sP;
            } else {
                mpn_sub_n(acc, T_val, P_val, N_LIMBS);
                sign = sT;
            }
        }

        coeffs_sign[deg] = sign;
        mpn_copyi(out[deg], acc + 1, N_LIMBS);
    }
}







































