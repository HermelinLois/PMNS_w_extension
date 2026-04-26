# include "reductions_interface.h"
# include <stdio.h>
# include <gmp.h>
#include <inttypes.h>

static mpz_t PHI;
static int PHI_ready = 0;
static const mp_limb_t phi_limbs[2] = {0ULL, 1ULL};


static inline void ensure_PHI() {
    if (PHI_ready)
        return;

    mpz_init(PHI);
    mpz_import(PHI, 2, -1, sizeof(mp_limb_t), 0, 0, phi_limbs);
    PHI_ready = 1;
}


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















# define SHOW 0











void reduction_montgomery_mpn_test(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t abs_polynomial[DEGREE][N_LIMBS], int coeffs_sign[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE], int it){
    int ACC_SIZE = N_LIMBS + 1;
    mp_limb_t acc[ACC_SIZE];
    int64_t Q[DEGREE];
    mp_limb_t T[DEGREE][N_LIMBS];
    int T_sign[DEGREE];
    mp_limb_t tmp[N_LIMBS];

    for (int j=0; j<DEGREE; j++) {
        mpn_zero(acc, ACC_SIZE);

        for (int i=0; i<DEGREE; i++) {
            int64_t mat_coeff = sublattice_inv[i][j];
            if (mat_coeff == 0 || mpn_zero_p(abs_polynomial[i], N_LIMBS)) continue;

            int mat_sign = (mat_coeff < 0);
            uint64_t abs_mat_coef = mat_sign ? (uint64_t)(-mat_coeff) : (uint64_t)mat_coeff;

            int prod_sign = mat_sign ^ coeffs_sign[i];

            if (prod_sign) {
                mp_limb_t borrow = mpn_submul_1(acc, abs_polynomial[i], N_LIMBS, abs_mat_coef);
                mpn_sub_1(acc + N_LIMBS, acc + N_LIMBS, 1, borrow);
            } else {
                mp_limb_t carry = mpn_addmul_1(acc, abs_polynomial[i], N_LIMBS, abs_mat_coef);
                mpn_add_1(acc + N_LIMBS, acc + N_LIMBS, 1, carry);
            }
        }
        Q[j] = (int64_t)acc[0];
    }

    if (it == SHOW) {
        printf("=== Q mpn sign-mag===\n");
        for (int j=0; j<DEGREE; j++)
            printf("  Q[%d] = %ld\n", j, Q[j]);
    }

    for (int j=0; j<DEGREE; j++) {
        mpn_zero(T[j], N_LIMBS);
        int sign = 0;

        for (int i=0; i<DEGREE; i++){
            int64_t mat_coeff = sublattice[i][j];
            int64_t Q_coeff   = Q[i];

            if (mat_coeff == 0 || Q_coeff == 0) continue;

            int Q_sign = (Q_coeff < 0);
            int mat_sign = (mat_coeff < 0);
            
            uint64_t mat_abs = mat_sign ? -(uint64_t)mat_coeff : (uint64_t)mat_coeff;
            uint64_t Q_abs = Q_sign ? -(uint64_t)Q_coeff   : (uint64_t)Q_coeff;

            int sign_prod = mat_sign ^ Q_sign;
            __int128 prod = (__int128)mat_abs * Q_abs;
            
            mpn_zero(tmp, N_LIMBS);

            tmp[0] = (mp_limb_t)prod;
            tmp[1] = (mp_limb_t)(prod >> GMP_NUMB_BITS);

            if (sign == sign_prod) {
                mpn_add_n(T[j], T[j], tmp, N_LIMBS);
            } else {
                if (mpn_cmp(T[j], tmp, N_LIMBS) >= 0) {
                    mpn_sub_n(T[j], T[j], tmp, N_LIMBS);
                } else {
                    mpn_sub_n(T[j], tmp, T[j], N_LIMBS);
                    sign = sign_prod;
                }
            }
        }
        T_sign[j] = sign;
    }

    if (it == SHOW) {
        printf("=== T mpn sign-mag===\n");
        for (int j=0; j<DEGREE; j++)
            printf("  T[%d] = %ld  (sign=%d)\n", j, (int64_t)T[j][0], T_sign[j]);
    }

    for (int deg = 0; deg < DEGREE; deg++) {
        mpn_zero(acc, ACC_SIZE);

        int sP = coeffs_sign[deg];
        int sT = T_sign[deg];

        mp_limb_t *P_val = abs_polynomial[deg];
        mp_limb_t *T_val = T[deg];

        int sign;

        if (sP == sT) {
            mpn_add_n(acc, P_val, T_val, N_LIMBS);
            sign = sP;
        } 
        else {
            if (mpn_cmp(P_val, T_val, N_LIMBS) >= 0) {
                mpn_sub_n(acc, P_val, T_val, N_LIMBS);
                sign = sP;
            } else {
                mpn_sub_n(acc, T_val, P_val, N_LIMBS);
                sign = sT;
            }
        }

        
        if (it == SHOW)
            printf("=== ADD SHIFT mpn sign-mag===\n");
            printf("  P[%d]=%ld (sign=%d)  T[%d]=%ld (sign=%d)  out=%ld (sign=%d)\n",
                   deg, (int64_t)P_val[0], sP,
                   deg, (int64_t)T_val[0], sT,
                   (int64_t)acc[1], sign);

        coeffs_sign[deg] = sign;
        mpn_copyi(out[deg], acc + 1, N_LIMBS);
    }
}




























static inline void mpn_asr_1(mp_limb_t *dst, const mp_limb_t *src, size_t n, int sign){
    mpn_copyi(dst, src + 1, n - 1);

    if (sign)
        dst[n - 1] = ~(mp_limb_t)0;
    else
        dst[n - 1] = 0;
}


// void reduction_montgomery_mpn(mp_limb_t out[DEGREE][N_LIMBS], mp_limb_t P[DEGREE][N_LIMBS], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]){
//     mp_limb_t T[DEGREE][N_LIMBS];
//     mp_limb_t tmp[N_LIMBS];
//     mp_limb_t acc[N_LIMBS + 1];
//     int64_t Q[DEGREE];

//     for (int j=0; j<DEGREE; j++){
//         int64_t acc = 0;
        
//         for (int i=0; i<DEGREE; i++) {
//             int64_t mat_coeff = sublattice_inv[i][j];
//             if (mat_coeff == 0) continue;
//             acc += mat_coeff * P[i][0];
//         }
//         Q[j] = acc;
//     }

//     for (int j = 0; j < DEGREE; j++) {
//         mpn_zero(T[j], N_LIMBS);

//         for (int i = 0; i < DEGREE; i++) {

//             int64_t mat_coeff = sublattice[i][j];
//             int64_t Q_coeff = Q[i];

//             if (mat_coeff == 0 || Q_coeff == 0) continue;

//             int sign = (mat_coeff < 0) ^ (Q_coeff < 0);

//             uint64_t umat_coeff = (mat_coeff < 0) ? -(uint64_t)mat_coeff : (uint64_t)mat_coeff;
//             uint64_t uQ_coeff = (Q_coeff < 0) ? -(uint64_t)Q_coeff : (uint64_t)Q_coeff;

//             __int128 prod = (__int128) umat_coeff * uQ_coeff;

//             mpn_zero(tmp, N_LIMBS);

//             tmp[0] = (mp_limb_t)prod;
//             if (N_LIMBS > 1)
//                 tmp[1] = (mp_limb_t)(prod >> GMP_NUMB_BITS);

//             if (sign){
//                 mpn_com(tmp, tmp, N_LIMBS);
//                 mpn_add_1(tmp, tmp, N_LIMBS, 1);
//             }

//             mpn_add_n(T[j], T[j], tmp, N_LIMBS);
//         }
//     }

//     for (int deg = 0; deg < DEGREE; deg++) {
//         mpn_zero(acc, N_LIMBS + 1);

//         acc[N_LIMBS] = mpn_add_n(acc, P_abs, T_abs, N_LIMBS);
//         int sign = (acc[N_LIMBS - 1] >> (GMP_NUMB_BITS - 1)) & 1;

//         mpn_copyi(out, acc + 1, N_LIMBS);

//         if (sign) {
//             mpn_com(out, out, N_LIMBS);
//             mpn_add_1(out, out, N_LIMBS, 1);
//         }
//     }
// }
















void reduction_montgomery_mpn(
    mp_limb_t out[DEGREE][N_LIMBS],
    mp_limb_t P[DEGREE][N_LIMBS],
    const int64_t sublattice[DEGREE][DEGREE],
    const int64_t sublattice_inv[DEGREE][DEGREE],
    int it)
{
    mp_limb_t T[DEGREE][N_LIMBS + 1];
    mp_limb_t tmp[N_LIMBS + 1];       // ← +1
    mp_limb_t acc[N_LIMBS + 1];
    int64_t Q[DEGREE];

    for (int j = 0; j < DEGREE; j++) {
        int64_t q = 0;
        for (int i = 0; i < DEGREE; i++)
            q += sublattice_inv[i][j] * (int64_t)P[i][0];
        Q[j] = q;
    }
    mp_limb_t T_pos[DEGREE][N_LIMBS + 1];  // somme des termes positifs
    mp_limb_t T_neg[DEGREE][N_LIMBS + 1];  // somme des termes négatifs (magnitudes)

    for (int j = 0; j < DEGREE; j++) {
        mpn_zero(T_pos[j], N_LIMBS + 1);
        mpn_zero(T_neg[j], N_LIMBS + 1);

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t Q_coeff   = Q[i];
            if (mat_coeff == 0 || Q_coeff == 0) continue;

            int mat_sign  = (mat_coeff < 0);
            int q_sign    = (Q_coeff   < 0);
            int prod_sign = mat_sign ^ q_sign;

            uint64_t umat = mat_sign ? -(uint64_t)mat_coeff : (uint64_t)mat_coeff;
            uint64_t uq   = q_sign   ? -(uint64_t)Q_coeff   : (uint64_t)Q_coeff;

            __int128 prod = (__int128)umat * uq;

            mp_limb_t tmp[N_LIMBS + 1];
            mpn_zero(tmp, N_LIMBS + 1);
            tmp[0] = (mp_limb_t)prod;
            tmp[1] = (mp_limb_t)(prod >> GMP_NUMB_BITS);

            if (prod_sign)
                mpn_add_n(T_neg[j], T_neg[j], tmp, N_LIMBS + 1);
            else
                mpn_add_n(T_pos[j], T_pos[j], tmp, N_LIMBS + 1);
        }
    }

    if (it == SHOW) {
        printf("=== T mpn com 2===\n");
        for (int j = 0; j < DEGREE; j++) {
            int t_sign = (T[j][N_LIMBS] >> (GMP_NUMB_BITS - 1)) & 1;  // ← N_LIMBS
            printf("  T[%d] = %ld  (sign=%d)\n", j, (int64_t)T[j][0], t_sign);
        }
    }

    for (int deg = 0; deg < DEGREE; deg++) {
        mp_limb_t P_ext[N_LIMBS + 2];
        mp_limb_t acc[N_LIMBS + 2];
        mp_limb_t tmp2[N_LIMBS + 2];

        int p_sign = (P[deg][N_LIMBS - 1] >> (GMP_NUMB_BITS - 1)) & 1;

        mpn_copyi(P_ext, P[deg], N_LIMBS);
        P_ext[N_LIMBS]     = p_sign ? ~(mp_limb_t)0 : 0;
        P_ext[N_LIMBS + 1] = p_sign ? ~(mp_limb_t)0 : 0;

        // étendre T_pos et T_neg sur N_LIMBS+2
        mp_limb_t Tp[N_LIMBS + 2], Tn[N_LIMBS + 2];
        mpn_copyi(Tp, T_pos[deg], N_LIMBS + 1); Tp[N_LIMBS + 1] = 0;
        mpn_copyi(Tn, T_neg[deg], N_LIMBS + 1); Tn[N_LIMBS + 1] = 0;

        // acc = P + T_pos - T_neg
        mpn_add_n(acc, P_ext, Tp, N_LIMBS + 2);
        mpn_sub_n(acc, acc, Tn, N_LIMBS + 2);

        int res_sign = (acc[N_LIMBS + 1] >> (GMP_NUMB_BITS - 1)) & 1;
        mp_limb_t tmp_out[N_LIMBS + 1];
        mpn_asr_1(tmp_out, acc, N_LIMBS + 2, res_sign);
        mpn_copyi(out[deg], tmp_out, N_LIMBS);

        if (it == SHOW) {
            printf("=== ADD SHIFT mpn com 2===\n");
            printf("  P[%d]=%ld (sign=%d)  T[%d]=%ld (sign=%d)  out=%ld (sign=%d)\n",
                deg, (int64_t)P[deg][0], p_sign,
                deg, (int64_t)T[deg][0], (T[deg][N_LIMBS + 1] >> 63)&1,
                (int64_t)out[deg][0], res_sign);
        }
    }
}