#ifndef PMNS_PARAMS_H
#define PMNS_PARAMS_H

#include <gmp.h>
#include <stdint.h>

# define DEGREE 3
# define RHO 12326944617285
# define N_LIMBS 2
# define EXTENSION_DEGREE 1
# define PHI_POW 64
# define N_TEST 100

static mp_limb_t GAMMA[EXTENSION_DEGREE][N_LIMBS] = {
    {4784070557864973851ULL, 9146308817568836342ULL}
};

static mp_limb_t P[N_LIMBS] = {1968545615899564403ULL, 10618920251238660050ULL};

#endif