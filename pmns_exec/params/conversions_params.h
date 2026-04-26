# include "pmns_params.h"

# define N_INT_RED_CLASSICAL 2
# define N_INT_RED_PSEUDO_FAST 2
# define N_INT_RED_FAST 1
# define THETA_POW 43
# define N_POL 3

static const mp_limb_t TRANSITION_MATRIX[EXTENSION_DEGREE][EXTENSION_DEGREE][N_LIMBS] = {
    {
        {1ULL, 0ULL}
    }
};

//~ pre-compute PMNS representation of theta^i * phi^r
static const int64_t PMNS_THETA_PSEUDO_FAST[1][N_POL][DEGREE] = {
{
    {230370849217ULL, 2904625845712ULL, -753150678476LL},
    {-120453392714LL, -1301687565111LL, 4225326649ULL},
    {54651646326ULL, -918032577727LL, -1575541337725LL}
}
};

//~ pre-compute PMNS representation of z^i where z is the root of the polynomial use to construct extension field
static const int64_t PMNS_FIELD_ROOTS[EXTENSION_DEGREE][DEGREE] = {
    {230370849217ULL, 2904625845712ULL, -753150678476LL}
};

//~ pre-compute PMNS representation of theta^i * phi^r * z**j
static const int64_t PMNS_THETA_FAST[EXTENSION_DEGREE][N_POL][DEGREE] = {
{
    {2895211771650ULL, 423443904370ULL, -3515315866999LL},
    {1475261082908ULL, 882162105513ULL, -565686007653LL},
    {-685892311643LL, 2111694920052ULL, -1207751981462LL}
}
};