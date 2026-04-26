#ifndef CONVERSION_API
#define CONVERSION_API

#include "../params/pmns_params.h"

void convert_element_to_pmns_classical(int64_t out[DEGREE], mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

void convert_element_to_pmns_pseudo_fast_nbc(int64_t out[DEGREE], mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

void convert_element_to_pmns_fast(int64_t out[DEGREE], mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

#endif