
#include "cooling_tables.h"

double colibre_cooling_rate(
    double log_u_cgs, double redshift, double n_H_cgs, float ZZsol,
    const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
    float d_n_H, int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling);
 
