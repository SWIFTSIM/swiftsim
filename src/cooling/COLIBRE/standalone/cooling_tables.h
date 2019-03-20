#include "cooling_struct.h"

#define SWIFT_STRUCT_ALIGNMENT 32

#define colibre_cooling_N_temperature 86
#define colibre_cooling_N_redshifts 10
#define colibre_cooling_N_density 71
#define colibre_cooling_N_metallicity 11
#define colibre_cooling_N_internalenergy 191
#define colibre_cooling_N_abundances 12
#define colibre_cooling_N_cooltypes 22
#define colibre_cooling_N_heattypes 24
#define colibre_cooling_N_electrontypes 14
#define colibre_cooling_N_elementtypes 12

#define element_H   0
#define element_He  1
#define element_C   2
#define element_N   3
#define element_O   4
#define element_Ne  5
#define element_Mg  6
#define element_Si  7
#define element_S   8
#define element_Ca  9
#define element_Fe 10
#define element_OA 11

#define cooltype_H2         12
#define cooltype_molecules  13
#define cooltype_HD         14
#define cooltype_NetFFH     15
#define cooltype_NetFFM     16
#define cooltype_eeBrems    17
#define cooltype_Compton    18
#define cooltype_Dust       19

#define heattype_H2         12
#define heattype_COdiss     13
#define heattype_CosmicRay  14
#define heattype_UTA        15
#define heattype_line       16
#define heattype_Hlin       17
#define heattype_ChaT       18
#define heattype_HFF        19
#define heattype_Compton    20
#define heattype_Dust       21

void allocate_part_data(struct particle_data *part);
float abundance_ratio_to_solar(struct particle_data *part, struct cooling_function_data *cooling,
    float ratio_solar[colibre_cooling_N_elementtypes]);
void read_cooling_header(struct cooling_function_data *cooling);
void read_cooling_tables(struct cooling_function_data *restrict cooling);

