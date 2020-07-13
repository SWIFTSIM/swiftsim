/* Local includes. */

#include "chemistry.h"
#include "dust_properties.h"

/**
 * @brief redistribute any dust mass back to element abundances
 * on star particle formation according to dust composition, to 
 * represent astration 
 *
 * @param p The gas particles.  
 * @param dp Global dust parameters for initialisation.
 */
void redistribute_dust_masses(const struct part* p, 
			      struct dustevo_props *dp) {
  /** 
   * iterate through grain species and element types and
   * redistribute dust abundances to element abundances,
   * according to composition array
   */
  ;
}

/**
 * @brief Prints the dust evolution model to stdout.
 *
 * @param dust #dustevo_props struct.
 */
void dustevo_print_backend(const struct dustevo_props *dp) {
    message("Running with a Trayford et al. (2020) dust evolution model.");
}
