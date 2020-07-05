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

#endif /* SWIFT_DUST_T20_PROPERTIES_H */
