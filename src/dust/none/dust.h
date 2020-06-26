#ifndef SWIFT_DUST_NONE_H
#define SWIFT_DUST_NONE_H


/**
 * @brief redistribute any dust mass back to element abundances
 * on star particle formation according to dust composition, to 
 * represent astration 
 *
 * Nothing here.
 *
 * @param p The gas particles.  
 * @param dp Global dust parameters for initialisation.
 */
void redistribute_dust_masses(const struct part* p, 
			      struct dustevo_props *dp) {};

#endif /* SWIFT_DUST_NONE_PROPERTIES_H */
