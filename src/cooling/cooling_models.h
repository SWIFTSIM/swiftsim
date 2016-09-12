#ifndef SWIFT_COOLING_MODELS_H
#define SWIFT_COOLING_MODELS_H

#include "cooling_models_struct.h"

#include "none/cooling.h"
#include "const_du/cooling.h"
#include "const_lambda/cooling.h"

static INLINE void cooling_init_backend(
    const struct swift_params* parameter_file, 
    const struct UnitSystem* us,
    const struct phys_const* phys_const,
    cooling_function_data_handle* cooling_handle_ptr) {

    switch (cooling_model) {
        case (cooling_none) : {
            // FIXME : misses corresponding free
            *cooling_handle_ptr = NULL;
            cooling_none_init_backend(parameter_file, us, phys_const, *cooling_handle_ptr);
            break;
        } 
        case (cooling_const_du) : {
            // FIXME : misses corresponding free
            *cooling_handle_ptr = malloc(sizeof(struct cooling_const_du_function_data));
            cooling_const_du_init_backend(parameter_file, us, phys_const, *cooling_handle_ptr);
            break;
        } 
        case (cooling_const_lambda) : {
            // FIXME : misses corresponding free
            *cooling_handle_ptr = malloc(sizeof(struct cooling_const_lambda_function_data));
            cooling_const_lambda_init_backend(parameter_file, us, phys_const, *cooling_handle_ptr);
            break;
        }
        default:
        {
            error("Selected a non-existing cooling model");
        }
    }
    
}

static INLINE void cooling_print_backend(
    cooling_function_data_handle cooling_handle) {
    switch (cooling_model) {
        case (cooling_none) : {
            cooling_none_print_backend(cooling_handle);
            break;
        }
        case (cooling_const_du) : {
            cooling_const_du_print_backend(cooling_handle);
            break;
        }
        case (cooling_const_lambda) : {
            cooling_const_lambda_print_backend(cooling_handle);
            break;
        }
        default:
        {
            error("Selected a non-existing cooling model");
        }
    }
}

__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct UnitSystem* restrict us,
    cooling_function_data_handle cooling_handle,
    struct part* restrict p, struct xpart* restrict xp, float dt) {
    switch (cooling_model) {
        case (cooling_const_du) : {
            cooling_const_du_cool_part(phys_const, us, cooling_handle, p, xp, dt);
            break;
        }
        case (cooling_const_lambda) : {
            cooling_const_lambda_cool_part(phys_const, us, cooling_handle, p, xp, dt);
            break;
        }
        default:
        {
            error("Selected a non-existing cooling model");
        }
    }
}

__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {
    switch (cooling_model) {
        case (cooling_none) : {
            return cooling_none_get_radiated_energy(xp);
            break;
        }
        case (cooling_const_du) : {
            return cooling_const_du_get_radiated_energy(xp);
            break;
        }
        case (cooling_const_lambda) : {
            return cooling_const_lambda_get_radiated_energy(xp);
            break;
        }
        default:
        {
            error("Selected a non-existing cooling model");
        }
    }
    return 0.0f;
}

__attribute__((always_inline)) INLINE static void cooling_init_part(
    const struct part* restrict p, struct xpart* restrict xp) {
    switch (cooling_model) {
        case (cooling_none) : {
            return cooling_none_init_part(p, xp);
            break;
        }
        case (cooling_const_du) : {
            return cooling_const_du_init_part(p, xp);
            break;
        }
        case (cooling_const_lambda) : {
            return cooling_const_lambda_init_part(p, xp);
            break;
        }
        default:
        {
            error("Selected a non-existing cooling model");
        }
    }    
}

__attribute__((always_inline)) INLINE static float cooling_timestep(
    cooling_function_data_handle cooling,
    const struct phys_const* restrict phys_const,
    const struct UnitSystem* restrict us, const struct part* restrict p) {
    switch (cooling_model) {
        case (cooling_none) : {
            return cooling_none_timestep(cooling, phys_const, us, p);
            break;
        }
        case (cooling_const_du) : {
            return cooling_const_du_timestep(cooling, phys_const, us, p);
            break;
        }
        case (cooling_const_lambda) : {
            return cooling_const_lambda_timestep(cooling, phys_const, us, p);
            break;
        }
        default:
        {
            error("Selected a non-existing cooling model");
        }        
    }    
    return FLT_MAX;
}
#endif /* SWIFT_COOLING_MODELS_H */

