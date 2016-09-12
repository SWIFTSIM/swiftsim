#ifndef SWIFT_COOLING_MODELS_STRUCT_H
#define SWIFT_COOLING_MODELS_STRUCT_H

enum cooling_models {
    cooling_none,
    cooling_const_du,
    cooling_const_lambda,
    cooling_grackle
};

static const int cooling_model = cooling_const_lambda;

// Use anonymous void * handles to permit selection at run-time
typedef void * cooling_function_data_handle;
typedef void * cooling_xpart_data_handle;

#endif /* SWIFT_COOLING_MODELS_STRUCT_H */

