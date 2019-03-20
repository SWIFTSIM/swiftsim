#include <string.h>
#include "cooling_function.h"

int main(int argc, char **argv) {

  struct cooling_function_data cooling;
  bzero(&cooling, sizeof(struct cooling_function_data));

  struct particle_data part;
  bzero(&part, sizeof(struct particle_data));

  /* Initialize cooling tables */
  cooling_init_backend(&cooling, &part);

  /* Print cooling function for the parameters set there */
  cooling_function(&cooling, &part);

}
