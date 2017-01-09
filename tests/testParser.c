/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 James Willis (james.s.willis@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "parser.h"

int main(int argc, char *argv[]) {
  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Print the contents of the structure to stdout. */
  parser_print_params(&param_file);

  /* Print the contents of the structure to a file in YAML format. */
  parser_write_params_to_file(&param_file, "parser_output.yml");

  /* Retrieve parameters and store them in variables defined above.
   * Have to specify the name of the parameter as it appears in the
   * input file: testParserInput.yaml.*/
  const int no_of_threads =
      parser_get_param_int(&param_file, "Scheduler:no_of_threads");
  const int no_of_time_steps =
      parser_get_param_int(&param_file, "Simulation:no_of_time_steps");
  const float max_h = parser_get_param_float(&param_file, "Simulation:max_h");
  const double start_time =
      parser_get_param_double(&param_file, "Simulation:start_time");
  const int kernel = parser_get_param_int(&param_file, "kernel");

  const int optional = parser_get_opt_param_int(&param_file, "optional", 1);

  char ic_file[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(&param_file, "IO:ic_file", ic_file);

  /* Print the variables to check their values are correct. */
  printf(
      "no_of_threads: %d, no_of_time_steps: %d, max_h: %f, start_time: %lf, "
      "ic_file: %s, kernel: %d optional: %d\n",
      no_of_threads, no_of_time_steps, max_h, start_time, ic_file, kernel,
      optional);

  assert(no_of_threads == 16);
  assert(no_of_time_steps == 10);
  assert(fabs(max_h - 1.1255) < 0.00001);
  assert(fabs(start_time - 1.23456789) < 0.00001);
  assert(strcmp(ic_file, "ic_file.ini") == 0); /*strcmp returns 0 if correct.*/
  assert(kernel == 4);
  assert(optional == 1);

  /* Test parsing boolean values. */
  assert(parser_get_param_bool(&param_file, "Booleans:true"));
  assert(parser_get_param_bool(&param_file, "Booleans:True"));
  assert(parser_get_param_bool(&param_file, "Booleans:yes"));
  assert(parser_get_param_bool(&param_file, "Booleans:1"));
  assert(!parser_get_param_bool(&param_file, "Booleans:false"));
  assert(!parser_get_param_bool(&param_file, "Booleans:False"));
  assert(!parser_get_param_bool(&param_file, "Booleans:no"));
  assert(!parser_get_param_bool(&param_file, "Booleans:0"));

  return 0;
}
