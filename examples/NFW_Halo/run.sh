#!/bin/bash

# self gravity G, external potential g, hydro s, threads t and high verbosity v
../swift -g -t 6 test.yml 2>&1 | tee output.log


