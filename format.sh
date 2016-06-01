#!/bin/bash

clang-format -style=file -i src/*.[ch] src/*/*.[ch] src/*/*/*.[ch] examples/main.c tests/*.[ch]
