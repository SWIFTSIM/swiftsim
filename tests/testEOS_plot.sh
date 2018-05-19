#!/bin/bash

echo ""

rm -f testEOS*.png

echo "Plotting testEOS output for each planetary material"

A1_mat_id=(
    100
    101
    102
    200
    201
    202
    401
    402
)

for mat_id in "${A1_mat_id[@]}"
do
    python ./testEOS.py "$mat_id"
done

exit $?
