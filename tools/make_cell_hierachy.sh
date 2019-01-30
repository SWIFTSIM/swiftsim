#!/bin/bash

output=cell_hierachy.html

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cp $DIR/data/cell_hierachy.html $output

echo $output has been generated
