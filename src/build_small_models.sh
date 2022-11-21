#!/usr/bin/bash

# Build the small models with order = {0, 5, 10, 20, 50, 100, 150, 200, 300}
# the executable shall be run from the source folder

# ORDERS=(0, 5, 10, 20, 50, 100, 150, 200, 300)
ORDERS=(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
SUP_LMAX=500

for o in ${ORDERS[@]}; do
    ./new_model_omp2 --small --lmodel $o --lmax $SUP_LMAX --txt --predict --diff
done