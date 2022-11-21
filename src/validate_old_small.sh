#!/usr/bin/bash

# Compute the MSE values for an old model

MODEL_LEVELS=(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

for l in ${MODEL_LEVELS[@]}; do
    ./validate --data ../../csv/ETOPO1_small.csv --npoint 64800 --lmax $l --model default_${l}_180_360.txt
done