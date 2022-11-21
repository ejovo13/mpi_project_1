#!/usr/bin/bash

# Build the old model working on the small data set to test error at different model
# levels

MODEL_LEVELS=(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

for l in ${MODEL_LEVELS[@]}; do
    ./model --data ../../csv/ETOPO1_small.csv --npoint 64800 --lmax $l --model default_${l}_180_360.txt
done

