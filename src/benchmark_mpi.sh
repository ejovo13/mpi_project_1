#!/usr/bin/bash

# usage: ./benchmark_mpi <degree_model> <number_trials>
#
# ex: ./benchmark 100 5 
#
#     - compute timing information and store it in mpi_100_5.csv to determine scalability of this
#       algorithm, run $number_trials times.

# Let's first create a function that will run an mpi model on different data sets
# mpi_s(np, lmodel, lmax)
function mpi_s {
    mpirun -np $1 ./new_model_mpi --small --lmodel $2 --lmax $3 --recompute
}

# mpi_s_cmd(np, lmodel, lmax) to get a _str_ of the
# command used to create a small model
function mpi_s_cmd {
    str="mpirun -np $1 ./new_model_mpi --small --lmodel $2 --lmax $3 --recompute"
    echo $str
}

# mpi_m(np, lmodel, lmax)
function mpi_m {
    mpirun -np $1 ./new_model_mpi --med --lmodel $2 --lmax $3 --recompute
}

# mpi_h(np, lmodel, lmax)
function mpi_h {
    mpirun -np $1 ./new_model_mpi --hi --lmodel $2 --lmax $3 --recompute
}

# extract_time_s(np, lmodel, lmax)
function extract_time {
    time=$(mpi_s $1 $2 $3 | grep Computed | awk '{print $NF}')
    echo ${time::-1}
}

# time_model(degree, n_trial)
# Create a csv with the timing data
# of computing a model with 1 to 12 mpi
# processes
function time_model {

    file_out="mpi_$1_$2.csv"
    echo "Saving timing information to $file_out"
    printf "run," > $file_out
    
    for np in {1..11}; do
        printf "np$np," >> $file_out
    done 
    printf "np12\n" >> $file_out

    echo "Timing model of degree $1 with $2 trials"

    for run in $(seq $2); do

        echo "=================================="
        echo "run: $run"
        echo "=================================="
        printf "$run," >> $file_out

        for np in {1..11}; do
            time=$(extract_time $np $1 $1) 
            echo $time
            printf "$time," >> $file_out
        done
        time=$(extract_time 12 $1 $1)
        echo $time
        printf "$time\n" >> $file_out
    done
}

#! //*                   == Main ==                *//

# time_model 50 5
# time_model 100 3
time_model $1 $2
echo "Result of mpi_s_cmd -> " $(mpi_s_cmd )"