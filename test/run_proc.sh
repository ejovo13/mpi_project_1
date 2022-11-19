#!/usr/bin/bash

# time and run this command 5 times when given the number of processes to launch
echo $0

mpirun -np $1 ./test/test_omp_proc_available