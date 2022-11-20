ORDERS=(0, 5, 10, 20, 50, 100, 150, 200, 300)
SUP_LMAX=500

for o in ${ORDERS[@]}; do
    ./new_model_omp2 --hi --lmodel $o --lmax $SUP_LMAX --txt --predict --diff
done