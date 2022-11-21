ORDERS=(0, 5, 10, 20, 50, 100, 150, 200, 300)
SMALL_ORDERS=(1, 2, 3, 4, 5, 10, 20, 30, 40, 50)
SUP_LMAX=500

for o in ${ORDERS[@]}; do
    ./new_model_omp2 --med --lmodel $o --lmax $SUP_LMAX --txt --predict --diff
done

for o in ${SMALL_ORDERS[@]}; do
    ./new_model_omp2 --med --lmodel $o --lmax $SUP_LMAX --txt --predict --diff
done