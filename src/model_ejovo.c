

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

#include "geodesy.h"

int lmax = -1;
int npoint = -1;
char * data_filename;
char * model_filename;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--data FILENAME              input file containing experimental data points\n");
        // printf("--model FILENAME             output file containing the model\n");
        printf("--npoint N                   number of points to read\n");
        printf("--lmax N                     order of the model\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[] = {
                {"data", required_argument, NULL, 'd'},
                {"npoint", required_argument, NULL, 'n'},
                {"lmax", required_argument, NULL, 'l'},
                // {"model", required_argument, NULL, 'm'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'd':
                        data_filename = optarg;
                        break;
                case 'm':
                        model_filename = optarg;
                        break;
                case 'n':
                        npoint = atoi(optarg);
                        break;
                case 'l':
                        lmax = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        // if (data_filename == NULL || model_filename == NULL || lmax < 0 || npoint <= 0)
        if (data_filename == NULL || lmax < 0 || npoint <= 0)
                usage(argv);
}

int main(int argc, char *argv[]) {

    process_command_line_options(argc, argv);

    // initialize rng
    ejovo_seed();

    /**========================================================================
     *!                           DATA_FILES.csv
     *========================================================================**/
    // const char filename_u[] = "ETOPO1_unit.csv";
    // const char filename_t[] = "ETOPO1_tiny.csv";
    // const char filename_s[] = "ETOPO1_small.csv";
    // const char filename_m[] = "ETOPO1_med.csv";

    // const int t_u = 6,     p_u = 12;
    // const int t_t = 18,    p_t = 36;
    // const int t_s = 180,   p_s = 360;
    // const int t_m = 540,   p_m = 1080;

    // // get lmax from the argument pasaed
    // int lmax = 0;

    // if (argc >= 2) {
    //     sscanf(argv[1], "%d", &lmax); 
    //     printf("lmax set to %d with %d args\n", lmax, argc);
    // } else {
    //     fprintf(stderr, "Please enter a value for lmax\n ");
    //     fprintf(stderr, "ex:  %s 10\n", argv[0]);
    //     exit(2);
    // }

    // iso_model *iso = compute_model(5, filename_s, 4, 5);
    // iso_model *iso = compute_model(lmax, filename_u, t_u, p_u);
    // iso_model *iso = compute_model(lmax, "ETOPO1_unif.csv", t_u, p_u);
    // iso_model *iso = compute_model(lmax, filename_t, t_t, p_t);
    // iso_model *iso = compute_model(lmax, filename_t, t_t, p_t);
    // iso_model *iso = compute_model(lmax, filename_s, t_s, p_s);
    // iso_model *iso = compute_model(lmax, filename_m, t_m, p_m);

    // iso_model *iso = compute_model(lmax, data_filename, npoint);
    iso_model *iso = compute_model_binary(lmax, data_filename, npoint, true);
    writeModifiedModel(iso->model, iso->data, "");
    // modelWrit


    write_iso(iso->data, "iso.csv");

    printf("============== Predictions ===============\n");
    update_prediction(iso);
    // print_predictions(iso);
    // print_coeff(iso);

    printf("MSE(iso) : %lf\n", compute_mse(iso));
    // printf("err(iso) : %lf\n", compute_average_error(iso)); 
    printf("ran with lmax := %d\n", lmax);

    /**========================================================================
     *!                           Gradient testing
     *========================================================================**/
    FILE *prediction_file = fopen("prediction.csv", "w");
    // update_prediction(iso);

    for (int i = 0; i < iso->data->N; i++) {

        int i_th = data_i_th(iso->data, i);
        int i_ph = data_i_ph(iso->data, i);
        
        fprintf(prediction_file, "%lf\t%lf\t%lf\n", phiw_to_lamp(vecat_d(iso->data->ph, i_ph)),
                                                    thew_to_phip(vecat_d(iso->data->th, i_th)),
                                                    vecat_d(iso->f_hat, i));

    }

    fclose(prediction_file);

    // int N_runs = 250;
    // double alpha = .01;
    // int N_draws = 50;

    // iso_model *iso_clone1 = copy_iso(iso);
    // iso_model *iso_clone2 = copy_iso(iso);

    // // iso_trained *trained = train_model(iso, N_runs, alpha);
    // // iso_trained *trained_stoch = train_model_stochastic(iso_clone1, N_runs, N_draws, alpha * 10);
    // iso_trained *trained = train_model(iso_clone2, N_runs, alpha);
    
    // report_training(trained_stoch, "test");
    // report_training(trained, "test");


    // printf("Initial error: \n");
    // print_el_d(Matrix_first_d(MSE));

    // printf("\nFinal error\n");
    // print_el_d(Matrix_last_d(MSE));
    // printf("\n");

    // double change_mse = Matrix_last_d(MSE) - Matrix_first_d(MSE);
    // printf("Total change in MSE after training: ");
    // if (change_mse > 0) {
    //     color_red();
    //     printf("%lf\n", change_mse);
    //     color_reset();
    // } else {
    //     printf("%lf\n", change_mse);
    // }

    // printf("rms: %lf\n", sqrt(Matrix_last_d(MSE)));

    // writeModel(iso->model, iso->data, "trained");

    // f_hat = modelPredictN(iso, iso->data->N);
    // reshape_d(f_hat, iso->data->p, iso->data->t);



    // // Matrix_print_d(f_hat);

    // // Print the difference between f_hat and the observed data point
    // Matrix_d *diff = Matrix_subtract_d(f_hat, iso->data->r);

    // // Matrix_print_d(diff);

    // printf("Max diff: %lf\n", max_d(diff));
    // printf("Min diff: %lf\n", min_d(diff));
    // printf("Max_abs diff: %lf\n", maxabs_d(diff));
    // printf("std diff: %lf\n", std_d(diff));
    // printf("total_error: %lf\n", sumabs_d(diff));
    // printf("Alternative mse: %lf\n", mean_squared_d(diff));
    // printf("Meanabs: %lf\n", sumabs_d(diff) / Matrix_size_d(diff));
    // // printf("Average ? error: %lf\n", mean_d(diff));
    // // printf("model dt and dp: %lf, %lf\n", iso->data->dt, iso->data->dp);

    // // printf("mse(manual) : %lf\n", compute_mse(iso));


    // // compute_model(5, filename_s, 5, 10);
}




