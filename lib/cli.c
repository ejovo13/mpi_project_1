#include "cli.h"

args_t *new_args(int argc, char **argv) {

    args_t *arg = (args_t *) malloc(sizeof(*arg));

    arg->argc = argc;
    arg->argv = argv;
    arg->lmodel = -1;
    arg->lbin = -1;
    arg->size_dataset = NULL;
    arg->data = NULL;
    arg->from = false;
    arg->a = 0;
    arg->txt = false;
    arg->predict = false;
    arg->diff = false;
    arg->recompute = false;
    arg->print_args = false;
    arg->help = false;
    arg->ascii = true;

    arg->rank = 0;
    arg->world_size = 1;

    return arg;
}

void usage(char ** argv)
{
    printf("%s [OPTIONS]\n\n", argv[0]);
    printf("Options:\n\n");
    printf("--[small | med | hi | ultra]    dataset to model\n");
    printf("--lmodel l                      degree of the model to compute\n");
    printf("[--lbin L]                      degree of the stored binary file\n");
    printf("[--from a]                      compute a model of degree lmodel starting from degree a\n");
    printf("                                the difference between the model and predicted value\n");
    printf("[--txt]                         output the model as a text file\n");
    printf("[--predict]                     predict the altitude values (f_hat)\n");
    printf("[--diff]                        predict the altitude values (f_hat) and compute\n");
    printf("[--recompute]                   compute the model coefficients no matter what files\n");
    printf("                                are present\n");
    printf("[--args]                        print the args_t object storing this runs arguments\n");
    printf("[--help]                        print this message\n");
    printf("[--no-ascii]                    dont print the world ascii art\n");
    printf("\n");
    printf("\n");
    printf("Examples:\n");
    // printf("");
    printf("\n");
    printf("  %s --small --lmodel 100\n", argv[0]);
    printf("  %s --med --lmodel 300 --lbin 1000\n", argv[0]);
    printf("  %s --med --lmodel 500 --lbin 1000 --from 300\n", argv[0]);
    printf("  %s --hi --lmodel 20 --txt --predict --diff\n", argv[0]);
    // exit(0);
}

void print_args(const args_t *args) {
    printf("args = {\n");
    printf("  .size_dataset\t  = %s\n", args->size_dataset);
    printf("  .lmodel\t  = %d\n", args->lmodel);
    printf("  .lbin\t\t  = %d\n", args->lmodel);
    printf("  .from\t\t  = %s\n", btos(args->from));
    printf("  .a\t\t  = %d\n", args->a);
    printf("  .txt\t\t  = %s\n", btos(args->txt));
    printf("  .predict\t  = %s\n", btos(args->predict));
    printf("  .diff\t\t  = %s\n", btos(args->diff));
    printf("  .recompute\t  = %s\n", btos(args->recompute));
    printf("  .plm_bin\t  = %s\n", args->plm_bin);
    printf("  .coeff_file_bin = %s\n", args->coeff_file_bin);
    printf("  .print_args\t  = %s\n", btos(args->print_args));
    printf("  .ascii\t  = %s\n", btos(args->ascii));
    printf("}\n");

}

args_t *process_command_line_options(int argc, char ** argv, bool __log, int rank)
{
    struct option longopts[] = {
        {"lmodel", required_argument, NULL, 'l'},
        {"lbin", required_argument, NULL, 'L'},
        {"small", no_argument, NULL, 's'},
        {"med", no_argument, NULL, 'm'},
        {"hi", no_argument, NULL, 'h'},
        {"ultra", no_argument, NULL, 'u'},
        {"txt", no_argument, NULL, 't'},
        {"predict", no_argument, NULL, 'p'},
        {"diff", no_argument, NULL, 'd'},
        {"from", required_argument, NULL, 'f'},
        {"recompute", no_argument, NULL, 'r'},
        {"args", no_argument, NULL, 'a'},
        {"no-ascii", no_argument, NULL, 'c'},
        {"help", no_argument, NULL, 'H'},
        {NULL, 0, NULL, 0}
    };

    args_t *args = new_args(argc, argv);

    char ch;
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {

        case 'L':
            args->lbin = atoll(optarg);
            break;
        case 'l':
            args->lmodel = atoll(optarg);
            args->lbin = args->lmodel;
            break;
        case 's':
            args->size_dataset = "small";
            args->data = get_data_small(__log);
            break;
        case 'm':
            args->size_dataset = "med";
            args->data = get_data_med(__log);
            break;
        case 'h':
            args->size_dataset = "hi";
            args->data = get_data_hi(__log);
            break;
        case 'H':
            args->help = true;
            break;
        case 'u':
            args->size_dataset = "ultra";
            args->data = get_data_ultra(__log);
            break;
        case 't':
            args->txt = true;
            break;
        case 'p':
            args->predict = true;
            break;
        case 'd':
            args->diff = true;
            args->predict = true;
            break;
        case 'f':
            args->from = true;
            args->a = atoi(optarg);
            break;
        case 'r':
            args->recompute = true;
            break;
        case 'a':
            args->print_args = true;
            break;
        case 'c':
            args->ascii = false;
            break;
        default:
            errx(1, "Unknown option\n");
        }
    }

    /* missing required args? */
    if (args->size_dataset == NULL || args->data == NULL || args->lbin < 0 || args->lmodel < 0) {
        if (rank == 0) {
            usage(argv);
        }

        free_data_iso(args->data);

        MPI_Finalize();
        exit(0);
    }

    char *plm_bin = calloc(100, sizeof(*plm_bin));
    char *coeff_file_bin = calloc(100, sizeof(*coeff_file_bin));

    sprintf(coeff_file_bin, "sph_%s_%d.bin", args->size_dataset, args->lmodel);
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", args->size_dataset, args->lbin);
    args->plm_bin = plm_bin;
    args->coeff_file_bin = coeff_file_bin;


    return args;
}
